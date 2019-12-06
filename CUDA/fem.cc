#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <time.h>

#define TEMP_BOT 0
#define TEMP_LEFT 0
#define TEMP_UP 200
#define TEMP_RIGHT 200
#define TEMP_BEGIN 20

#define ANIMATION_FRAME_DELAY 10
#define TIME 21
#define ALFA1 1
#define ALFA2 1

#define EPS 0.001

float delta;
float eps = INT_MAX;
long long unsigned int x_nodes_num;
long long unsigned int y_nodes_num;
long long unsigned int nodes_num;
FILE *fp;
FILE *gnuPlot;

long long unsigned int N_threads = 512;
long long unsigned int N_blocks;

float get_element(float *x, long long unsigned int  row_size, long long unsigned int  i, long long unsigned int  j)
{
    return x[i * row_size + j];
}

void set_element(float value, float *x, long long unsigned int row_size, long long unsigned int  i, long long unsigned int  j)
{
    x[i * row_size + j] = value;
}

void initAB(float *A, float *B, float *x)
{
    long long unsigned int i = 0, j = 0, offset = 0, step = 0, index = 0;
    // Заполнѝем нижний левый угол плаѝтины
    set_element(1, A, nodes_num, 0, 0);
    set_element(-1, A, nodes_num, 0, 1);
    B[0] = 0;
    // Заполнѝем верхний край плаѝтины
    offset = x_nodes_num * (y_nodes_num - 1);
    for (i = 0; i < x_nodes_num; i++)
    {
        set_element(1, A, nodes_num, i + offset, i + offset);
        B[i + offset] = TEMP_UP;
    }

    // Заполнѝем правый край плаѝтины
    offset = x_nodes_num - 1;
    index = offset;
    step = x_nodes_num;
    for (i = 0; i < y_nodes_num; i++)
    {
        set_element(1, A, nodes_num, index, index);
        B[index] = TEMP_RIGHT;
        index += step;
    }

    // Заполнѝем левый край плаѝтины
    offset = x_nodes_num;
    index = offset;
    step = x_nodes_num;
    for (i = 0; i < y_nodes_num - 2; i++)
    {
        if (i == 2)
        {
            set_element(1 + delta, A, nodes_num, index, index);
            set_element(-1, A, nodes_num, index, index + 1);
            B[index] = 0;
        }
        else
        {
            set_element(1, A, nodes_num, index, index);
            set_element(-1, A, nodes_num, index, index + 1);
            B[index] = 0;
        }
        index += step;
    }

    // Заполнѝем нижний край плаѝтины
    offset = 1;
    index = offset;
    step = 1;
    for (i = 0; i < x_nodes_num - 2; i++)
    {
        set_element(1, A, nodes_num, index, index);
        set_element(-1, A, nodes_num, index, index + x_nodes_num);
        B[index] = 0;
        index += step;
    }
    // Заполнѝем внутреннюю чаѝть плаѝтины
    offset = x_nodes_num + 1;
    step = x_nodes_num;
    for (i = 0; i < y_nodes_num - 2; i++)
    {
        for (j = 0; j < x_nodes_num - 2; j++)
        {
            index = offset + i * x_nodes_num + j;
            set_element(5, A, nodes_num, index, index);
            set_element(-1, A, nodes_num, index, index - 1);
            set_element(-1, A, nodes_num, index, index + 1);
            set_element(-1, A, nodes_num, index, index - x_nodes_num);
            set_element(-1, A, nodes_num, index, index + x_nodes_num);
            B[index] = x[index];
        }
    }
}

void initB(float *B, float *x)
{
    long long unsigned int i = 0, j = 0, offset = 0, index = 0, step = 0;
    // Заполнѝем внутреннюю чаѝть плаѝтины
    offset = x_nodes_num + 1;
    step = x_nodes_num;
    for (i = 0; i < y_nodes_num - 2; i++)
    {
        for (j = 0; j < x_nodes_num - 2; j++)
        {
            index = offset + i * x_nodes_num + j;
            B[index] = x[index];
        }
    }
}

void initX(float *x)
{
    for (long long unsigned int  i = 0; i < nodes_num; i++)
    {
        x[i] = (float)TEMP_BEGIN;
    }
}

void initX0(float *A, float *x, float *B)
{
    for (long long unsigned int  i = 0; i < nodes_num; i++)
    {
        float aii = get_element(A, nodes_num, i, i);
        x[i] = B[i] / aii;
    }
}

/*
void initA(float *A, float *B, float *x) {
    long long unsigned int i = 0, j = 0, offset = 0, step = 0, index = 0;

    // Заполнѝем нижний левый угол плаѝтины
    set_element(1, A, nodes_num, 0, 0);
    set_element(-1, A, nodes_num, 0, 1);

    // Заполнѝем верхний край плаѝтины
    offset = x_nodes_num * (y_nodes_num - 1);
    for (i = 0; i < x_nodes_num; i++)
    {
        set_element(1, A, nodes_num, i + offset, i + offset);
    }
    // Заполнѝем правый край плаѝтины
    offset = x_nodes_num - 1;
    index = offset;
    step = x_nodes_num;
    for (i = 0; i < y_nodes_num; i++)
    {
        set_element(1, A, nodes_num, index, index);
        index += step;
    }
    // Заполнѝем левый край плаѝтины
    offset = x_nodes_num;
    index = offset;
    step = x_nodes_num;
    for (i = 0; i < y_nodes_num - 2; i++)
    {
        if (i == 2)
        {
            set_element(1 + delta, A, nodes_num, index, index);
            set_element(-1, A, nodes_num, index, index + 1);
        }
        else
        {
            set_element(1, A, nodes_num, index, index);
            set_element(-1, A, nodes_num, index, index + 1);
        }
        index += step;
    }
    // Заполнѝем нижний край плаѝтины
    offset = 1;
    index = offset;
    step = 1;
    for (i = 0; i < x_nodes_num - 2; i++)
    {
        set_element(1, A, nodes_num, index, index);
        set_element(-1, A, nodes_num, index, index + x_nodes_num);
        index += step;
    }
    // Заполнѝем внутреннюю чаѝть плаѝтины
    offset = x_nodes_num + 1;
    step = x_nodes_num;
    for (i = 0; i < y_nodes_num - 2; i++)
    {
        for (j = 0; j < x_nodes_num - 2; j++)
        {
            index = offset + i * x_nodes_num + j;
            set_element(5, A, nodes_num, index, index);
            set_element(-1, A, nodes_num, index, index - 1);
            set_element(-1, A, nodes_num, index, index + 1);

            set_element(-1, A, nodes_num, index, index - x_nodes_num);
            set_element(-1, A, nodes_num, index, index + x_nodes_num);
        }
    }
}
*/

// void printToFile(float *temp)
// {
//     for (long long unsigned int i = 0; i < x_nodes_num; i++)
//     {
//         for (long long unsigned int j = 0; j < y_nodes_num; j++)
//         {
//             float xCoor = i * delta;
//             float y_coor = j * delta;
//             fprintf(fp, "%.1lf %.1lf %.1lf\n", xCoor, y_coor, temp[j * x_nodes_num + i]);
//         }
//         fprintf(fp, "\n");
//     }
//     fprintf(fp, "\n");
// }

void printToFile(float *temp)
{
    for (int j = y_nodes_num - 1; j >= 0; j--)
    {
        for (int i = 0; i < x_nodes_num; i++)
        {
            fprintf(fp, "%.1f\t", temp[j * x_nodes_num + i]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

void zeroingA(float *A)
{
    long long unsigned int i = 0, j = 0;
    for (i = 0; i < nodes_num; i++)
    {
        for (j = 0; j < nodes_num; j++)
        {
            set_element(0, A, nodes_num, i, j);
        }
    }
}

void zeroingB(float *vectorB)
{
    long long unsigned int i = 0;
    for (i = 0; i < nodes_num; i++)
    {
        vectorB[i] = 0;
    }
}

void Solve(float *A, float *F, float *X0, float *X1, long long unsigned int N) 
{
    for (long long unsigned int i = 0; i < N; i++) {
        float aa, sum = 0.;
        for (long long unsigned int j = 0; j < N; j++)
        {
            sum += A[j + i * N] * X0[j];
        }
        aa = A[i + i * N];
        X1[i] = X0[i] + (F[i] - sum) / aa;
    }
}

void Eps(float *X0, float *X1, float *delta, long long unsigned int N)
{
    for (long long unsigned int i = 0; i < N; i++) 
    {
        delta[i] = fabs(X0[i] - X1[i]);
        X0[i] = X1[i];
    }
}


int main(int argc, char **argv)
{
    if (argc != 3) {
        printf("Error data. Please, enter these aruments:\n1st argument - x_nodes;\n2nd argument - y_nodes.\n");
        return 1;
    }

    x_nodes_num = atoi(argv[1]);
    y_nodes_num = atoi(argv[2]);

    if (x_nodes_num == 0 || y_nodes_num == 0) {
        printf("Error data. Please, enter 2 integers.\n");
        return 1;
    }

    nodes_num = x_nodes_num * y_nodes_num;
    delta = 1;
    long long unsigned int  i = 0, j = 0;
    float CPUstart, CPUstop;
    float CPUtime = 0.0f;

    printf("x_nodes_num = %llu\ny_nodes_num=%llu\n", x_nodes_num, y_nodes_num);

    CPUstart = clock();

    float *A = (float *)calloc(nodes_num * nodes_num, sizeof(float));
    float *F = (float *)calloc(nodes_num, sizeof(float));
    float *X = (float *)calloc(nodes_num, sizeof(float));
    float *X0 = (float *)calloc(nodes_num, sizeof(float));
    float *X1 = (float *)calloc(nodes_num, sizeof(float));
    float *Delta = (float *)calloc(nodes_num, sizeof(float));

    fp = fopen("result_CPP.txt", "w");
    if (fp == NULL)
    {
        printf("Open failed\n");
        return -1;
    }

    long long unsigned int  times = (TIME - 1) / delta;
    printf("Times: %llu\n", times);

    initX(X);
    initAB(A, F, X);
    initX0(A, X0, F);

    for (i = 1; i <= times; i++)
    {
        fprintf(fp, "%llu sec\n", i);
        printf("%llu\n", i);
        long long unsigned int k = 0;
        eps = INT_MAX;
        while (eps > EPS)
        {
            k++;
            Solve(A, F, X0, X1, nodes_num);
            Eps(X0, X1, Delta, nodes_num);

            eps = 0.;
            for (j = 0; j < nodes_num; j++)
            {
                eps += Delta[j];
            }
            eps = eps / nodes_num;
            // printf("\n Eps[%i]=%e ", k, eps);
            // memcpy(X0, X1, nodes_num * sizeof(float));
        }

        printToFile(X1);
        memcpy(X, X1, nodes_num * sizeof(float));
        initB(F, X);
        initX0(A, X0, F);
    }

    CPUstop = clock();
	CPUtime = 1000.*(CPUstop - CPUstart) / CLOCKS_PER_SEC;
	printf("CPU time : %.3f ms\n", CPUtime);

    fclose(fp);

    // Оѝвобождение памѝи
    free(A);
    free(F);
    free(X);
    free(X0);
    free(X1);
    free(Delta);

}
