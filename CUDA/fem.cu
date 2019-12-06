#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>

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

__global__ void Solve(float *dA, float *dF, float *dX0, float *dX1, long long unsigned int N)
{
    float aa, sum = 0.;
    long long unsigned int t = blockIdx.x * blockDim.x + threadIdx.x;

    for (long long unsigned int j = 0; j < N; j++)
    {
        sum += dA[j + t * N] * dX0[j];
    }
    aa = dA[t + t * N];
    dX1[t] = dX0[t] + (dF[t] - sum) / aa;
}

__global__ void Eps(float *dX0, float *dX1, float *delta, long long unsigned int N)
{
    long long unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    delta[i] = fabs(dX0[i] - dX1[i]);
    dX0[i] = dX1[i];
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
    cudaEvent_t GPUstart, GPUstop;
    float GPUtime = 0.0f;

    N_threads = 6;

    if ((nodes_num % N_threads) == 0)
    {
        N_blocks = (nodes_num / N_threads);
    }
    else
    {
        N_blocks = (nodes_num / N_threads) + 1;
    }

    printf("x_nodes_num = %llu\ny_nodes_num=%llu\nN_threads=%llu\nN_blocks=%llu\n", x_nodes_num, y_nodes_num, N_threads, N_blocks);

    dim3 Threads(N_threads);
    dim3 Blocks(N_blocks);

    cudaEventCreate(&GPUstart);
	cudaEventCreate(&GPUstop);

	cudaEventRecord(GPUstart, 0);

    float *hA = (float *)calloc(nodes_num * nodes_num, sizeof(float));
    float *hF = (float *)calloc(nodes_num, sizeof(float));
    float *hX = (float *)calloc(nodes_num, sizeof(float));
    float *hX0 = (float *)calloc(nodes_num, sizeof(float));
    float *hX1 = (float *)calloc(nodes_num, sizeof(float));
    float *hDelta = (float *)calloc(nodes_num, sizeof(float));

    float *dA, *dF, *dX0, *dX1, *dDelta;

    cudaMalloc((void **)&dA, nodes_num * nodes_num * sizeof(float)); // матрица A
    cudaMalloc((void **)&dF, nodes_num * sizeof(float));             // правая часть F
    cudaMalloc((void **)&dX0, nodes_num * sizeof(float));            // решение X(n)
    cudaMalloc((void **)&dX1, nodes_num * sizeof(float));            // решение X(n+1)
    cudaMalloc((void **)&dDelta, nodes_num * sizeof(float));         // разница |X(n+1)- X(n)|

    fp = fopen("result_CUDA.txt", "w");
    if (fp == NULL)
    {
        printf("Open failed\n");
        return -1;
    }

    long long unsigned int  times = (TIME - 1) / delta;
    printf("Times: %llu\n", times);

    initX(hX);
    initAB(hA, hF, hX);
    initX0(hA, hX0, hF);
    cudaMemcpy(dA, hA, nodes_num * nodes_num * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dX0, hX0, nodes_num * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dF, hF, nodes_num * sizeof(float), cudaMemcpyHostToDevice);

    for (i = 1; i <= times; i++)
    {
        fprintf(fp, "%llu sec\n", i);
        printf("%llu\n", i);
        long long unsigned int k = 0;
        eps = INT_MAX;
        while (eps > EPS)
        {
            k++;
            Solve<<<Blocks, Threads>>>(dA, dF, dX0, dX1, nodes_num);
            Eps<<<Blocks, Threads>>>(dX0, dX1, dDelta, nodes_num);

            cudaMemcpy(hX1, dX1, nodes_num * sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(hDelta, dDelta, nodes_num * sizeof(float), cudaMemcpyDeviceToHost);
            eps = 0.;
            for (j = 0; j < nodes_num; j++)
            {
                eps += hDelta[j];
            }
            eps = eps / nodes_num;
            // printf("\n Eps[%i]=%e ", k, eps);
            // cudaMemcpy(dX0, dX1, nodes_num * sizeof(float), cudaMemcpyDeviceToDevice);
        }

        printToFile(hX1);
        cudaMemcpy(hX, dX1, nodes_num * sizeof(float), cudaMemcpyDeviceToHost);
        initB(hF, hX);
        cudaMemcpy(dF, hF, nodes_num * sizeof(float), cudaMemcpyHostToDevice);
        initX0(hA, hX0, hF);
        cudaMemcpy(dX0, hX0, nodes_num * sizeof(float), cudaMemcpyHostToDevice);
    }

    cudaEventRecord(GPUstop, 0);
	cudaEventSynchronize(GPUstop);

	cudaEventElapsedTime(&GPUtime, GPUstart, GPUstop);
	printf("GPU time : %.3f ms\n", GPUtime);


    fclose(fp);

    // Оѝвобождение памѝи
    free(hA);
    free(hF);
    free(hX);
    free(hX0);
    free(hX1);
    free(hDelta);

    cudaFree(dA);
    cudaFree(dF);
    cudaFree(dX0);
    cudaFree(dX1);
    cudaFree(dDelta);
}
