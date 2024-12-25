#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include <math.h>

#define N 1024  // 矩阵的维度
#define P 4     // 进程数，假设为 4 (2x2 网格)
#define BLOCK_SIZE (N / P)  // 每个进程的子矩阵块大小

void matmul_block(double* A, double* B, double* C, int block_size) {
    // 矩阵 A、B、C 的矩阵块乘法
    for (int i = 0; i < block_size; i++) {
        for (int j = 0; j < block_size; j++) {
            for (int k = 0; k < block_size; k++) {
                C[i * block_size + j] += A[i * block_size + k] * B[k * block_size + j];
            }
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 假设进程总数为 P*P
    int sqrt_p = sqrt(size);  // 进程网格的边长

    // 生成本地的 A、B 和 C 矩阵块
    double* A_local = (double*)malloc(sizeof(double) * BLOCK_SIZE * BLOCK_SIZE);
    double* B_local = (double*)malloc(sizeof(double) * BLOCK_SIZE * BLOCK_SIZE);
    double* C_local = (double*)malloc(sizeof(double) * BLOCK_SIZE * BLOCK_SIZE);

    // 初始化 A_local 和 B_local
    // 在实际应用中，可能需要从文件或其他进程获取数据
    for (int i = 0; i < BLOCK_SIZE; i++) {
        for (int j = 0; j < BLOCK_SIZE; j++) {
            A_local[i * BLOCK_SIZE + j] = rank + i + j;
            B_local[i * BLOCK_SIZE + j] = rank + i - j;
            C_local[i * BLOCK_SIZE + j] = 0.0;
        }
    }

    // FOX 算法
    for (int phase = 0; phase < sqrt_p; phase++) {
        int source = (rank + phase) % sqrt_p;
        int dest = (rank + sqrt_p - phase) % sqrt_p;

        // 每一轮，进程交换数据
        double* temp_A = (double*)malloc(sizeof(double) * BLOCK_SIZE * BLOCK_SIZE);
        double* temp_B = (double*)malloc(sizeof(double) * BLOCK_SIZE * BLOCK_SIZE);
        
        // 数据交换 (A 与 B)
        MPI_Sendrecv(A_local, BLOCK_SIZE * BLOCK_SIZE, MPI_DOUBLE, dest, 0, 
                     temp_A, BLOCK_SIZE * BLOCK_SIZE, MPI_DOUBLE, source, 0, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        MPI_Sendrecv(B_local, BLOCK_SIZE * BLOCK_SIZE, MPI_DOUBLE, dest, 1, 
                     temp_B, BLOCK_SIZE * BLOCK_SIZE, MPI_DOUBLE, source, 1, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // 计算矩阵块乘法
        matmul_block(A_local, temp_B, C_local, BLOCK_SIZE);

        // 清理
        free(temp_A);
        free(temp_B);

        // 在每轮结束后，更新 A_local 和 B_local
        for (int i = 0; i < BLOCK_SIZE; i++) {
            for (int j = 0; j < BLOCK_SIZE; j++) {
                A_local[i * BLOCK_SIZE + j] = temp_A[i * BLOCK_SIZE + j];
                B_local[i * BLOCK_SIZE + j] = temp_B[i * BLOCK_SIZE + j];
            }
        }
    }

    // 收集结果
    // 可以通过 MPI_Gather 或其他方式将结果合并到根进程
    if (rank == 0) {
        // 输出最终结果
        printf("Matrix multiplication completed.\n");
    }

    // 清理资源
    free(A_local);
    free(B_local);
    free(C_local);

    MPI_Finalize();
    return 0;
}
