#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 8  // 假设矩阵大小为 N*N

void calculate_block(double *A, double *B, int rows, int cols) {
    for (int i = 1; i < rows - 1; i++) {
        for (int j = 1; j < cols - 1; j++) {
            B[i * cols + j] = (A[(i - 1) * cols + j] + A[i * cols + (j + 1)] +
                               A[(i + 1) * cols + j] + A[i * cols + (j - 1)]) / 4.0;
        }
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rows_per_proc = N / size;
    int cols = N;

    double *A = (double *)malloc((rows_per_proc + 2) * cols * sizeof(double)); // 包含边界
    double *B = (double *)malloc((rows_per_proc + 2) * cols * sizeof(double));

    // 初始化数据
    for (int i = 0; i < (rows_per_proc + 2) * cols; i++) {
        A[i] = rank + 1;  // 示例数据
    }

    // 边界交换
    if (rank > 0) {
        MPI_Send(A + cols, cols, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(A, cols, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank < size - 1) {
        MPI_Send(A + rows_per_proc * cols, cols, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(A + (rows_per_proc + 1) * cols, cols, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // 局部计算
    calculate_block(A, B, rows_per_proc + 2, cols);

    // 汇总数据
    double *result = NULL;
    if (rank == 0) {
        result = (double *)malloc(N * N * sizeof(double));
    }
    MPI_Gather(B + cols, rows_per_proc * cols, MPI_DOUBLE, result, rows_per_proc * cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // 打印结果
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%6.2f ", result[i * N + j]);
            }
            printf("\n");
        }
        free(result);
    }

    free(A);
    free(B);
    MPI_Finalize();
    return 0;
}