#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>

#define N 4  // 总进程数
#define P 2   // 参数服务器进程数
#define Q (N - P)  // 工作进程数

// 每个工作进程和参数服务器进程的互动过程
void worker(int rank, int p, MPI_Comm comm) {
    int server_rank = rank % p;  // 每个工作进程发送给哪个参数服务器
    double rand_num = rand() % 100;  // 生成一个随机数

    // 发送给对应的参数服务器
    MPI_Send(&rand_num, 1, MPI_DOUBLE, server_rank, 0, comm);

    // 接收更新后的参数（平均值）
    double updated_value;
    MPI_Recv(&updated_value, 1, MPI_DOUBLE, server_rank, 1, comm, MPI_STATUS_IGNORE);
    
    printf("Worker %d sent %f to Server %d and received updated value %f\n", rank, rand_num, server_rank, updated_value);
}

void parameter_server(int rank, int q, MPI_Comm comm) {
    double values[q];  // 存储从工作进程接收到的数值

    // 接收来自工作进程的数据
    for (int i = 0; i < q; i++) {
        MPI_Recv(&values[i], 1, MPI_DOUBLE, P + q * rank + i, 0, comm, MPI_STATUS_IGNORE);
    }

    // 计算平均值
    double sum = 0;

    for (int i = 0; i < q; i++) {
        sum += values[i];
    }
    double avg = sum / q;

    // 将平均值返回给对应的工作进程
    for (int i = 0; i < q; i++) {
        MPI_Send(&avg, 1, MPI_DOUBLE, P + rank * q + i, 1, comm);
    }

    printf("Parameter Server %d computed average: %f and sent to workers\n", rank, avg);
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != N) {
        printf("The number of processes must be %d.\n", N);
        MPI_Finalize();
        return 1;
    }

    if (rank < P) {
        // 参数服务器进程
        parameter_server(rank, Q / P, MPI_COMM_WORLD);
    } else {
        // 工作进程
        worker(rank, P, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
