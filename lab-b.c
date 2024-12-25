#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include "string.h"

void printRecv(const int* sendbuf, int size) {
    for (int i = 0; i < size; i++) {
        printf("%d ",sendbuf[i]);
    }
    printf("\ndone!\n");
}

void custom_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                     void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int *send_data = (int *)sendbuf;
    int *recv_data = (int *)recvbuf;
    MPI_Status status;

    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Send(send_data + i * sendcount, sendcount, sendtype, i, 0, comm);
        } else {
            memcpy(recv_data + i * recvcount, send_data + i * sendcount, sendcount * sizeof(int));
        }
        if (i != rank) {
            MPI_Recv(recv_data + i * recvcount, recvcount, recvtype, i, 0, comm, &status);
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int data_per_proc = 4; // 每个进程发送的数据量
    int *sendbuf = (int *)malloc(size * data_per_proc * sizeof(int));
    int *recvbuf_custom = (int *)malloc(size * data_per_proc * sizeof(int));
    int *recvbuf_alltoall = (int *)malloc(size * data_per_proc * sizeof(int));

    // 初始化发送缓冲区
    for (int i = 0; i < size * data_per_proc; i++) {
        sendbuf[i] = rank * 10 + i;
    }
    
    
    printRecv(sendbuf, size * data_per_proc);


    // 测试 custom_alltoall
    double start_custom = MPI_Wtime();
    custom_alltoall(sendbuf, data_per_proc, MPI_INT, recvbuf_custom, data_per_proc, MPI_INT, MPI_COMM_WORLD);
    double end_custom = MPI_Wtime();

    // 测试 MPI_Alltoall
    double start_alltoall = MPI_Wtime();
    MPI_Alltoall(sendbuf, data_per_proc, MPI_INT, recvbuf_alltoall, data_per_proc, MPI_INT, MPI_COMM_WORLD);
    double end_alltoall = MPI_Wtime();

    printRecv(recvbuf_alltoall,size * data_per_proc);
    printRecv(recvbuf_custom,size * data_per_proc);

    // 验证结果是否一致
    int is_correct = 1;
    for (int i = 0; i < size * data_per_proc; i++) {
        if (recvbuf_custom[i] != recvbuf_alltoall[i]) {
            is_correct = 0;
            break;
        }
    }
    if (rank == 0) {
        printf("Results match: %s\n", is_correct ? "Yes" : "No");
        printf("Custom Alltoall Time: %f seconds\n", end_custom - start_custom);
        printf("MPI_Alltoall Time: %f seconds\n", end_alltoall - start_alltoall);
    }

    free(sendbuf);
    free(recvbuf_custom);
    free(recvbuf_alltoall);

    MPI_Finalize();
    return 0;
}
