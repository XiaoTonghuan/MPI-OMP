#include <string.h>
#include <stdio.h>
#include <mpi/mpi.h>

#define MAX_SIZE 256
#define MAX_MSG 256

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size, len;
    char hostName[MPI_MAX_PROCESSOR_NAME];

    MPI_Get_processor_name(hostName, &len);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Number of tasks=%d, My rank=%d Running on %s.\n", size, rank, hostName);

    // 收集所有进程的主机名
    char all_hostnames[size][MAX_SIZE];
    MPI_Allgather(hostName, MAX_SIZE, MPI_CHAR, all_hostnames, MAX_SIZE, MPI_CHAR, MPI_COMM_WORLD);

    // 确定分组
    int color = rank;
    for (int i = 0; i < size; ++i) {
        if (strcmp(all_hostnames[i], hostName) == 0) {
            color = i;
            break;
        }
    }

    // 按分组创建子通信器
    MPI_Comm node_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &node_comm);

    // 获取子通信器中的 rank 和 size
    int subrank, subsize;
    MPI_Comm_rank(node_comm, &subrank);
    MPI_Comm_size(node_comm, &subsize);

    printf("Grouped by node: Rank %d belongs to node %s, Its new rank is %d, size of subComm is %d\n", 
           rank, hostName, subrank, subsize);
    // msg
    char msg[MAX_MSG];
    if(rank == 0) {
        snprintf(msg, sizeof(msg), "HelloMPI!");
    }
    MPI_Bcast(msg, sizeof(msg), MPI_CHAR, 0, node_comm);

    printf("node: %s, receive msg %s\n",hostName,msg);
    
    // 释放子通信器
    MPI_Comm_free(&node_comm);
    MPI_Finalize();

    return 0;
}
