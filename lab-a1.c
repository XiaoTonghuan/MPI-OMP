#include <string.h>
#include <stdio.h>
#include <mpi/mpi.h>

#define MAX_SIZE 256

int main(void)
{
    MPI_Init(NULL, NULL);
    int rank, size, len;
    char hostName[MPI_MAX_PROCESSOR_NAME];
    
    MPI_Get_processor_name(hostName, &len);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Number of task=%d, My rank=%d Running on %s.\n", size, rank, hostName);

    // all gather 所有进程的名称
    char all_hostnames[size][MAX_SIZE];
    MPI_Allgather(hostName, len, MPI_CHAR,
                  all_hostnames, len, MPI_CHAR, MPI_COMM_WORLD);

    //确定分组
    int color = size;
    for(int i = 0; i < size; ++i) {
        if(strcmp(all_hostnames[i], hostName) == 0){
            color = i;
            break;
        }
    }
    //rank 如果相同，按照全局 rank 大小排序确定新的 rank
    MPI_Comm node_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &node_comm);

    int subrank = 0;
    MPI_Comm_rank(node_comm, &subrank);

    int sub_size = 0;
    MPI_Comm_size(node_comm, &sub_size);

    printf("groupd by node, Rank %d belongs to node %s, Its new rank is %d, size of subComm is %d\n",rank,hostName,subrank,sub_size);

    MPI_Finalize();
    return 0;
}