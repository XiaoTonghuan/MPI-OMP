#include <stdio.h>
#include <mpi/mpi.h>

void binary_sum_up(int * sendbuf, int rank, int size) {
    
    rank += 1; //从1开始更好计算
    MPI_Status status;

    int lson = rank * 2;
    int rson = rank * 2 + 1;
    int recvbuf = 0;
    if(lson <= size) {
        MPI_Recv(&recvbuf, 1, MPI_INT , lson - 1, 0, MPI_COMM_WORLD, &status);
        *sendbuf += recvbuf;
    }
    if(rson <= size) {
        MPI_Recv(&recvbuf, 1, MPI_INT, rson - 1, 0, MPI_COMM_WORLD, &status);
        *sendbuf += recvbuf;
    }
    int fa = rank / 2;
    if(fa <= 0) return;
    MPI_Send(sendbuf, 1, MPI_INT, fa - 1, 0, MPI_COMM_WORLD);
    
}

void binary_sum_down(int * sendbuf, int rank, int size) {
    rank += 1;
    MPI_Status status;
    
    if(rank != 1)  {
        int recvbuf = 0;
        MPI_Recv(&recvbuf, 1, MPI_INT, rank / 2 - 1, 0, MPI_COMM_WORLD, &status);
        *sendbuf = recvbuf;
    }
    if(rank <= size / 2) {
        if (rank * 2 - 1 < size) {
            MPI_Send(sendbuf, 1, MPI_INT, rank * 2 - 1, 0, MPI_COMM_WORLD);
        }
        if (rank * 2 < size) {
            MPI_Send(sendbuf, 1, MPI_INT, rank * 2, 0, MPI_COMM_WORLD);
        }    
    }
}

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //设每个进程持有的数据恰好等于rank + 1
    int local_data = rank + 1;
    int res = 0;
    // printf("Process %d local data: %d\n", rank, local_data);

    // 调用蝶式全和
    binary_sum_up(&local_data,rank, size);
    res = local_data;
    binary_sum_down(&res,rank,size);
    // 打印每个进程的总和
    printf("Process %d final sum: %d\n", rank, res);

    MPI_Finalize();
    return 0;
}