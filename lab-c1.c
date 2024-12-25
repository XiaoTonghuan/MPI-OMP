#include <stdio.h>
#include <mpi/mpi.h>

void butterfly_sum(int *sendbuf,int rank, int size, MPI_Comm comm) {
    
    int step = 1;
    
    while (step < size) {
        // 计算要和哪个进程交换数据
        int partner = rank ^ step;  // XOR operation to get the partner
        if (partner < size) { // 确保 partner 是有效的
            int recv_data = 0;
            MPI_Sendrecv(sendbuf, 1, MPI_INT, partner, 0,
                         &recv_data, 1, MPI_INT, partner, 0,
                         comm, MPI_STATUS_IGNORE);
            *sendbuf += recv_data;
        }
        step <<= 1;  // 翻倍步长
    }
    
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //设每个进程持有的数据恰好等于rank + 1
    int local_data = rank + 1;

    // printf("Process %d local data: %d\n", rank, local_data);

    // 调用蝶式全和
    butterfly_sum(&local_data,rank, size, MPI_COMM_WORLD);

    // 打印每个进程的总和
    printf("Process %d final sum: %d\n", rank, local_data);

    MPI_Finalize();
    return 0;
}
