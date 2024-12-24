#include <stdio.h>
#include <mpi/mpi.h>
int main(void)
{
    MPI_Init(NULL, NULL);
    int rank, size, len;
    char hostName[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(hostName, &len);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Number of task=%d, My rank=%d Running on %s.\n", size, rank, hostName);
    MPI_Finalize();
    return 0;
}