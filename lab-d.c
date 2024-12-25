#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 128 //矩阵A, B的维度
#define P 4
#define BLOCK_SIZE (N / P)  // 每个子矩阵的大小

void Print_Mat(int *A, int n){//打印矩阵
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%d  ",A[i * n + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int int_sqrt(int p){//计算p的根号，返回整数结果
    for(int i = 1; i <= p; i++){
        if(i * i == p){
            return i;
        }
    }
    return -1;
}

//获得A的第i,j个n*n方块,储存在a中
void Get_Block(int *A, int *a, int i, int j, int n){
    for(int k = 0; k < n; k++){
        for(int l = 0; l < n; l++){
            a[k * n + l] = A[i * n * N + j * n + k * N + l];
        }
    }
}

//矩阵相乘的串行算法，计算矩阵乘法 a*b，计算结果储存在c中
void Multi_Mat(int *a, int *b, int *c, int n){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
                c[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
}

//fox矩阵乘法
void Fox(int *a, int *b, int *c, int p, int n, int rank){
    int *temp_a = malloc(n * n * sizeof(int));//用来接收a
    int *temp_b = malloc(n * n * sizeof(int));//用来接收b

    //将处理器记为P_ij
    int j = rank % p;
    int i = rank  / p;

    //对方块 c 初始化
    for(int k = 0; k < n * n; k++){
        c[k] = 0;
    }
    int senddest = 0, recvdest = 0;
    // fox循环
    for(int round = 0; round < p; round++){

        int ir = (i + round) % p;

        if(j == ir){ //
            for(int k = 0; k < p; k++) {// Transfer
                if(k != j){ //向这一行所有不是自己的处理器转发数据
                    MPI_Send(a, n * n, MPI_INT, i * p + k, (round + 1) * p, MPI_COMM_WORLD);
                }
                else {
                    for(int l = 0; l < n * n; l++){
                        temp_a[l] = a[l];
                    }
                }
            }
        }
        else {
            recvdest = i * p + ir;
            MPI_Recv(temp_a, n * n, MPI_INT, recvdest, (round + 1) * p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        //做矩阵乘法
        Multi_Mat(temp_a, b, c, n);

        if(round + 1 == p) break;
        
        //将块B_ij上移
        senddest = i > 0 ? i - 1:i - 1 + p;
        recvdest = i < p - 1 ? i + 1:i + 1 - p;
        MPI_Sendrecv(b, n * n, MPI_INT, senddest * p + j, 0, 
                    temp_b, n * n, MPI_INT, recvdest * p + j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(int i = 0; i < n * n; i++){
            b[i] = temp_b[i];
        }
    }
}

//将C的第i,j个n*n方块用c赋值
void Copy_Block(int *c, int *C, int i, int j, int n){
    for(int k = 0; k < n; k++){
        for(int l = 0; l < n; l++){
            C[i * n * N + j * n + k * N + l] = c[k * n + l];
        }
    }
}

// processor0 接收其他处理器的计算结果
void Recv_Block(int *c, int *C, int sp, int n){
    for(int i = 0; i < sp; i++){
        for(int j = 0; j < sp; j++){
            if(i + j != 0) 
                MPI_Recv(c, n * n, MPI_INT, i * sp + j, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            Copy_Block(c, C, i, j, n);
        }
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *A, *B, *C;
    if(rank == 0) { 
        //随机生成一个矩阵
        srand(time(NULL));
        A = (int *)malloc(N * N * sizeof(int));
        B = (int *)malloc(N * N * sizeof(int));
        C = (int *)malloc(N * N * sizeof(int));
        for(int i = 0; i < N * N; i++){
            A[i] = rand() % 10;
            B[i] = rand() % 10;
        }
    }

    int p = int_sqrt(P);
    
    int *a = malloc(BLOCK_SIZE * BLOCK_SIZE * sizeof(int));
    int *b = malloc(BLOCK_SIZE * BLOCK_SIZE * sizeof(int));
    int *c = malloc(BLOCK_SIZE * BLOCK_SIZE * sizeof(int));

    //计时开始
    MPI_Barrier(MPI_COMM_WORLD);
    double start_fox = MPI_Wtime();

    if(rank == 0){
        //由 rank0 把A，B的 BLOCK * BLOCK 方块发送给另外p * p - 1 个处理器
        for(int i = 0; i < p; i++){
            for(int j = 0; j < p; j++){
                if(i == 0 && j == 0) continue;
                Get_Block(A, a, i, j, BLOCK_SIZE);
                Get_Block(B, b, i, j, BLOCK_SIZE);
                //Print_Mat(a, n);
                MPI_Send(a, BLOCK_SIZE * BLOCK_SIZE, MPI_INT, i * p + j, 'a', MPI_COMM_WORLD);
                MPI_Send(b, BLOCK_SIZE * BLOCK_SIZE, MPI_INT, i * p + j, 'b', MPI_COMM_WORLD);
            }
        }
        // processor0 储存方块A_00, B_00
        Get_Block(A, a, 0, 0, BLOCK_SIZE);
        Get_Block(B, b, 0, 0, BLOCK_SIZE);
    }
    else{// 其他处理器接收方块A_ij, B_ij
        MPI_Recv(a, BLOCK_SIZE * BLOCK_SIZE, MPI_INT, 0, 'a', MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, BLOCK_SIZE * BLOCK_SIZE, MPI_INT, 0, 'b', MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    // Fox矩阵乘法
    Fox(a, b, c, p, BLOCK_SIZE, rank);

    if(rank != 0){//其他处理器向 processor0 发送计算结果
        MPI_Send(c, BLOCK_SIZE * BLOCK_SIZE, MPI_INT, 0, 3, MPI_COMM_WORLD);
    }
    else{//processor0 接收方块结果，并赋值给矩阵C
        Recv_Block(c, C, p, BLOCK_SIZE);
    }

    //计时结束
    double finish_fox = MPI_Wtime();

    // 由 processor0 打印计算结果
    if (rank == 0)
    {   
        // 检验fox乘法结果正确性
        #if 0
        printf("Matrix A = \n");Print_Mat(A, N);
        printf("Matrix B = \n");Print_Mat(B, N);
        printf("fox计算结果 C = \n");Print_Mat(C, N);
        #endif
        printf("fox并行乘法耗时: %f\n", finish_fox - start_fox);
        
        //归零
        for(int i=0; i < N * N; i++){
            C[i]=0;
        }

        //串行乘法
        double start =  MPI_Wtime();
        //Multi_Mat(A, B, C, N);
        double finish = MPI_Wtime();

        #if 0
        printf("串行计算结果 C = \n");Print_Mat(C, N);
        #endif
        printf("串行乘法耗时: %f\n", finish - start);

        //计算加速比
        double rate = (finish - start)/(finish_fox - start_fox);
        printf("并行加速比为: %f\n", rate);
    }

    MPI_Finalize();
    return 0;
}

