# 并行计算实验
## MPI 上机实验

### lab-a

- 写个将 MPI 进程按其所在节点分组的程序
- 在 1.1 的基础上，写个广播程序，主要思想是：按节点分组后，广播的 root 进程将消息
“发送”给各组的“0 号”，再由这些“0”号进程在其小组内执行 MPI_Bcast

### lab-b

使用 MPI_Send 和 MPI_Recv 来模拟 MPI_Alltoall。将你的实验与相关 MPI 通信函数做评测和对比

### lab-c

N 个进程求 N 个数全和, 要求每个进程保持全和

- 使用 MPI 实现碟式全和
- 使用 MPI 实现二叉树方式求全和

### lab-d

使用 MPI 实现矩阵相乘的 MPI 算法

### lab-e

参数服务器系统 MPI 模拟

设系统中总计有 N 个进程，其中 P 个进程作为参数服务器进程，而 Q 个进程作
为工作进程（N = P + Q， 且 0 < P < < Q）。工作进程和服务器进程的互动过程
如下：

1. 第 i 个工作进程首先产生一个随机数，发送给第 i%P 个参数服务器进程。然
后等待并接收它对应的参数服务器进程发送更新后的数值
2. 每个参数服务器进程等待并接收来自它对应的所有工作进程的数据，在此之
后，经通信，使所有的参数服务器获得所有工作进程发送数据的平均值
3. 每个参数服务器发送该平均值给它对应的所有工作进程

试给出上述互动过程的 MPI 程序实现


### lab-f

矩阵 A 和 B 均为 N*N 的双精度数矩阵，有 P 个处理器。针对以下程序片
段，分别采用按行块连续划分以及棋盘式划分方式，给出相应的 MPI 并
行实现。

```cpp
for(i=1; i<N-1; i++)
    for(j=1;j<N-1; j++)
        B[i][j] = (A[i-1][j] + A[i][j+1] + A[i+1][j] + A[i][j-1]) / 4.0
```cpp