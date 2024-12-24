cd ex1 && mpicc -Wall 1-1.c -o ../build/1-1
cd ../build && mpirun -np 4 ./1-1