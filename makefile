N = 12000
P = 12
T = 2 

all: jacobi-mpi-c jacobiseq-c

jacobi-mpi-c:
	mpicc jacobi-mpi.c -o jacobi-mpi -lm -Wall -fopenmp

jacobiseq-c:
	gcc -fopenmp jacobiseq.c -o jacobiseq

jacobi-mpi-run:
	mpirun -np $(P) -host hal02,hal03,hal04,hal08,hal09,hal10 --oversubscribe ./jacobi-mpi $(N) $(P) $(T)

jacobiseq-run:
	./jacobiseq $(N) 

clean:
	rm -f jacobi-mpi.o jacobiseq.o jacobi-mpi jacobiseq