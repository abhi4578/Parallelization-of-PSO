gcc -fopenmp omp.c -lm -lgsl -lgslcblas -o omp
gcc -fopenmp omp_parallel.c -lm -lgsl -lgslcblas -o omp_parallel
mpicc -fopenmp mpiomp.c -lm -lgsl -lgslcblas -o mpiomp
mpicc mpipso.c -lm -lgsl -lgslcblas -o mpi
gcc serialpso.c -lm -lgsl -lgslcblas -o serial