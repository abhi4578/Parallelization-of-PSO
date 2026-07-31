

# Paralelización-de-PSO
Paralelización de la Optimización de Enjambre de Partículas utilizando los frameworks OpenMp, MPI y CUDA, y comparación de su rendimiento.
PSO serial: https://github.com/m-ahsen/pso/blob/master/pso_serial.c

    comando para compilar omp.c    : gcc -fopenmp omp.c -lm -lgsl -lgslcblas
    comando para compilar mpiomp.c : mpicc -fopenmp mpiomp.c -lm -lgsl -lgslcblas
    comando para compilar mpipso.c : mpicc mpipso.c -lm -lgsl -lgslcblas
