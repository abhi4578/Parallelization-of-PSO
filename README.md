# Parallelization-of-PSO
Paralleling  Particle Swarm Optimization using OpenMp, MPI and CUDA frameworks and comparing the performance
Serial PSO: https://github.com/m-ahsen/pso/blob/master/pso_serial.c

    command to compile omp.c    : gcc -fopenmp omp.c -lgsl
    command to compile mpiomp.c : mpicc -fopenmp mpiomp.c -lgsl
    command to compile mpipso.c : mpicc mpipso.c -lgsl
