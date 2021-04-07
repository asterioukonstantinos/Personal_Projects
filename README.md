Implementatios of the following projects :

1) Heat Dissipation on Cylinder surface. The application is written in C and the OpenMP standard is used to achieve data-parallelism. The number of threads is specified with the enviroment variable OMP_NUM_THREADS.

2) Heat Dissipation on Cylinder surface. The application is written in C and POSIX threads are used to achieve parallelism. The number of threads is specified through the command line via the -p flag.

3) Historgram of Grayscale image using Pthreads. The application is written in C and each threads gets a unique address space. When all threads join The master thread, the master thread traverses through all the address spaces and accumulates the results. This approach is more efficient compared to atomic operations, mutex locks and semaphores.

4) Mergesort Algorithm. The applcation is written in C and the OpenMP standard is used to achieve task-parallelism. The number of threads is specified with the enviroment variable OMP_NUM_THREADS. This approach is faster compared to using OpenMP sections.

5) N-Body Simulation. The application is written in C and the OpenMPI standard is used. The application run on the DAS-5 supercomputer. Each node has 32  Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz. A block distribution of the workload is performed across the nodes.
