### Task2
For compilation in ```/usr/local/Cellar/open-mpi/3.0.0/share/openmpi/mpic++-wrapper-data.txt```
change compiler from clang++ to g++-7.

Run ```mpic++ -fopenmp main.cpp -o main.out``` to compile.

If "mpi.h" not fount use ```export PATH=$PATH:/usr/local/Cellar/open-mpi/3.0.0```.

Use ```mpiexec -n 4 ./main.out``` to run program on n processes.
