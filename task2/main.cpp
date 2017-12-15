#include <cmath>
#include <mpi.h>

#include "block.h"
//#include <omp.h>

// OMP flag
#define OMP_enabled False

// Hyperparameters
#define Length 1.0, 1.0, 1.0
#define Time 0.1
#define TimeSteps 200


int main(int argc, char** argv) {

    int size = 0, rank = 0, print_rank = -1;
    // Init MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get args
    int grid_size = 0;
    if (argc > 1) {
        grid_size = atoi(argv[1]);
        if (argc > 2) print_rank = atoi(argv[2]); }
    if (argc == 1 or grid_size == 0) {
        if (rank == 0) std::cout << "Invalid usage. Use: ./prog <int grid_size> [< int print_rank >]\n";
        MPI_Finalize();
        return 0; }

    // Count grid parameters
    dbl_triplet d = dbl_triplet {Length} / grid_size - 1;
    double tau = Time / TimeSteps;

    // Create block object
    Block block(size, rank, grid_size);

    if (rank == print_rank) {
        std::cout << "Process rank = " << rank << std::endl;
        block.print();
        printf("_____________________\n");

        Array3D tst {int_triplet {3, 3, 3}};
        tst.print();
    }

    MPI_Finalize();
    return 0;
}
