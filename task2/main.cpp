#include <mpi.h>
#include <tuple>
#include "block.h"

// Hyperparameters
#define Length 1.0, 1.0, 1.0
#define Time 0.01
#define TimeSteps 20


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
                if (argc > 2) print_rank = atoi(argv[2]);
        }
        if (argc == 1 or grid_size == 0) {
                if (rank == 0) std::cout << "Invalid usage. Use: ./prog <int grid_size> [< int print_rank >]\n";
                MPI_Finalize();
                return 0;
        }

        // Count grid parameters
        triplet<double> L {Length};
        triplet<double> d = L / (grid_size - 1);
        double tau = Time / TimeSteps;

        // Create block object
        Block block(size, rank, grid_size, L, d, tau);

        if (rank == print_rank) {
                std::cout << "Process rank = " << rank << std::endl;
                block.info();
                printf("_____________________\n");
        }

        MPI_Finalize();
        return 0;
}
