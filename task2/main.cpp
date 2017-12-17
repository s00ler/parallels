#include "block.h"

int main(int argc, char** argv) {

        int size = 0, rank = 0, print_error = 0;
        // Init MPI
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        // Get args
        int grid_size = 0;
        if (argc > 1) {
                grid_size = atoi(argv[1]);
                if (argc > 2) print_error = atoi(argv[2]);
        }
        if (argc == 1 or grid_size == 0) {
                if (rank == 0) std::cout << "Invalid usage. Use: ./prog <int grid_size> [< int print_error >]\n";
                MPI_Finalize();
                return 0;
        }
        std::cout << "Processes: " << size
                  << " Grid size:" << grid_size <<'\n';

        // Create block object
        Block block(size, rank, grid_size);

        if (rank == 0) {
                std::cout << "Process rank = " << rank << std::endl;
                block.info();
                printf("_____________________\n");
        }

        // Compute
        block.compute(print_error);

        MPI_Finalize();
        return 0;
}
