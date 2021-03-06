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
                if (rank == 0) std::cout << "Invalid usage. Use: ./prog <int grid_size> [< int print_rank >]\n";
                MPI_Finalize();
                return 0;
        }
        if (rank == 0) {
            std::cout << "Processes: " << size << "  "
                      << "Grid size: " << grid_size << "\n";
        }
        // Create block object
        Block block(size, rank, grid_size);

        for (int proc_rank = 0; proc_rank < size; proc_rank++)
                if (rank == proc_rank and print_error == 1) {
                        block.info();
                        printf("_____________________\n");
                }

        // Compute
        block.compute(print_error);

        MPI_Finalize();
//        Array3D tmpx = {1,3,4};
//        Array3D tmpy = {3,1,4};
//        Array3D tmpz = {3,4,1};
//
//        tmpx.print();
//        tmpy.print();
//        tmpz.print();
        return 0;
}
