#include <cmath>
#include <iostream>
#include <vector>
#include <mpi.h>

#include "triplet.h"
//#include <omp.h>

// OMP flag
#define OMP_enabled False

// Hyperparameters
#define Length 1.0, 1.0, 1.0
#define Time 0.1
#define TimeSteps 200


triplet::_double L = {Length};

class Block {

triplet::_int shape {};
triplet::_int position {};

void print_block_param(triplet::_int data, const char *description = "Undescribed");

triplet::_int count_block_partition(int proc_num);

triplet::_int count_block_shape(triplet::_int block_partition, int mesh_size);

triplet::_int count_block_position(triplet::_int block_partition, int rank);

public:
Block (int proc_num, int rank, int grid_size);

void print();

};


int main(int argc, char** argv) {

        int size = 0, rank = 0, print_rank = -1;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // Get mesh size if possible
        int grid_size = 0;
        if (argc > 1) {
                grid_size = atoi(argv[1]);
                if (argc > 2) print_rank = atoi(argv[2]);
        }
        if (argc == 1 or grid_size == 0) {
                if (rank == 0) std::cout << "Invalid usage. Use: ./prog <int mesh_size> [< int rank_to_print >]\n";
                MPI_Finalize();
                return 0;
        }


        // Count grid parameters
        triplet::_double d = L / grid_size - 1;
        double dt = Time / TimeSteps;

        // Create block object
        Block block(size, rank, grid_size);

        if (rank == print_rank) {
                printf("rank=%d\n", rank);
                block.print();
                printf("_____________________\n");
        }

        MPI_Finalize();
        return 0;
}


void Block::print_block_param(triplet::_int data, const char *description) {
        std::cout << description << ": "
                  << "x=" << data.x << " | "
                  << "y=" << data.y << " | "
                  << "z=" << data.z << std::endl;
}

triplet::_int Block::count_block_partition(int proc_num) {
        triplet::_int block_partition = {1,1,1};
        int axis = 0;
        while (proc_num != 1)
        {
                switch (axis) {
                case 0: {
                        block_partition.x <<= 1;
                        break;
                }
                case 1: {
                        block_partition.y <<= 1;
                        break;
                }
                case 2: {
                        block_partition.z <<= 1;
                        break;
                }
                default: {}
                }
                proc_num >>= 1;
                axis = (axis + 1) % 3;
        }
        return block_partition;
}

triplet::_int Block::count_block_shape(triplet::_int block_partition, int mesh_size) {
        triplet::_int block_shape = {};
        block_shape.x = (mesh_size / block_partition.x);
        block_shape.y = (mesh_size / block_partition.y);
        block_shape.z = (mesh_size / block_partition.z);
        return block_shape;
}

triplet::_int Block::count_block_position(triplet::_int block_partition, int rank) {
        triplet::_int block_position = { rank % block_partition.x,
                                         (rank % (block_partition.x * block_partition.y)) / block_partition.x,
                                         rank / (block_partition.x * block_partition.y)};
        return block_position;
}

Block::Block(int proc_num, int rank, int grid_size) {
        triplet::_int partition = count_block_partition(proc_num);
        this->shape = count_block_shape(partition, grid_size);
        this->position =  count_block_position(partition, rank);
}

void Block::print() {
        print_block_param(this->shape, "Shape");
        print_block_param(this->position, "Position");
}
