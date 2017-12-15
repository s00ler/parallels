#ifndef TASK2_BLOCK_H
#define TASK2_BLOCK_H

#include "arrays.h"

class Block {

    int_triplet shape {};
    int_triplet position {};
    int_triplet partition {};

    double *prev_values{};
    double *curr_values{};
    double *next_values{};

    void print_block_param(int_triplet data, const char *description = "Undescribed");

    int_triplet count_block_partition(int proc_num);

    int_triplet count_block_shape(int_triplet block_partition, int grid_size);

    int_triplet count_block_position(int_triplet block_partition, int rank);

public:
    Block (int proc_num, int rank, int grid_size);

    void print();

};

//**************************************************************************
// Implementations
//**************************************************************************


void Block::print_block_param(int_triplet data, const char *description) {
    std::cout << description << ": " << data << std::endl;
}

int_triplet Block::count_block_partition(int proc_num) {
    int_triplet block_partition = {1,1,1};
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

int_triplet Block::count_block_shape(int_triplet block_partition, int grid_size) {
    int_triplet block_shape = {};
    block_shape.x = (grid_size / block_partition.x);
    block_shape.y = (grid_size / block_partition.y);
    block_shape.z = (grid_size / block_partition.z);
    return block_shape;
}

int_triplet Block::count_block_position(int_triplet block_partition, int rank) {
    int_triplet block_position = {rank % block_partition.x,
                                 (rank % (block_partition.x * block_partition.y)) / block_partition.x,
                                  rank / (block_partition.x * block_partition.y)};
    return block_position;
}

Block::Block(int proc_num, int rank, int grid_size) {
    this->partition = count_block_partition(proc_num);
    this->shape = count_block_shape(this->partition, grid_size);
    this->position =  count_block_position(this->partition, rank);
}

void Block::print() {
    print_block_param(this->partition, "Partition");
    print_block_param(this->shape, "Shape");
    print_block_param(this->position, "Position");
}


#endif //TASK2_BLOCK_H
