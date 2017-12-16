#ifndef TASK2_BLOCK_H
#define TASK2_BLOCK_H

#include "arrays.h"
#include <cmath>
#include <omp.h>

// OMP flag
//#define OMP_enabled true

class Block {
typedef struct {Array3D *prev; Array3D *curr; Array3D *next;} Values;

typedef struct {Array3D *x_prev; Array3D *y_prev; Array3D *z_prev;
                Array3D *x_next; Array3D *y_next; Array3D *z_next;} Bounds;

typedef struct {int x_prev; int y_prev; int z_prev;
                int x_next; int y_next; int z_next;} Neighbours;

triplet<int> shape {};
triplet<int> position {};
triplet<int> partition {};
triplet<double> L {};
triplet<double> d {};

double tau = 0.0;

Values values {};

Bounds bound {};
Bounds bound_buf {};

Neighbours neighbour {};

void print_block_param(triplet<int> data, const char *description = "Undescribed");

triplet<int> count_block_partition(int proc_num);

triplet<int> count_block_shape(triplet<int> block_partition, int grid_size);

triplet<int> count_block_position(triplet<int> block_partition, int rank);

void init_bounds();

void init_bounds_buf();

int get_process_id(triplet<int> position);

int get_process_id(int x, int y, int z);

void set_neighbours_process_id();

void init_values();

double phi (double x, double y, double z);

void first_step();

triplet<double> count_laplassian(int i, int j, int k);

void second_step();

public:

Block (int proc_num, int rank, int grid_size, triplet<double> L, triplet<double> d, double tau);

void info();

};





//**************************************************************************
// Implementations
//**************************************************************************

Block::Block(int proc_num, int rank, int grid_size, triplet<double> L, triplet<double> d, double tau) {
        partition = count_block_partition(proc_num);
        shape = count_block_shape(partition, grid_size);
        position =  count_block_position(partition, rank);
        this->L = L;
        this->d = d;
        this->tau = tau;

        init_values();
        init_bounds();
        init_bounds_buf();
        set_neighbours_process_id();

        first_step();
        values.prev = values.curr;
        values.curr = values.next;
        second_step();
        values.curr->print();
}

void Block::init_bounds() {
        bound.x_prev = new Array3D(shape.y, shape.z);
        bound.x_next = new Array3D(shape.y, shape.z);
        bound.y_prev = new Array3D(shape.x, shape.z);
        bound.y_next = new Array3D(shape.x, shape.z);
        bound.z_prev = new Array3D(shape.x, shape.y);
        bound.z_next = new Array3D(shape.x, shape.y);
}

void Block::init_bounds_buf() {

        bound_buf.x_prev = new Array3D(shape.y, shape.z);
        bound_buf.x_next = new Array3D(shape.y, shape.z);
        bound_buf.y_prev = new Array3D(shape.x, shape.z);
        bound_buf.y_next = new Array3D(shape.x, shape.z);
        bound_buf.z_prev = new Array3D(shape.x, shape.y);
        bound_buf.z_next = new Array3D(shape.x, shape.y);
}

triplet<int> Block::count_block_partition(int proc_num) {
        triplet<int> block_partition {1,1,1};
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


triplet<int> Block::count_block_shape(triplet<int> block_partition, int grid_size) {
        triplet<int> block_shape {};
        block_shape.x = (grid_size / block_partition.x);
        block_shape.y = (grid_size / block_partition.y);
        block_shape.z = (grid_size / block_partition.z);
        return block_shape;
}

triplet<int> Block::count_block_position(triplet<int> block_partition, int rank) {
        triplet<int> block_position = {rank % block_partition.x,
                                       (rank % (block_partition.x * block_partition.y)) / block_partition.x,
                                       rank / (block_partition.x * block_partition.y)};
        return block_position;
}

void Block::print_block_param(triplet<int> data, const char *description) {
        std::cout << description << ": " << data << std::endl;
}

void Block::info() {
        print_block_param(partition, "Partition");
        print_block_param(shape, "Shape");
        print_block_param(position, "Position");
        std::cout << "Proc_id: " << get_process_id(position) << std::endl;
        std::cout << "Neighbour processes: \n -------------------\n"
                  << "x prev = " << neighbour.x_prev
                  << ", " << " x next = " << neighbour.x_next << std::endl
                  << "y prev = " << neighbour.y_prev
                  << ", " << " y next = " << neighbour.y_next << std::endl
                  << "z prev = " << neighbour.z_prev
                  << ", " << " z next = " << neighbour.z_next << std::endl;

}

int Block::get_process_id(triplet<int> position) {
        return get_process_id(position.x, position.y, position.z);
}

int Block::get_process_id(int x, int y, int z) {
        return x + y * partition.x + z * partition.x * partition.y;
}

void Block::init_values() {
        values.prev = new Array3D(shape+2);
        values.curr = new Array3D(shape+2);
        values.next = new Array3D(shape+2);
}

void Block::set_neighbours_process_id() {
        if (position.x == 0) {
                neighbour.x_prev = get_process_id(partition.x - 1, position.y, position.z);
                neighbour.x_next = get_process_id(position.x + 1, position.y, position.z);
        }
        neighbour.x_prev = get_process_id(position.x - 1, position.y, position.z);
        neighbour.x_next = get_process_id(position.x == (partition.x - 1) ? 0 : position.x + 1, position.y, position.z);

        if (position.y == 0) {
                neighbour.y_prev = get_process_id(position.x, partition.y - 1, position.z);
                neighbour.y_next = get_process_id(position.x, position.y + 1, position.z);
        }
        neighbour.y_prev  = get_process_id(position.x, position.y - 1, position.z);
        neighbour.y_next = get_process_id(position.x, position.y == (partition.y - 1) ? 0 : position.y + 1, position.z);

        if (position.z == 0) {
                neighbour.z_prev = get_process_id(position.x, position.y, partition.z - 1);
                neighbour.z_next = get_process_id(position.x, position.y, position.z + 1);
        }
        neighbour.z_prev  = get_process_id(position.x, position.y, position.z - 1);
        neighbour.z_next = get_process_id(position.x, position.y, position.z == (partition.z - 1) ? 0 : position.z + 1);
}

double Block::phi(double x, double y, double z) {
        return sin(2* M_PI * x / L.x)
               * sin(2 * M_PI * y / L.y)
               * sin(2 * M_PI * z / L.z);
}

void Block::first_step() {
//#pragma omp parallel for if (OMP_enabled)
        for (int k = position.z * shape.z + 1; k < position.z * shape.z + shape.z + 1; ++k)
                for (int j = position.y * shape.y + 1; j < position.y * shape.y + shape.y + 1; ++j)
                        for (int i = position.x * shape.x + 1; i < position.x * shape.x + shape.x + 1; ++i) {
                                double t = phi((i - 1) * d.x, (j - 1) * d.y, (k - 1) * d.z);
                                *values.next->iloc(i, j, k) = phi((i - 1) * d.x, (j - 1) * d.y, (k - 1) * d.z);
                        }

}

triplet<double> Block::count_laplassian(int i, int j, int k) {
        triplet<double> laplassian {};

        laplassian.x = (*values.curr->iloc(i - 1, j,k)
                        - 2 * *values.curr->iloc(i, j,k)
                        + *values.curr->iloc(i + 1, j,k))
                       / (d.x * d.x);

        laplassian.y = (*values.curr->iloc(i, j-1,k)
                        - 2 * *values.curr->iloc(i, j,k)
                        + *values.curr->iloc(i, j + 1,k))
                       / (d.y * d.y);

        laplassian.z = (*values.curr->iloc(i, j,k - 1)
                        - 2 * *values.curr->iloc(i, j,k)
                        + *values.curr->iloc(i, j,k + 1))
                       / (d.z * d.z);

        return laplassian;
}

void Block::second_step() {
//#pragma omp parallel for if (OMP_enabled)
        for (int k = position.z * shape.z + 1; k < position.z * shape.z + shape.z + 1; ++k)
                for (int j = position.y * shape.y + 1; j < position.y * shape.y + shape.y + 1; ++j)
                        for (int i = position.x * shape.x + 1; i < position.x * shape.x + shape.x + 1; ++i) {
                                triplet<double> laplassian = count_laplassian(i, j, k);
                                *values.next->iloc(i, j, k) = *values.curr->iloc(i, j, k)
                                                              + (laplassian.x + laplassian.y + laplassian.z) * tau * tau /2;
                        }

}

#endif //TASK2_BLOCK_H
