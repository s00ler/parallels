#ifndef TASK2_BLOCK_H
#define TASK2_BLOCK_H

#include "arrays.h"
#include <cmath>
#include <omp.h>
#include <mpi.h>


// OMP flag
#define OMP_enabled true

// Hyperparameters
#define Length 2.0, 1.0, 2.0
#define Time 0.1
#define TimeSteps 200

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

int rank {};

double tau {};

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

void next_step();

void set_periodic_boundaries();

void transfer_data_between_processes();

void apply_periodic_boundaries();

public:

Block (int proc_num, int rank, int grid_size);

void info();

void compute(int print_proc);

};





//**************************************************************************
// Implementations
//**************************************************************************

Block::Block(int proc_num, int rank, int grid_size) {
        partition = count_block_partition(proc_num);
        shape = count_block_shape(partition, grid_size);
        position =  count_block_position(partition, rank);
        this->rank = rank;
        this->L = triplet<double> {Length};
        this->d = L / (grid_size - 1);;
        tau = Time / TimeSteps;

        init_values();
        init_bounds();
        init_bounds_buf();
        set_neighbours_process_id();
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

        neighbour.x_prev = get_process_id(position.x == 0 ? partition.x - 1 : position.x - 1, position.y, position.z);
        neighbour.x_next = get_process_id(position.x == (partition.x - 1) ? 0 : position.x + 1, position.y, position.z);

        neighbour.y_prev = get_process_id(position.x, position.y == 0 ? partition.y - 1 : position.y - 1, position.z);
        neighbour.y_next = get_process_id(position.x, position.y == (partition.y - 1) ? 0 : position.y + 1, position.z);

        neighbour.z_prev = get_process_id(position.x, position.y, position.z == 0 ? partition.z - 1 : position.z - 1);
        neighbour.z_next = get_process_id(position.x, position.y, position.z == (partition.z - 1) ? 0 : position.z + 1);

}

double Block::phi(double x, double y, double z) {
        return sin(2* M_PI * x / L.x)
               * sin(2 * M_PI * y / L.y)
               * sin(2 * M_PI * z / L.z);
}

void Block::first_step() {
#pragma omp parallel for if (OMP_enabled)
        for (int k = 1; k < shape.z + 1; ++k)
                for (int j = 1; j < shape.y + 1; ++j)
                        for (int i = 1; i < shape.x + 1; ++i)
                                *values.next->iloc(i, j, k) = phi((position.x * shape.x + i - 1) * d.x,
                                                                  (position.y * shape.y + j - 1) * d.y,
                                                                  (position.z * shape.z + k - 1) * d.z);

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
#pragma omp parallel for if (OMP_enabled)
        for (int k = 1; k < shape.z + 1; ++k)
                for (int j = 1; j < shape.y + 1; ++j)
                        for (int i = 1; i < shape.x + 1; ++i) {
                                triplet<double> laplassian = count_laplassian(i, j, k);
                                *values.next->iloc(i, j, k) = *values.curr->iloc(i, j, k)
                                                              + (laplassian.x + laplassian.y + laplassian.z)
                                                              * tau * tau / 2;
                        }

}

void Block::compute(int print_proc = -1) {
        double start_time = MPI_Wtime();
        double iteration_finish_time {}, iteration_start_time {};
        for (int ts = 0; ts <= TimeSteps; ++ts) {
                iteration_start_time = MPI_Wtime();

                if (ts == 0) first_step();
                if (ts == 1) second_step();
                if (ts >= 2) next_step();
                set_periodic_boundaries();
                transfer_data_between_processes();
                apply_periodic_boundaries();

                values.prev = values.curr;
                values.curr = values.next;

                iteration_finish_time = MPI_Wtime();
                if (rank == -print_proc) {
                        std::cout << "Time step " << ts << "/" << TimeSteps << " done.\n"
                                  << "Execution time: " <<  iteration_finish_time - iteration_start_time
                                  << "\n-----------------------\n";
                }
        }
        if (rank == 0) {
                std::cout << "Total execution time: " << MPI_Wtime() - start_time << std::endl;
        }

}

void Block::next_step() {
#pragma omp parallel for if (OMP_enabled)
        for (int k = 1; k < shape.z + 1; ++k)
                for (int j = 1; j < shape.y + 1; ++j)
                        for (int i = 1; i < shape.x + 1; ++i) {
                                triplet<double> laplassian = count_laplassian(i, j, k);
                                *values.next->iloc(i, j, k) = *values.curr->iloc(i, j, k) * 2
                                                              - *values.prev->iloc(i, j, k)
                                                              + (laplassian.x + laplassian.y + laplassian.z)
                                                              * tau * tau / 2;
                        }
}

void Block::set_periodic_boundaries() {
#pragma omp parallel for if (OMP_enabled)
        for (int k = 1; k < shape.z; ++k)
                for (int j = 1; j < shape.y; ++j) {
                        *bound.x_prev->iloc(j - 1, k - 1) = *values.next->iloc(position.x == 0 ? 1 : 2, j, k);
                        *bound.x_next->iloc(j - 1, k - 1) = *values.next->iloc(position.x == partition.x - 1
                                                                               ? shape.x : shape.x -1, j, k);
                }
#pragma omp parallel for if (OMP_enabled)
        for (int k = 1; k < shape.z; ++k)
                for (int i = 1; i < shape.x; ++i) {
                        *bound.y_prev->iloc(i - 1, k - 1) = *values.next->iloc(i, position.y == 0 ? 1 : 2, k);
                        *bound.y_next->iloc(i - 1, k - 1) = *values.next->iloc(1, position.y == partition.y - 1
                                                                               ? shape.y : shape.y -1, k);
                }
#pragma omp parallel for if (OMP_enabled)
        for (int j = 1; j < shape.y; ++j)
                for (int i = 1; i < shape.x; ++i) {
                        *bound.z_prev->iloc(i - 1, j - 1) = *values.next->iloc(i, j, position.z == 0 ? 1 : 2);
                        *bound.z_next->iloc(i - 1, j - 1) = *values.next->iloc(i, j, position.z == partition.z - 1
                                                                               ? shape.z : shape.z -1);
                }
}

void Block::transfer_data_between_processes() {
        MPI_Status status {};
        MPI_Request requests[12];

        MPI_Irecv(bound_buf.x_prev->data, shape.y * shape.z, MPI_DOUBLE,
                  neighbour.x_prev, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Isend(bound.x_prev->data, shape.y * shape.z, MPI_DOUBLE,
                  neighbour.x_prev, 1, MPI_COMM_WORLD, &requests[1]);

        MPI_Irecv(bound_buf.x_next->data, shape.y * shape.z, MPI_DOUBLE,
                  neighbour.x_next, 1, MPI_COMM_WORLD, &requests[2]);
        MPI_Isend(bound.x_next->data, shape.y * shape.z, MPI_DOUBLE,
                  neighbour.x_next, 0, MPI_COMM_WORLD, &requests[3]);

        MPI_Irecv(bound_buf.y_prev->data, shape.x * shape.z, MPI_DOUBLE,
                  neighbour.y_prev, 2, MPI_COMM_WORLD, &requests[4]);
        MPI_Isend(bound.y_prev->data, shape.x * shape.z, MPI_DOUBLE,
                  neighbour.y_prev, 3, MPI_COMM_WORLD, &requests[5]);


        MPI_Irecv(bound_buf.y_next->data, shape.x * shape.z, MPI_DOUBLE,
                  neighbour.y_next, 3, MPI_COMM_WORLD, &requests[6]);
        MPI_Isend(bound.y_next->data, shape.x * shape.z, MPI_DOUBLE,
                  neighbour.y_next, 2, MPI_COMM_WORLD, &requests[7]);


        MPI_Irecv(bound_buf.z_prev->data, shape.x * shape.y, MPI_DOUBLE,
                  neighbour.z_prev, 4, MPI_COMM_WORLD, &requests[8]);
        MPI_Isend(bound.z_prev->data, shape.x * shape.y, MPI_DOUBLE,
                  neighbour.z_prev, 5, MPI_COMM_WORLD, &requests[9]);


        MPI_Irecv(bound_buf.z_next->data, shape.x * shape.y, MPI_DOUBLE,
                  neighbour.z_next, 5, MPI_COMM_WORLD, &requests[10]);
        MPI_Isend(bound.z_next->data, shape.x * shape.y, MPI_DOUBLE,
                  neighbour.z_next, 4, MPI_COMM_WORLD, &requests[11]);

        for (auto &request : requests)
                MPI_Wait(&request, &status);
}

void Block::apply_periodic_boundaries() {
#pragma omp parallel for if (OMP_enabled)
        for (int k = 1; k < shape.z + 1; ++k)
                for (int j = 1; j < shape.y + 1; ++j) {
                        *values.next->iloc(0, j, k)           = *bound_buf.x_prev->iloc(j - 1, k - 1);
                        *values.next->iloc(shape.x + 1, j, k) = *bound_buf.x_next->iloc(j - 1, k - 1);
                }
#pragma omp parallel for if (OMP_enabled)
        for (int k = 1; k < shape.z + 1; ++k)
                for (int i = 1; i < shape.x + 1; ++i) {
                        *values.next->iloc(i, 0, k)           = *bound_buf.y_prev->iloc(i - 1, k - 1);
                        *values.next->iloc(i, shape.y + 1, k) = *bound_buf.y_next->iloc(i - 1, k - 1);
                }
#pragma omp parallel for if (OMP_enabled)
        for (int j = 1; j < shape.y + 1; ++j)
                for (int i = 1; i < shape.x + 1; ++i) {
                        *values.next->iloc(i, j, 0)           = *bound_buf.z_prev->iloc(i - 1, j - 1);
                        *values.next->iloc(i, j, shape.z + 1) = *bound_buf.z_next->iloc(i - 1, j - 1);
                }
}

#endif //TASK2_BLOCK_H
