#ifndef TASK2_ARRAY3D_H
#define TASK2_ARRAY3D_H

#include "triplet.h"

#define OMP_enabled true

struct Array3D {

        double *data {};

        triplet<int> shape {};

        Array3D(triplet<int> shape);

        Array3D(int x, int y, int z);

        double operator[] (triplet<int> idx);

        double *iloc(int x, int y, int z);

        void print();

        Array3D& operator= (Array3D other);
};

//**************************************************************************
// Implementations
//**************************************************************************

Array3D::Array3D(triplet<int> shape) {
        data = new double [shape.x * shape.y * shape.z] {};
        this->shape = shape;
}

Array3D::Array3D(int x = 1, int y = 1, int z = 1) {
        data = new double [x * y * z] {};
        shape = triplet<int> {x, y, z};
}

double Array3D::operator[](triplet<int> idx) {
        return data[idx.x + (idx.y * shape.x) + (idx.z * shape.x * shape.y)];
}

double *Array3D::iloc(int x = 0, int y = 0, int z = 0) {
        return &data[x + (y * shape.x) + (z * shape.x * shape.y)];
}

Array3D& Array3D::operator= (Array3D other) {
#pragma omp parallel for if (OMP_enabled)
        for (size_t i = 0; i < other.shape.x * other.shape.y * other.shape.z; ++i)
                this->data[i] = other.data[i];
        return *this;
}

void Array3D::print() {
    std::cout << "Array shape: " << shape << std::endl;
        for (int i = 0; i < shape.x; i++ ) {
                for (int j = 0; j < shape.y; j++ ) {
                        for (int k = 0; k < shape.z; k++ ) {
                                std::cout << *iloc(i, j, k) << " ";
                        }
                        std::cout << std::endl;
                }
                std::cout << "_____________________" << std::endl;
        }
}

#endif //TASK2_ARRAY3D_H
