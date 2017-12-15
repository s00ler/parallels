#ifndef TASK2_ARRAY3D_H
#define TASK2_ARRAY3D_H

#include "triplet.h"

class Array3D {

    double *data {};
    int_triplet shape {};

public:
    Array3D(int_triplet shape);

    double operator[] (int_triplet idx);

    void print();
};

//**************************************************************************
// Implementations
//**************************************************************************

Array3D::Array3D(int_triplet shape) {
    this->data = new double [shape.x * shape.y * shape.z] {};
    this->shape = shape;
}

double Array3D::operator[](int_triplet idx) {
    return this->data[idx.x + (idx.y * this->shape.x) + (idx.z * this->shape.x * this->shape.y)];
}

void Array3D::print() {
    for (int i = 0; i < this->shape.x; i++ ) {
        for (int j = 0; j < this->shape.y; j++ ) {
            for (int k = 0; k < this->shape.z; k++ ) {
                std::cout << this->operator[](int_triplet{i, j ,k}) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "_____________________" << std::endl;
    }
}

#endif //TASK2_ARRAY3D_H
