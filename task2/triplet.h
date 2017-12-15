#ifndef TASK2_TRIPLET_H
#define TASK2_TRIPLET_H

#include <iostream>

struct int_triplet{
    int x, y, z;

    int_triplet(int x = 0, int y = 0, int z = 0) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    int_triplet operator* (int value) {
        return {this->x * value,
                this->y * value,
                this->z * value};
    }

    int_triplet operator/ (int value) {
        return {this->x / value,
                this->y / value,
                this->z / value};
    }

    int_triplet operator+ (int value) {
        return {this->x + value,
                this->y + value,
                this->z + value};
    }

    int_triplet operator- (int value) {
        return {this->x - value,
                this->y - value,
                this->z - value};
    }

    int_triplet *incx() {
        this->x++;
        return this;
    }

    int_triplet *incy() {
        this->y++;
        return this;
    }

    int_triplet *incz() {
        this->z++;
        return this;
    }

    friend std::ostream& operator<< (std::ostream& os, const int_triplet& obj);
};


struct dbl_triplet {
    double x, y, z;

    dbl_triplet(double x = 0.0, double y = 0.0, double z = 0.0) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    dbl_triplet operator* (double value) {
        return {this->x * value,
                this->y * value,
                this->z * value};
    }

    dbl_triplet operator/ (double value) {
        return {this->x / value,
                this->y / value,
                this->z / value};
    }

    dbl_triplet operator+ (double value) {
        return {this->x + value,
                this->y + value,
                this->z + value};
    }

    dbl_triplet operator- (double value) {
        return {this->x - value,
                this->y - value,
                this->z - value};
    }

    friend std::ostream& operator<< (std::ostream& os, const dbl_triplet& obj);
};

std::ostream& operator<<(std::ostream& os, const int_triplet& obj) {
    os << "x = " << obj.x << ", "
       << "y = " << obj.y << ", "
       << "z = " << obj.z;
    return os;
}

std::ostream& operator<<(std::ostream& os, const dbl_triplet& obj) {
    os << "x = " << obj.x << ", "
       << "y = " << obj.y << ", "
       << "z = " << obj.z;
    return os;
}
#endif //TASK2_TRIPLET_H
