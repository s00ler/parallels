#ifndef TASK2_TRIPLET_H
#define TASK2_TRIPLET_H

#include <iostream>

template <typename T>
struct triplet {
    T x {}, y {}, z {};

    triplet(T x = 0, T y = 0, T z = 0) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    triplet operator* (T value) {
        return {x * value,
                y * value,
                z * value};
    }

    triplet operator/ (T value) {
        return {x / value,
                y / value,
                z / value};
    }

    triplet operator+ (T value) {
        return {x + value,
                y + value,
                z + value};
    }

    triplet operator- (T value) {
        return {x - value,
                y - value,
                z - value};
    }

    void print() {
        std::cout << *this << std::endl;
    }

    template <typename CT>
    friend std::ostream& operator<< (std::ostream& os, const triplet<CT>& obj);
};

template <typename CT>
std::ostream& operator<< (std::ostream& os, const triplet<CT>& obj) {
    os << "x = " << obj.x << ", "
       << "y = " << obj.y << ", "
       << "z = " << obj.z;
    return os;
}

#endif //TASK2_TRIPLET_H
