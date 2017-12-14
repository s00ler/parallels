#ifndef TASK2_TRIPLET_H
#define TASK2_TRIPLET_H
namespace triplet {
struct _int {
        int x, y, z;

        _int(int x = 0, int y = 0, int z = 0) {
                this->x = x;
                this->y = y;
                this->z = z;
        }

        _int operator*(int value) {
                return {this->x * value,
                        this->y * value,
                        this->z * value};
        }

        _int operator/(int value) {
                return {this->x / value,
                        this->y / value,
                        this->z / value};
        }

        _int operator+(int value) {
                return {this->x + value,
                        this->y + value,
                        this->z + value};
        }

        _int operator-(int value) {
                return {this->x - value,
                        this->y - value,
                        this->z - value};
        }
};

struct _double {
        double x, y, z;

        _double(double x = 0.0, double y = 0.0, double z = 0.0) {
                this->x = x;
                this->y = y;
                this->z = z;
        }

        _double operator*(double value) {
                return {this->x * value,
                        this->y * value,
                        this->z * value};
        }

        _double operator/(double value) {
                return {this->x / value,
                        this->y / value,
                        this->z / value};
        }

        _double operator+(double value) {
                return {this->x + value,
                        this->y + value,
                        this->z + value};
        }

        _double operator-(double value) {
                return {this->x - value,
                        this->y - value,
                        this->z - value};
        }
};
}
#endif //TASK2_TRIPLET_H
