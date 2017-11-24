#include <iostream>
#include <vector>


int main(int argc, char const *argv[]) {
        std::cout << "Hello world!" << '\n';
        std::vector<int> v;
        v.push_back(1);
        v.push_back(1);
        v.push_back(1);
        v.push_back(1);
        for (auto i: v) {
                std::cout << i + 1 << '\n';
        }

        return 0;
}
