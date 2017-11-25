#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <ctime>

int main( int argc, char* argv[] )
{
        omp_set_num_threads( 8 );

        double pi = acos( -1.0 );

        std::cout << "Allocating memory ..." << std::endl;
        std::vector<double> my_vector(128000000, 0.0);
        std::cout << "Done!" << std::endl << std::endl;

        std::cout << "Entering main loop ... " << std::endl;
        std::time_t start = std::time(nullptr);
     #pragma omp parallel for
        for (size_t i=0; i < my_vector.size(); ++i )
        {
                my_vector[i] = exp( -sin( i*i + pi*log(i+1) ) );
        }
        std::cout << "Done!" << std::endl;
        std::time_t end = std::time(nullptr);
        std::cout << "Time: " << end - start << '\n';

        std::cout << "Entering main loop ... " << std::endl;
        start = std::time(nullptr);
        for (size_t i=0; i < my_vector.size(); ++i )
        {
                my_vector[i] = exp( -sin( i*i + pi*log(i+1) ) );
        }
        std::cout << "Done!" << std::endl;
        end = std::time(nullptr);
        std::cout << "Time: " << end - start << '\n';
        return 0;
}
