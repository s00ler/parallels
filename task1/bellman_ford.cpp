#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <ctime>
using namespace boost;

struct EdgeProperties {
        float weight;
};

int main(int argc, char const *argv[]) {
        std::cout << "Start!" << std::endl<< std::flush;
        typedef std::pair < int, int > E;
        int N = 262144;
        const long long n_edges = 8388608;
        long long phdr = 0;
        static E edge_array[n_edges];
        static float weight_array[n_edges];
        std::fstream source(argv[1], std::ios::in | std::ios::binary);
        std::cout << "File init" << std::endl<< std::flush;

        source.read((char*) &N, sizeof(int));
        source.read((char*) &phdr, sizeof(long long));
        std::cout << "Read Begin!" << std::endl;
        for(long long i = 0; i < n_edges; ++i) {
                int src_id = 0, dst_id = 0;
                float weight = 0;
                // read i-th edge data
                source.read((char *) (&src_id), sizeof(int));
                source.read((char *) (&dst_id), sizeof(int));
                source.read((char *) (&weight), sizeof(float));
                edge_array[i] = E(src_id, dst_id);
                weight_array[i] = weight;
        }
        std::cout << "Read done!" << std::endl<< std::flush;

        typedef adjacency_list < vecS, vecS, directedS,
                                 no_property, EdgeProperties> Graph;

        static Graph g(edge_array, edge_array + n_edges, N);

        graph_traits < Graph >::edge_iterator ei, ei_end;

        property_map < Graph, float EdgeProperties::* >::type weight_pmap = get(&EdgeProperties::weight, g);

        int i = 0;
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei, ++i)
                weight_pmap[*ei] = weight_array[i];

        std::vector<int> distance(N, (std::numeric_limits < short >::max)());
        std::vector<std::size_t> parent(N);
        for (i = 0; i < N; ++i)
                parent[i] = i;
        distance[262143] = 0;
        std::vector<float> Teps;
        std::cout << "Bellman-Ford start!" << std::endl<< std::flush;
        for (int i = 0; i < 10; i++) {
                clock_t start = clock();
                bool r = bellman_ford_shortest_paths
                                 (g, N, weight_map(weight_pmap).distance_map(&distance[0]).
                                 predecessor_map(&parent[0]));
                clock_t end = clock() - start;
                Teps.push_back( (float)n_edges / ((float)end / CLOCKS_PER_SEC));
        }
        float sum = 0.0;
        for (auto a: Teps) {
                sum += a;
        }
        std::cout << "avg Teps: " << sum / Teps.size() << "\n";
        return 0;
}
