#include <fstream>
#include <iostream>
#include <string>

using namespace std;


int main(int argc, char const *argv[]) {
        if (argc < 3) {
                cout << "Not enough arguments." << '\n';
                cout << "Right way to use program:" << '\n';
                cout << argv[0] << " " << "source_graph_name " << "resulting_graph_name"<< '\n';
                return 0;
        }
        int vertices_count = 0;
        long long edges_count = 0;
        string new_name = argv[2];
        new_name = new_name + ".weigths";
        fstream source(argv[1], ios::in | ios::binary);
        fstream edges(argv[2], ios::out | ios::binary);
        fstream weights(new_name, ios::out | ios::binary);

        source.read((char*)(&vertices_count), sizeof(int));
        source.read((char*)(&edges_count), sizeof(long long));
        int ph = 0;
        for(long long i = 0; i < edges_count; i++)
        {
                int src_id = 0, dst_id = 0;
                float weight = 0;
                // read i-th edge data
                source.read((char*)(&src_id), sizeof(int));
                source.read((char*)(&dst_id), sizeof(int));
                source.read((char*)(&weight), sizeof(float)); // remove it for unweighed graph
                // write i-th edge data
                edges.write((char*)(&src_id), sizeof(int));
                edges.write((char*)(&dst_id), sizeof(int));
                edges.write((char*)(&ph), sizeof(int));
                weights.write((char*)(&weight), sizeof(float));  // remove it for unweighed graph
        }
        source.close();
        edges.close();
        weights.close();
        return 0;
}
