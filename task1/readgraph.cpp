#include <fstream>
#include <iostream>
#include <string>

using namespace std;


int main(int argc, char const *argv[]) {
        int vertices_count = 0;
        long long edges_count = 0;
        int kind = atoi(argv[2]); // 0 if generator graph used, 1 if Graph 500 graph
        if (kind == 0) {
                // read header
                fstream file(argv[1], ios::in | ios::binary);
                file.read((char*)(&vertices_count), sizeof(int));
                file.read((char*)(&edges_count), sizeof(long long));
                // print graph
                // get & print graph data for WEIGHTED graph
                for(long long i = 0; i < edges_count; i++)
                {
                        int src_id = 0, dst_id = 0;
                        float weight = 0;
                        // read i-th edge data
                        file.read((char*)(&src_id), sizeof(int));
                        file.read((char*)(&dst_id), sizeof(int));
                        file.read((char*)(&weight), sizeof(float)); // remove it for unweighed graph
                        //print edge data
                        cout << src_id << " " << dst_id << " | " << weight << endl;
                }

                cout << "Graph has " << vertices_count << " vertices" << endl;
                cout << "Graph has " << edges_count << " edges" << endl;
                file.close();
        }
        if (kind == 1) {
                long long counter = 0;
                string new_name = argv[1];
                new_name = new_name + ".weigths";
                int tmp = 0;
                fstream file(argv[1], ios::in | ios::binary);
                fstream fileweights(new_name, ios::in | ios::binary);
                do {
                        counter++;
                        int src_id = 0, dst_id = 0;
                        float weight = 0;

                        // // read i-th edge data
                        file.read((char*)(&src_id), sizeof(int32_t));
                        file.read((char*)(&dst_id), sizeof(int32_t));
                        file.read((char*)(&tmp), sizeof(int32_t));
                        fileweights.read((char*)(&weight), sizeof(float)); // remove it for unweighed graph
                        cout << counter << ": "<< src_id << "  " << dst_id << " | " << weight << endl;

                } while(!file.eof());
                file.close();
                fileweights.close();
                cout << "Graph has " << counter - 1 << " edges" << endl;
        }
        return 0;
}
