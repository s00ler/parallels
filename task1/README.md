## parallels

### Task 1
Компиляция файла generator.cpp выполняется при помощи команды:
```
g++ generator.cpp -O3 -o generator -fopenmp
```
На MacBook использовалась команда:
```
g++-7 generator.cpp -O3 -o generator -fopenmp
```

В результате будет создан исполняемый файл ./generator, который может быть использован для генерации графа со следующими параметрами:

Пример 1:
```
./generator -s 20 -e 32 -directed -weighted -t 14 -file graph.bin -type RMAT
```
Данный генератор использует следующие флаги:

-s N - определяет масштаб графа, задает общее число вершин в генерируемом графе, равное 2^N. По умолчанию N установлено равным 15, в результате чего генерируется граф с 32 тыс. вершин.
-e M - определяет среднюю степень связанности каждой вершины, определяет суммарное число ребер в графе, которое будет равно M * (2^N). По умолчанию M равно 32.
-directed и -undirected - определяет, генерируется ориентированный или неориентированный граф. В случае ориентированного графа все дуги произвольны, в случае неориентированного дуги имеют вид (src_id, dst_id, weight), причем для каждой дуги src_id <= dst_id.
-weighted и -unweighted - определяет наличие весов для каждой дуги в графе. В случае взвешенного графа (weighted), каждая дуга определяется тройкой чисел (src_id, dst_id, weight), в случае невзвешенного - парой (src_id, dst_id).
-t T - число openMP потоков, используемых для генерации графа. По умолчанию данное значение устанавливается равным результату вызова функции omp_get_max_threads().
-file NAME - задает имя выходного файла, куда будет сохранен сгенерированный граф. По умолчанию выходной файл называется graph_N.bin, где N - масштаб графа.
-type TYPE - задает тип генерируемого графа. Значение TYPE может быть быть равно RMAT или SSCA2.
В результате работы исполняемый файл ./generator создаст бинарный файл следующего формата: в начале файла следует заголовок, состоящий из двух чисел int vertices_count (4 байта), long long int edges_count(8 байт) - чисел вершин и ребер графа соответственно. Затем следует граф, хранимый в формате списка ребер: edges_count троек чисел (int src_id, int dst_id, float weight), соответствующих каждому ребру. В случае невзвешенного графа вместо троек будут пары (int src_id, int dst_id).

Кроме того, на стандартный поток вывода будет выведена вспомогательная информация о числе вершин и ребер графа, а также количестве используемых потоков.

Графы должны генерироваться на целевом компьютере, где будут производиться вычисления. Это необходимо, чтобы для каждого запуска не перегонять по сети большие объёмы данных.

Пример считывания взвешенного графа из файла приведен ниже.

Пример 2. Считывание и вывод на стандартный поток взвешенного графа из сгенерированного файла
```C++
// open file
fstream file(file_name, ios::in | ios::binary);

int vertices_count = 0;
long long edges_count = 0;

// read header
file.read((char*)(&vertices_count), sizeof(int));
file.read((char*)(&edges_count), sizeof(long long));

// print graph
cout << "Graph has " << vertices_count << " vertices" << endl;
cout << "Graph has " << edges_count << " edges" << endl;

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
file.close();
```
### Using Graph500

use ```export SKIP_BFS=1``` to run only SSSP with Graph500
use ```export TMPFILE=file_name``` to set up file name
use ```export REUSEFILE=1``` to enable file reuse
