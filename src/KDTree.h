#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include <set>
#include <queue>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

// Estructura que representa una grilla
struct Grid {
    vector<pair<uint64_t, uint64_t>> ranges;  // Rangos en cada dimension [start, end]
    set<vector<uint64_t>> points;        // Conjunto de puntos en esta grilla
    uint8_t next;                       // Eje en el cual se realiza el proceso de construcción (para kdtree bi)
    uint8_t level;
    Grid(uint8_t d, uint64_t S, uint8_t level);
    Grid(); // Declaracion del constructor
};

// Estructura que representa la respuesta de la grilla
struct Response {
    bit_vector bv;                 // Bitvector con la respuesta de la grilla
    vector<Grid> subgrids;
    uint8_t level;// Subgrillas resultantes con puntos

    Response(uint8_t d, uint8_t level);               // Declaracion del constructor
};

class KDTree {
public:
    vector<vector<uint64_t>> points;  // Lista de puntos
    uint64_t S;                       // Tamaño de la grilla
    uint8_t d;                       // Dimensiones del espacio
    queue<Grid> grids;          // Arreglo de grillas
    vector<uint64_t> attr_set;
    vector<bit_vector> bitvector;


    Response (KDTree::*get_response)(Grid&);  // Mtodo de creacion de las subgrillas
    
    KDTree(const vector<vector<uint64_t>>& pts, uint64_t S, uint8_t d, const vector<uint64_t> &attr_set); // Declaración del constructor
    
    bit_vector representation; // Vector de bits que representa al arbol
    
    void build_tree(); // Constructor del vector de bits

    Response get_resp(Grid& g);

};

vector<bit_vector> join(vector<vector<bit_vector>> bitrees,vector<vector<uint64_t>> attr, uint8_t h, uint8_t d);
uint8_t intersect(vector<bit_vector>* bitvectors,vector<vector<bit_vector>> bitrees, vector<pair<uint8_t, uint8_t>> posiciones, vector<vector<uint64_t>> attr_set,uint64_t att, uint8_t p, uint8_t h);

#endif

