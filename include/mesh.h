#ifndef MESH_H
#define MESH_H

typedef struct Mesh {
    int num_rows;
    int num_cols;

    double *x;
    double *y;
    double *z;
} Mesh;

#endif