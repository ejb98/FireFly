#include <stdio.h>
#include <stdlib.h>

#include "mesh.h"

Mesh construct_mesh(int num_rows, int num_cols) {
    size_t num_bytes = sizeof(double) * num_rows * num_cols;

    Mesh mesh = {.num_rows = num_rows,
                 .num_cols = num_cols,
                 .x = (double *) malloc(num_bytes),
                 .y = (double *) malloc(num_bytes),
                 .z = (double *) malloc(num_bytes)};

    return mesh;
}