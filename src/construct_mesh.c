#include <stdio.h>
#include <stdlib.h>

#include "mesh.h"
#include "check_null.h"

Mesh construct_mesh(int num_rows, int num_cols) {
    size_t num_bytes = sizeof(double) * num_rows * num_cols;

    Mesh mesh = {.num_rows = num_rows,
                 .num_cols = num_cols,
                 .x = (double *) malloc(num_bytes),
                 .y = (double *) malloc(num_bytes),
                 .z = (double *) malloc(num_bytes)};

    check_null("construct_mesh", "x", mesh.x);
    check_null("construct_mesh", "y", mesh.y);
    check_null("construct_mesh", "z", mesh.z);

    return mesh;
}