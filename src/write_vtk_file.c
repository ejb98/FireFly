#include <stdio.h>
#include <stdlib.h>

#include "mesh.h"
#include "isempty.h"
#include "sub2ind.h"

void write_vtk_file(const Mesh *mesh, const char *file_name) {
    FILE *fp = fopen(file_name, "w");

    if (fp == NULL) {
        fprintf(stderr, "failed to open file");

        return;
    }

    if (isempty(mesh)) {
        fprintf(stderr, "write_to_vtk_file: unable to write empty mesh to file");
        fclose(fp);

        return;
    }

    size_t num_rows = (size_t) mesh->num_rows;
    if (num_rows < 2) {
        fprintf(stderr, "write_to_vtk_file: mesh must have at least 2 rows");
        fclose(fp);

        return;
    }

    size_t num_cols = (size_t) mesh->num_cols;
    if (num_cols < 2) {
        fprintf(stderr, "write_to_vtk_file: mesh must have at least 2 columns");
        fclose(fp);

        return;
    }

    size_t i0, i1, i2, i3;
    size_t num_quads = (num_rows - 1) * (num_cols - 1);
    size_t num_points = num_rows * num_cols;

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Mesh Surface\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %zu float\n", num_points);

    for (size_t i = 0; i < num_points; i++) {
        fprintf(fp, "%f %f %f\n", mesh->x[i], mesh->y[i], mesh->z[i]);
    }

    fprintf(fp, "CELLS %zu %zu\n", num_quads, 5 * num_quads);
    for (size_t i = 0; i < num_rows - 1; i++) {
        for (size_t j = 0; j < num_cols - 1; j++) {
            i0 = i * num_cols + j;
            i1 = i * num_cols + j + 1;
            i2 = (i + 1) * num_cols + j + 1;
            i3 = (i + 1) * num_cols + j;

            fprintf(fp, "4 %zu %zu %zu %zu\n", i0, i1, i2, i3);
        }
    }

    fprintf(fp, "CELL_TYPES %zu\n", num_quads);
    for (size_t i = 0; i < num_quads; i++) {
        fprintf(fp, "9\n");
    }

    fclose(fp);
}