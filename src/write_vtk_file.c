#include <stdio.h>
#include <stdlib.h>

#include "mesh.h"
#include "sub2ind.h"
#include "get_size.h"

void write_vtk_file(const Mesh *mesh, const char *file_name) {
    FILE *fp = fopen(file_name, "w");

    if (fp == NULL) {
        fprintf(stderr, "write_to_vtk_file: failed to open %s", file_name);

        return;
    }

    if (mesh->x == NULL || mesh->y == NULL || mesh->z == NULL) {
        fprintf(stderr, "write_to_vtk_file: unable to write empty mesh to %s", file_name);
        fclose(fp);

        return;
    }

    if (mesh->num_rows < 2) {
        fprintf(stderr, "write_to_vtk_file: mesh must have at least 2 rows");
        fclose(fp);

        return;
    }

    if (mesh->num_cols < 2) {
        fprintf(stderr, "write_to_vtk_file: mesh must have at least 2 columns");
        fclose(fp);

        return;
    }

    size_t i0, i1, i2, i3;
    size_t num_quads = ((size_t) mesh->num_rows - 1) * (mesh->num_cols - 1);
    size_t num_points = get_size(mesh);

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Mesh Surface\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp, "POINTS %zu float\n", num_points);
    for (size_t i = 0; i < num_points; i++) {
        fprintf(fp, "%f %f %f\n", mesh->x[i], mesh->y[i], mesh->z[i]);
    }

    fprintf(fp, "CELLS %zu %zu\n", num_quads, 5 * num_quads);
    for (int i = 0; i < mesh->num_rows - 1; i++) {
        for (int j = 0; j < mesh->num_cols - 1; j++) {
            i0 = sub2ind(i, j, mesh->num_cols);
            i1 = sub2ind(i, j + 1, mesh->num_cols);
            i2 = sub2ind(i + 1, j + 1, mesh->num_cols);
            i3 = sub2ind(i + 1, j, mesh->num_cols);

            fprintf(fp, "4 %zu %zu %zu %zu\n", i0, i1, i2, i3);
        }
    }

    fprintf(fp, "CELL_TYPES %zu\n", num_quads);
    for (size_t i = 0; i < num_quads; i++) {
        fprintf(fp, "9\n");
    }

    fclose(fp);
}