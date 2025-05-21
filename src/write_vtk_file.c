#include <stdio.h>

#include "mesh.h"
#include "isempty.h"
#include "sub2ind.h"

void write_vtk_file(const Mesh *mesh, const char *file_name) {
    FILE *fp = fopen(file_name, "w");

    if (fp == NULL) {
        perror("failed to open file");

        return;
    }

    if (isempty(mesh)) {
        perror("write_to_vtk_file: unable to write empty mesh to file");
        fclose(fp);

        return;
    }

    int num_rows = mesh->num_rows;
    if (num_rows < 2) {
        perror("write_to_vtk_file: mesh must have at least 2 rows");
        fclose(fp);

        return;
    }

    int num_cols = mesh->num_cols;
    if (num_cols < 2) {
        perror("write_to_vtk_file: mesh must have at least 2 columns");
        fclose(fp);

        return;
    }

    int i0, i1, i2, i3;
    int num_quads = (num_rows - 1) * (num_cols - 1);
    int num_points = num_rows * num_cols;

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Mesh Surface\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %d float\n", num_points);

    for (int i = 0; i < num_points; i++) {
        fprintf(fp, "%f %f %f\n", mesh->x[i], mesh->y[i], mesh->z[i]);
    }

    fprintf(fp, "CELLS %d %d\n", num_quads, 5 * num_quads);
    for (int i = 0; i < num_rows - 1; i++) {
        for (int j = 0; j < num_cols - 1; j++) {
            i0 = sub2ind(i, j, num_cols);
            i1 = sub2ind(i, j + 1, num_cols);
            i2 = sub2ind(i + 1, j + 1, num_cols);
            i3 = sub2ind(i + 1, j, num_cols);

            fprintf(fp, "4 %d %d %d %d\n", i0, i1, i2, i3);
        }
    }

    fprintf(fp, "CELL_TYPES %d\n", num_quads);
    for (int i = 0; i < num_quads; i++) {
        fprintf(fp, "9\n");
    }

    fclose(fp);
}