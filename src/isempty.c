#include "mesh.h"

int isempty(const Mesh *mesh) {
    int shape_undefined = !mesh->num_rows || !mesh->num_cols;
    int any_pointers_null = !mesh->x || !mesh->y || !mesh->z;

    return (shape_undefined || any_pointers_null);
}