#include <stdlib.h>

#include "mesh.h"

int isempty(const Mesh *mesh) {
    int shape_undefined = (mesh->num_rows == 0) || (mesh->num_cols == 0);
    int any_pointers_null = (mesh->x == NULL) || (mesh->y == NULL) || (mesh->z == NULL);

    return (shape_undefined || any_pointers_null);
}