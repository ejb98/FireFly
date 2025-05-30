#include <stdlib.h>

#include "mesh.h"
#include "vector.h"

void mesh_to_vector(const Mesh *mesh, size_t index, Vector3D *vector) {
    vector->x = mesh->x[index];
    vector->y = mesh->y[index];
    vector->z = mesh->z[index];
}