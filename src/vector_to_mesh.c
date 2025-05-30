#include <stdlib.h>

#include "mesh.h"
#include "vector.h"

void vector_to_mesh(const Vector3D *vector, Mesh *mesh, size_t index) {
    mesh->x[index] = vector->x;
    mesh->y[index] = vector->y;
    mesh->z[index] = vector->z;
}