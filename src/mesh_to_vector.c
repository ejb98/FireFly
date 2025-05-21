#include "mesh.h"
#include "vector.h"

void mesh_to_vector(const Mesh *mesh, int index, Vector *v) {
    v->x = mesh->x[index];
    v->y = mesh->y[index];
    v->z = mesh->z[index];
}