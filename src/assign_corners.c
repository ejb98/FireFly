#include "mesh.h"
#include "vector.h"
#include "sub2ind.h"
#include "mesh_to_vector.h"

void assign_corners(const Mesh *mesh, int i, int j, Vector *corners) {
    mesh_to_vector(mesh, sub2ind(i, j, mesh->num_cols), corners);
    mesh_to_vector(mesh, sub2ind(i, j + 1, mesh->num_cols), corners + 1);
    mesh_to_vector(mesh, sub2ind(i + 1, j + 1, mesh->num_cols), corners + 2);
    mesh_to_vector(mesh, sub2ind(i + 1, j, mesh->num_cols), corners + 3);
}