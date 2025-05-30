#ifndef MESH_TO_VECTOR_H
#define MESH_TO_VECTOR_H

#include <stdlib.h>

#include "mesh.h"
#include "vector.h"

void mesh_to_vector(const Mesh *mesh, size_t index, Vector3D *vector);

#endif