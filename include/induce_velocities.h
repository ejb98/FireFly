#ifndef INDUCE_VELOCITIES_H
#define INDUCE_VELOCITIES_H

#include "mesh.h"
#include "vector.h"

void induce_velocities(Vector3D *point, Mesh *rings, int i, int j, Vector3D *vertical_buffer,
                       Vector3D *horizontal_buffer, Vector3D *velocities, double cutoff);

#endif