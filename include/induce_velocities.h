#ifndef INDUCE_VELOCITIES_H
#define INDUCE_VELOCITIES_H

#include "mesh.h"
#include "vector.h"

void induce_velocities(Vector *point, Mesh *rings, int i, int j, Vector *vertical_buffer,
                       Vector *horizontal_buffer, Vector *velocities, double cutoff);

#endif