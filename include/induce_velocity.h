#ifndef INDUCE_VELOCITY_H
#define INDUCE_VELOCITY_H

#include "mesh.h"
#include "vector.h"

void induce_velocity(Mesh *points,
                     Mesh *rings,
                     Mesh *normals_or_velocities,
                     double *a,
                     double *b_or_gammas,
                     Vector *velocity_buffer,
                     double cutoff,
                     int mirror,
                     int append);
#endif