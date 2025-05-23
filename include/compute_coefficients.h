#ifndef COMPUTE_COEFFICIENTS_H
#define COMPUTE_COEFFICIENTS_H

#include "mesh.h"

void compute_coefficients(Mesh *points, Mesh *normals, Mesh *rings,
                          double *a, double *b, double cutoff);

#endif