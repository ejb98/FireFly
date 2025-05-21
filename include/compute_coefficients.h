#ifndef COMPUTE_COEFFICIENTS_H
#define COMPUTE_COEFFICIENTS_H

#include "mesh.h"
#include "vector.h"

void compute_coefficients(Mesh *points, Mesh *rings, Mesh *normals,
                          double *a, double *b, Vector *buffer,
                          double cutoff, int mirror, int append);
#endif