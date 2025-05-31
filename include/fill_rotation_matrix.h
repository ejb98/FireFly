#ifndef FILL_ROTATION_MATRIX_H
#define FILL_ROTATION_MATRIX_H

#include <math.h>

#include "vector3d.h"

void fill_rotation_matrix(const Vector3D *rotation, double rotation_matrix[3][3]) {
    double cg = cos(rotation->z);
    double cb = cos(rotation->y);
    double ca = cos(rotation->x);
    double sg = sin(rotation->z);
    double sb = sin(rotation->y);
    double sa = sin(rotation->x);

    rotation_matrix[0][0] = cb * cg;
    rotation_matrix[0][1] = sa * sb * cg - ca * sg;
    rotation_matrix[0][2] = ca * sb * cg + sa * sg;
    rotation_matrix[1][0] = cb * sg;
    rotation_matrix[1][1] = sa * sb * sg + ca * cg;
    rotation_matrix[1][2] = ca * sb * sg - sa * cg;
    rotation_matrix[2][0] = -sb;
    rotation_matrix[2][1] = sa * cb;
    rotation_matrix[2][2] = ca * cb;
}

#endif