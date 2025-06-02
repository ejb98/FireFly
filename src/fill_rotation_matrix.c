#include <math.h>

#include "vector3d.h"

void FillRotationMatrix(const Vector3D *rotation, double *rotation_matrix) {
    double cg = cos(rotation->z);
    double cb = cos(rotation->y);
    double ca = cos(rotation->x);
    double sg = sin(rotation->z);
    double sb = sin(rotation->y);
    double sa = sin(rotation->x);

    rotation_matrix[0] = cb * cg;
    rotation_matrix[1] = sa * sb * cg - ca * sg;
    rotation_matrix[2] = ca * sb * cg + sa * sg;
    rotation_matrix[3] = cb * sg;
    rotation_matrix[4] = sa * sb * sg + ca * cg;
    rotation_matrix[5] = ca * sb * sg - sa * cg;
    rotation_matrix[6] = -sb;
    rotation_matrix[7] = sa * cb;
    rotation_matrix[8] = ca * cb;
}