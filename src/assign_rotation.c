#include <math.h>

void assign_rotation(double rotation_matrix[3][3],
                     double pitch, double roll, double yaw) {
    double cg = cos(yaw);
    double cb = cos(pitch);
    double ca = cos(roll);
    double sg = sin(yaw);
    double sb = sin(pitch);
    double sa = sin(roll);

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