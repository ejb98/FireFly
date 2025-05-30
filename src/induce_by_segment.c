#include <math.h>

#include "vector.h"

#define FOUR_PI 12.566370614

void induce_by_segment(Vector3D *point, Vector3D *point1, Vector3D *point2,
                       Vector3D *velocity_normalized, double cutoff) {
    Vector3D r1 = {point->x - point1->x, point->y - point1->y, point->z - point1->z};
    Vector3D r2 = {point->x - point2->x, point->y - point2->y, point->z - point2->z};
    Vector3D r1_c_r2 = {r1.y * r2.z - r1.z * r2.y, r1.z * r2.x - r1.x * r2.z, r1.x * r2.y - r1.y * r2.x};

    double mag_sq = r1_c_r2.x * r1_c_r2.x + r1_c_r2.y * r1_c_r2.y + r1_c_r2.z * r1_c_r2.z;
    double mag_r1 = sqrt(r1.x * r1.x + r1.y * r1.y + r1.z * r1.z);
    double mag_r2 = sqrt(r2.x * r2.x + r2.y * r2.y + r2.z * r2.z);

    if (mag_r1 < cutoff || mag_r2 < cutoff || mag_sq < cutoff) {
        velocity_normalized->x = 0.0;
        velocity_normalized->y = 0.0;
        velocity_normalized->z = 0.0;

        return;
    }

    Vector3D r0 = {point2->x - point1->x, point2->y - point1->y, point2->z - point1->z};

    double r0_d_r1 = r0.x * r1.x + r0.y * r1.y + r0.z * r1.z;
    double r0_d_r2 = r0.x * r2.x + r0.y * r2.y + r0.z * r2.z;
    double constant = (r0_d_r1 / mag_r1 - r0_d_r2 / mag_r2) / (FOUR_PI * mag_sq);

    velocity_normalized->x = r1_c_r2.x * constant;
    velocity_normalized->y = r1_c_r2.y * constant;
    velocity_normalized->z = r1_c_r2.z * constant;
}