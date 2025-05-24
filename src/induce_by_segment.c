#include <math.h>

#include "vector.h"

#define FOUR_PI 12.566370614

void induce_by_segment(Vector *point, Vector *point1, Vector *point2,
                       Vector *velocity, double gamma, double cutoff) {
    Vector r1 = {point->x - point1->x, point->y - point1->y, point->z - point1->z};
    Vector r2 = {point->x - point2->x, point->y - point2->y, point->z - point2->z};
    Vector r1_c_r2 = {r1.y * r2.z - r1.z * r2.y, r1.z * r2.x - r1.x * r2.z, r1.x * r2.y - r1.y * r2.x};

    double mag_sq = r1_c_r2.x * r1_c_r2.x + r1_c_r2.y * r1_c_r2.y + r1_c_r2.z * r1_c_r2.z;
    double mag_r1 = sqrt(r1.x * r1.x + r1.y * r1.y + r1.z * r1.z);
    double mag_r2 = sqrt(r2.x * r2.x + r2.y * r2.y + r2.z * r2.z);

    if (mag_r1 < cutoff || mag_r2 < cutoff || mag_sq < cutoff) {
        velocity->x = 0.0;
        velocity->y = 0.0;
        velocity->z = 0.0;

        return;
    }

    Vector r0 = {point2->x - point1->x, point2->y - point1->y, point2->z - point1->z};

    double r0_d_r1 = r0.x * r1.x + r0.y * r1.y + r0.z * r1.z;
    double r0_d_r2 = r0.x * r2.x + r0.y * r2.y + r0.z * r2.z;
    double constant = gamma / (FOUR_PI * mag_sq) * (r0_d_r1 / mag_r1 - r0_d_r2 / mag_r2);

    velocity->x = r1_c_r2.x * constant;
    velocity->y = r1_c_r2.y * constant;
    velocity->z = r1_c_r2.z * constant;
}