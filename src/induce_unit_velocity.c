#include <math.h>

#include "vector3d.h"
#include "constants.h"

void InduceUnitVelocity(const Vector3D *point, const Vector3D *point1, 
                        const Vector3D *point2, Vector3D *unit_velocity, double cutoff) {
    Vector3D r1 = {point->x - point1->x, point->y - point1->y, point->z - point1->z};
    Vector3D r2 = {point->x - point2->x, point->y - point2->y, point->z - point2->z};
    Vector3D r1_c_r2 = {r1.y * r2.z - r1.z * r2.y, r1.z * r2.x - r1.x * r2.z, r1.x * r2.y - r1.y * r2.x};

    double mag_sq = r1_c_r2.x * r1_c_r2.x + r1_c_r2.y * r1_c_r2.y + r1_c_r2.z * r1_c_r2.z;
    double mag_r1 = sqrt(r1.x * r1.x + r1.y * r1.y + r1.z * r1.z);
    double mag_r2 = sqrt(r2.x * r2.x + r2.y * r2.y + r2.z * r2.z);

    if (mag_r1 < cutoff || mag_r2 < cutoff || mag_sq < cutoff) {
        unit_velocity->x = 0.0;
        unit_velocity->y = 0.0;
        unit_velocity->z = 0.0;

        return;
    }

    Vector3D r0 = {point2->x - point1->x, point2->y - point1->y, point2->z - point1->z};

    double r0_d_r1 = r0.x * r1.x + r0.y * r1.y + r0.z * r1.z;
    double r0_d_r2 = r0.x * r2.x + r0.y * r2.y + r0.z * r2.z;
    double constant = (r0_d_r1 / mag_r1 - r0_d_r2 / mag_r2) / (FOUR_PI * mag_sq);

    unit_velocity->x = r1_c_r2.x * constant;
    unit_velocity->y = r1_c_r2.y * constant;
    unit_velocity->z = r1_c_r2.z * constant;
}