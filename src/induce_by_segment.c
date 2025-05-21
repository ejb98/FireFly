#include <math.h>

#include "dot.h"
#include "cross.h"
#include "vector.h"

void induce_by_segment(Vector *p, Vector *p1, Vector *p2,
                       Vector *velocity, double gamma, double cutoff) {
    Vector r1 = {p->x - p1->x, p->y - p1->y, p->z - p1->z};
    Vector r2 = {p->x - p2->x, p->y - p2->y, p->z - p2->z};
    Vector r1_c_r2;

    cross(&r1, &r2, &r1_c_r2);

    double mag_sq = dot(&r1_c_r2, &r1_c_r2);
    double mag_r1 = sqrt(dot(&r1, &r1));
    double mag_r2 = sqrt(dot(&r2, &r2));

    if (mag_r1 < cutoff || mag_r2 < cutoff || mag_sq < cutoff) {
        velocity->x = 0.0;
        velocity->y = 0.0;
        velocity->z = 0.0;

        return;
    }

    Vector r0 = {p2->x - p1->x, p2->y - p1->y, p2->z - p1->z};

    const double pi = 3.141592654;
    double r0_d_r1 = dot(&r0, &r1);
    double r0_d_r2 = dot(&r0, &r2);
    double c = gamma / (4.0 * pi * mag_sq) * (r0_d_r1 / mag_r1 - r0_d_r2 / mag_r2);

    velocity->x = r1_c_r2.x * c;
    velocity->y = r1_c_r2.y * c;
    velocity->z = r1_c_r2.z * c;
}