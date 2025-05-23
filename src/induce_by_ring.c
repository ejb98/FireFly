#include "vector.h"
#include "induce_by_segment.h"

Vector induce_by_ring(Vector *point, Vector *corners, Vector *v_induced, double gamma, double cutoff) {
    Vector u_induced;
    Vector vel_induced[4];
    
    induce_by_segment(point, corners, corners + 1, vel_induced, gamma, cutoff);
    induce_by_segment(point, corners + 1, corners + 2, vel_induced + 1, gamma, cutoff);
    induce_by_segment(point, corners + 2, corners + 3, vel_induced + 2, gamma, cutoff);
    induce_by_segment(point, corners + 3, corners, vel_induced + 3, gamma, cutoff);

    u_induced.x = vel_induced[1].x + vel_induced[3].x;
    u_induced.y = vel_induced[1].y + vel_induced[3].y;
    u_induced.z = vel_induced[1].z + vel_induced[3].z;

    v_induced->x = vel_induced[0].x + vel_induced[2].x + u_induced.x;
    v_induced->y = vel_induced[0].y + vel_induced[2].y + u_induced.y;
    v_induced->z = vel_induced[0].z + vel_induced[2].z + u_induced.z;

    return u_induced;
}