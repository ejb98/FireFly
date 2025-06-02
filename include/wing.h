#ifndef WING_H
#define WING_H

#include "vector3d.h"
#include "geometry.h"

typedef struct Wing {
    int naca_m;
    int naca_p;
    int num_spanwise_panels;
    int num_chordwise_panels;
    double semi_span;
    double root_chord;
    double angle_of_attack;
    double leading_sweep_angle;
    double trailing_sweep_angle;

    Vector3D position;
    Vector3D rotation;
    Vector3D last_position;
    Vector3D last_rotation;
} Wing;

void Wing_Print(const Wing *wing);

#endif