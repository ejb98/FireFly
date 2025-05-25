#ifndef WING_H
#define WING_H

#include "mesh.h"
#include "arena.h"
#include "vector.h"

typedef struct Wing {
    int num_spanwise_panels;
    int num_chordwise_panels;
    int num_wake_point_rows;
    int num_wake_deforming_rows;
    int iteration;
    int naca_m;
    int naca_p;

    int *pivot_vector;

    double cutoff;
    double semi_span;
    double root_chord;
    double angle_of_attack;
    double sweep_angle_leading;
    double sweep_angle_trailing;
    double wake_offset;
    double surface_area;
    double aspect_ratio;

    double *wake_vorticity;
    double *bound_vorticity;
    double *freestream_velocities;
    double *normal_velocities;
    double *b_wing_on_wing;
    double *b_wake_on_wing;
    double *a_wing_on_wing;
    double *a_wake_on_wing;
    double *right_hand_side;

    Mesh tangent_vectors;
    Mesh normal_vectors;
    Mesh surface_panels;
    Mesh control_points;
    Mesh bound_rings;
    Mesh wake_rings;
    Mesh wake_displacements;

    Vector position;
    Vector rotation;
    Vector position_prev;
    Vector rotation_prev;

    Vector *horizontal_buffer;

    Arena memory;
} Wing;

#endif