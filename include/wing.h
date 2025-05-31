#ifndef WING_H
#define WING_H

#include "vector3d.h"
#include "geometry.h"

typedef struct Wing {
    int naca_m;
    int naca_p;
    int iteration;
    int num_time_steps;
    int num_spanwise_panels;
    int num_chordwise_panels;

    int *pivot_vector;

    double semi_span;
    double root_chord;
    double cutoff_radius;
    double angle_of_attack;
    double starting_vortex_offset;
    double leading_edge_sweep_angle;
    double trailing_edge_sweep_angle;

    double *pressures;
    double *a_wing_on_wing;
    double *b_wing_on_wing;
    double *right_hand_side;
    double *wake_vortex_strengths;
    double *bound_vortex_strengths;
    double *previous_bound_vortex_strengths;

    Vector3D position;
    Vector3D rotation;
    Vector3D previous_position;
    Vector3D previous_rotation;
    Vector3D vertical_velocity_buffer;

    Vector3D *normal_vectors;
    Vector3D *surface_points;
    Vector3D *control_points;
    Vector3D *wake_ring_points;
    Vector3D *bound_ring_points;
    Vector3D *kinematic_velocities;
    Vector3D *wake_induced_velocities;
    Vector3D *wake_point_displacements;
    Vector3D *spanwise_tangent_vectors;
    Vector3D *chordwise_tangent_vectors;
    Vector3D *horizontal_velocity_buffer;
} Wing;

Wing *Wing_construct(int naca_m,
                     int naca_p,
                     int num_time_steps,
                     int num_spanwise_panels,
                     int num_chordwise_panels,
                     double semi_span,
                     double root_chord,
                     double cutoff_radius,
                     double angle_of_attack,
                     double starting_vortex_offset,
                     double leading_edge_sweep_angle,
                     double trailing_edge_sweep_angle);

void Wing_free(Wing *wing);
void Wing_compute_surface_points(Wing *wing);
void Wing_compute_surface_vectors(Wing *wing);
void Wing_print_attributes(const Wing *wing);
void Wing_write_points_to_vtk(const Wing *wing, Geometry geometry, const char *file_path);
void Wing_get_corners(const Wing *wing, Geometry geometry, int i, int j, Vector3D **corners);
Vector3D *Wing_get_points(const Wing *wing, Geometry geometry, int *num_rows, int *num_cols);

#endif