#ifndef WING_H
#define WING_H

#include "vector3d.h"
#include "geometry.h"

typedef struct Wing {
    int naca_m;
    int naca_p;
    int iteration;
    int num_time_steps;
    int geometry_changed;
    int num_spanwise_panels;
    int num_chordwise_panels;

    int *pivot_vector;

    size_t num_control_points;

    double semi_span;
    double root_chord;
    double surface_area;
    double cutoff_radius;
    double angle_of_attack;
    double starting_vortex_offset;
    double leading_edge_sweep_angle;
    double trailing_edge_sweep_angle;

    double *pressures;
    double *a_wing_on_wing;
    double *b_wing_on_wing;
    double *surface_areas;
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
    Vector3D *previous_control_points;
    Vector3D *wake_induced_velocities;
    Vector3D *wake_point_displacements;
    Vector3D *spanwise_tangent_vectors;
    Vector3D *chordwise_tangent_vectors;
    Vector3D *horizontal_velocity_buffer;
} Wing;

typedef struct WingProperties {
    int naca_m;
    int naca_p;
    int num_spanwise_panels;
    int num_chordwise_panels;
    double semi_span;
    double root_chord;
    double angle_of_attack;
    double leading_edge_sweep_angle;
    double trailing_edge_sweep_angle;
} WingProperties;

void Wing_Deallocate(Wing *wing);
void Wing_ComputeSurfacePoints(Wing *wing);
void Wing_ComputeSurfaceVectors(Wing *wing);
void Wing_ComputeSurfaceAreas(Wing *wing);
void Wing_ComputeKinematicVelocities(Wing *wing, double delta_time);
void Wing_ComputeBoundRingPoints(Wing *wing);
void Wing_ComputeControlPoints(Wing *wing);
void Wing_PrintAttributes(const Wing *wing);
void Wing_WritePoints2VTK(const Wing *wing, Geometry geometry, const char *file_path);
void Wing_GetCorners(const Wing *wing, Geometry geometry, int i, int j, Vector3D **corners);
void Wing_Process(Wing *wing, double delta_time);
Vector3D *Wing_GetPoints(const Wing *wing, Geometry geometry, int *num_rows, int *num_cols);
Wing *Wing_Construct(const WingProperties *wing_properties, int num_time_steps,
                     double starting_vortex_offset, double cutoff_radius);
                     
#endif