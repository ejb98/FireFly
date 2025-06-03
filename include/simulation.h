#ifndef SIMULATION_H
#define SIMULATION_H

#include <stdbool.h>

#include "wing.h"
#include "vector3d.h"

typedef struct Simulation {
    int iteration;
    int num_time_steps;

    int *pivot_vector;

    bool is_complete;
    bool has_updated_geometry;

    double delta_time;
    double cutoff_radius;
    double starting_vortex_offset;

    double *pressures;
    double *a_wing_on_wing;
    double *b_wing_on_wing;
    double *surface_areas;
    double *right_hand_side;
    double *wake_vortex_strengths;
    double *bound_vortex_strengths;
    double *last_bound_vortex_strengths;

    Wing *wing;

    Vector3D unit_velocities[4];
    Vector3D chordwise_velocity_buffer;

    Vector3D *normals;
    Vector3D *surface_points;
    Vector3D *control_points;
    Vector3D *corner_buffer[4];
    Vector3D *wake_ring_points;
    Vector3D *bound_ring_points;
    Vector3D *spanwise_tangents;
    Vector3D *chordwise_tangents;
    Vector3D *last_control_points;
    Vector3D *kinematic_velocities;
    Vector3D *wake_induced_velocities;
    Vector3D *wake_point_displacements;
    Vector3D *spanwise_velocity_buffer;
} Simulation;

void Simulation_ComputeSurfacePoints(Simulation *sim);
void Simulation_ComputeSurfaceVectors(Simulation *sim);
void Simulation_ComputeSurfaceAreas(Simulation *sim);
void Simulation_ComputeKinematicVelocities(Simulation *sim);
void Simulation_ComputeBoundRingPoints(Simulation *sim);
void Simulation_ComputeControlPoints(Simulation *sim);
void Simulation_ComputeCoefficients(Simulation *sim);
void Simulation_ComputeWakeInducedVelocities(Simulation *sim);
void Simulation_RollupWake(Simulation *sim);
void Simulation_ShedWake(Simulation *sim);
void Simulation_Process(Simulation *sim);
void Simulation_Solve(Simulation *sim);
void Simulation_Deallocate(Simulation *sim);
void Simulation_InduceUnitVelocities(Simulation *sim, Geometry geometry, int i, int j, const Vector3D *point);
void Simulation_WritePoints2VTK(const Simulation *sim, Geometry geometry, const char *file_path);
void Simulation_GetCorners(const Simulation *sim, Geometry geometry, int i, int j, Vector3D *corners[]);
int Simulation_GetNumRows(const Simulation* sim, Geometry geometry);
int Simulation_GetNumColumns(const Simulation* sim, Geometry geometry);
size_t Simulation_GetNumPoints(const Simulation* sim, Geometry geometry);
size_t Simulation_GetNumQuads(const Simulation* sim, Geometry geometry);
Vector3D *Simulation_GetPoints(const Simulation *sim, Geometry geometry);
Simulation *Simulation_Init(Wing *wing, int num_time_steps, double delta_time,
                            double starting_vortex_offset, double cutoff_radius);

#endif