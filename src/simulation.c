#include <stdbool.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "wing.h"
#include "sub2ind.h"
#include "vector3d.h"
#include "constants.h"
#include "simulation.h"
#include "allocate_doubles.h"
#include "fill_rotation_matrix.h"

Simulation *Simulation_Init(const Wing *wing, int num_time_steps, double delta_time,
                            double starting_vortex_offset, double cutoff_radius) {
    Simulation *sim = (Simulation *) malloc(sizeof(Simulation));

    if (sim == NULL) {
        fprintf(stderr, "Simulation_Init: malloc returned NULL");

        return NULL;
    }

    Vector3D zeros = {0.0, 0.0, 0.0};

    size_t num_control_points = wing->num_chordwise_panels * wing->num_spanwise_panels;
    size_t num_surface_points = (wing->num_chordwise_panels + 1) * (wing->num_spanwise_panels + 1);
    size_t num_wake_rings_max = (num_time_steps - 1) * (wing->num_spanwise_panels);
    size_t num_wake_points_max = num_time_steps * (wing->num_spanwise_panels + 1);

    sim->iteration = -1;
    sim->num_time_steps = num_time_steps;

    int *pivot_vector = (int *) malloc(sizeof(int) * num_control_points);

    if (pivot_vector == NULL) {
        fprintf(stderr, "Simulation_Init: malloc returned NULL");

        return NULL;
    }

    sim->is_complete = false;
    sim->has_updated_geometry = false;

    sim->delta_time = delta_time;
    sim->cutoff_radius = cutoff_radius;
    sim->starting_vortex_offset = starting_vortex_offset;

    sim->pressures = AllocateDoubles(num_control_points);
    sim->a_wing_on_wing = AllocateDoubles(num_control_points * num_control_points);
    sim->b_wing_on_wing = AllocateDoubles(num_control_points * num_control_points);
    sim->surface_areas = AllocateDoubles(num_control_points);
    sim->right_hand_side = AllocateDoubles(num_control_points);
    sim->wake_vortex_strengths = AllocateDoubles(num_wake_rings_max);
    sim->bound_vortex_strengths = AllocateDoubles(num_control_points);
    sim->last_bound_vortex_strengths = AllocateDoubles(num_control_points);

    sim->wing = *wing;

    sim->chordwise_velocity_buffer = zeros;

    sim->normals = Vector3D_Allocate(num_control_points);
    sim->surface_points = Vector3D_Allocate(num_surface_points);
    sim->control_points = Vector3D_Allocate(num_control_points);
    sim->wake_ring_points = Vector3D_Allocate(num_wake_points_max);
    sim->bound_ring_points = Vector3D_Allocate(num_surface_points);
    sim->spanwise_tangents = Vector3D_Allocate(num_control_points);
    sim->chordwise_tangents = Vector3D_Allocate(num_control_points);
    sim->last_control_points = Vector3D_Allocate(num_control_points);
    sim->kinematic_velocities = Vector3D_Allocate(num_control_points);
    sim->wake_induced_velocities = Vector3D_Allocate(num_control_points * num_wake_rings_max);
    sim->wake_point_displacements = Vector3D_Allocate(num_wake_points_max);
    sim->spanwise_velocity_buffer = Vector3D_Allocate(wing->num_spanwise_panels);

    return sim;
}

void Simulation_Deallocate(Simulation *sim) {    
    free(sim->pivot_vector);

    free(sim->pressures);
    free(sim->a_wing_on_wing);
    free(sim->b_wing_on_wing);
    free(sim->surface_areas);
    free(sim->right_hand_side);
    free(sim->wake_vortex_strengths);
    free(sim->bound_vortex_strengths);
    free(sim->last_bound_vortex_strengths);

    free(sim->normals);
    free(sim->surface_points);
    free(sim->control_points);
    free(sim->wake_ring_points);
    free(sim->bound_ring_points);
    free(sim->spanwise_tangents);
    free(sim->chordwise_tangents);
    free(sim->last_control_points);
    free(sim->kinematic_velocities);
    free(sim->wake_induced_velocities);
    free(sim->wake_point_displacements);
    free(sim->spanwise_velocity_buffer);

    free(sim);
}

void Simulation_ComputeSurfacePoints(Simulation *sim) {
    double constant;
    double local_chord;
    double normalized_x;
    double normalized_z;
    double leading_offset;
    double trailing_offset;
    double rotation_matrix[9];
    double p = sim->wing.naca_p / 10.0;
    double m = sim->wing.naca_m / 100.0;
    double leading_tangent = tan((90.0 - sim->wing.leading_sweep_angle) * PI / 180.0);
    double trailing_tangent = tan((90.0 - sim->wing.trailing_sweep_angle) * PI / 180.0);

    Vector3D rotation = {0.0, sim->wing.angle_of_attack * PI / 180.0, 0.0};
    Vector3D *point;

    FillRotationMatrix(&rotation, rotation_matrix);

    for (int j = 0; j < sim->wing.num_spanwise_panels + 1; j++) {
        for (int i = 0; i < sim->wing.num_chordwise_panels + 1; i++) {
            point = sim->surface_points + Sub2Ind(i, j, sim->wing.num_spanwise_panels + 1);

            point->y = sim->wing.semi_span * j / sim->wing.num_spanwise_panels;
            normalized_x = ((double) i) / sim->wing.num_chordwise_panels;

            constant = 2.0 * p * normalized_x - normalized_x * normalized_x;

            if (normalized_x >= p || !sim->wing.naca_p) {
                normalized_z = (m / ((1.0 - p) * (1.0 - p))) * (1.0 - 2.0 * p + constant);
            } else {
                normalized_z = (m / (p * p)) * constant;
            }

            leading_offset = point->y * leading_tangent;
            trailing_offset = point->y * trailing_tangent;

            local_chord = sim->wing.root_chord + trailing_offset - leading_offset;

            point->x = leading_offset + normalized_x * local_chord;
            point->z = local_chord * normalized_z;

            Vector3D_Rotate(point, rotation_matrix);
        }
    }
}

void Simulation_ComputeSurfaceVectors(Simulation *sim) {
    size_t ipanel;
    Vector3D left;
    Vector3D back;
    Vector3D front;
    Vector3D right;
    Vector3D vectora;
    Vector3D vectorb;
    Vector3D *normal;
    Vector3D *corners[4];

    for (int i = 0; i < sim->wing.num_chordwise_panels; i++) {
        for (int j = 0; j < sim->wing.num_spanwise_panels; j++) {
            ipanel = Sub2Ind(i, j, sim->wing.num_spanwise_panels);
            normal = sim->normals + ipanel;

            Simulation_GetCorners(sim, SURFACE_POINTS, i, j, corners);
            Vector3D_Lerp(corners[0], corners[3], 0.5, &left);
            Vector3D_Lerp(corners[2], corners[3], 0.5, &back);
            Vector3D_Lerp(corners[0], corners[1], 0.5, &front);
            Vector3D_Lerp(corners[1], corners[2], 0.5, &right);
            Vector3D_Subtract(corners[2], corners[0], &vectora);
            Vector3D_Subtract(corners[1], corners[3], &vectorb);
            Vector3D_Cross(&vectora, &vectorb, normal);
            Vector3D_Normalize(normal);
            Vector3D_Direction(&front, &back, sim->chordwise_tangents + ipanel);
            Vector3D_Direction(&left, &right, sim->spanwise_tangents + ipanel);
        }
    }
}

Vector3D *Simulation_GetPoints(const Simulation *sim, Geometry geometry) {
    Vector3D *points;

    switch (geometry) {
        case CONTROL_POINTS:
            points = sim->control_points;

            break;
        case SURFACE_POINTS:
            points = sim->surface_points;

            break;
        case WAKE_RING_POINTS:
            points = sim->wake_ring_points;

            break;
        case BOUND_RING_POINTS:
            points = sim->bound_ring_points;
    }

    return points;
}

int Simulation_GetNumRows(const Simulation* sim, Geometry geometry) {
    switch (geometry) {
        case CONTROL_POINTS:
            return sim->wing.num_chordwise_panels;

        case SURFACE_POINTS:
        case BOUND_RING_POINTS:
            return sim->wing.num_chordwise_panels + 1;

        case WAKE_RING_POINTS:
            return sim->iteration + 1;
    }
}

int Simulation_GetNumColumns(const Simulation* sim, Geometry geometry) {
    switch (geometry) {
        case SURFACE_POINTS:
        case BOUND_RING_POINTS:
            return sim->wing.num_spanwise_panels + 1;

        case CONTROL_POINTS:
            return sim->wing.num_spanwise_panels;

        case WAKE_RING_POINTS:
            return sim->iteration < 0 ? 0 : sim->wing.num_spanwise_panels + 1;
    }
}

size_t Simulation_GetNumPoints(const Simulation* sim, Geometry geometry) {
    size_t num_rows = Simulation_GetNumRows(sim, geometry);
    size_t num_cols = Simulation_GetNumColumns(sim, geometry);

    return num_rows * num_cols;
}

size_t Simulation_GetNumQuads(const Simulation* sim, Geometry geometry) {
    size_t num_rows = Simulation_GetNumRows(sim, geometry);
    size_t num_cols = Simulation_GetNumColumns(sim, geometry);

    return (num_rows - 1) * (num_cols - 1);
}

void Simulation_GetCorners(const Simulation *sim, Geometry geometry, int i, int j, Vector3D *corners[]) {
    int num_cols = Simulation_GetNumColumns(sim, geometry);

    Vector3D *points = Simulation_GetPoints(sim, geometry);

    corners[0] = points + Sub2Ind(i, j, num_cols);
    corners[1] = points + Sub2Ind(i, j + 1, num_cols);
    corners[2] = points + Sub2Ind(i + 1, j + 1, num_cols);
    corners[3] = points + Sub2Ind(i + 1, j, num_cols);
}

void Simulation_WritePoints2VTK(const Simulation *sim, Geometry geometry, const char *file_path) {
    int num_rows = Simulation_GetNumRows(sim, geometry);
    int num_cols = Simulation_GetNumColumns(sim, geometry);

    size_t num_quads = Simulation_GetNumQuads(sim, geometry);
    size_t num_points = Simulation_GetNumPoints(sim, geometry);

    char full_path[175];
    char description[20];

    Vector3D *points = Simulation_GetPoints(sim, geometry);

    switch (geometry) {
        case CONTROL_POINTS:
            strcpy(description, "control_points");

            break;
        case SURFACE_POINTS:
            strcpy(description, "surface_points");

            break;
        case WAKE_RING_POINTS:
            strcpy(description, "wake_ring_points");

            break;
        case BOUND_RING_POINTS:
            strcpy(description, "bound_ring_points");
    }
    
    if (num_rows < 2) {
        fprintf(stderr, "Simulation_WritePoints2VTK: points must have at least two rows");

        return;
    }

    if (num_cols < 2) {
        fprintf(stderr, "Simulation_WritePoints2VTK: points must have at least two columns");

        return;
    }

    if (sim->iteration < 0) {
        snprintf(full_path, sizeof(full_path), "%s%s.vtk", file_path, description);
    } else {
        snprintf(full_path, sizeof(full_path), "%s%s.vtk.%d", file_path, description, sim->iteration);
    }

    FILE *file = fopen(full_path, "w");

    if (file == NULL) {
        fprintf(stderr, "Simulation_WritePoints2VTK: failed to open %s", full_path);

        return;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Mesh Surface\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(file, "POINTS %zu float\n", num_points);
    for (size_t i = 0; i < num_points; i++) {
        fprintf(file, "%f %f %f\n", points[i].x, points[i].y, points[i].z);
    }

    fprintf(file, "CELLS %zu %zu\n", num_quads, 5 * num_quads);
    for (int i = 0; i < num_rows - 1; i++) {
        for (int j = 0; j < num_cols - 1; j++) {
            fprintf(file, "4 %zu %zu %zu %zu\n", Sub2Ind(i, j, num_cols),
                                                 Sub2Ind(i, j + 1, num_cols),
                                                 Sub2Ind(i + 1, j + 1, num_cols),
                                                 Sub2Ind(i + 1, j, num_cols));
        }
    }

    fprintf(file, "CELL_TYPES %zu\n", num_quads);
    for (size_t i = 0; i < num_quads; i++) {
        fprintf(file, "9\n");
    }

    fclose(file);
}

void Simulation_ComputeControlPoints(Simulation *sim) {
    size_t ipanel;
    Vector3D left;
    Vector3D right;
    Vector3D *corners[4];

    for (int i = 0; i < sim->wing.num_chordwise_panels; i++) {
        for (int j = 0; j < sim->wing.num_spanwise_panels; j++) {
            ipanel = Sub2Ind(i, j, sim->wing.num_spanwise_panels);

            Simulation_GetCorners(sim, SURFACE_POINTS, i, j, corners);
            Vector3D_Lerp(corners[0], corners[3], 0.75, &left);
            Vector3D_Lerp(corners[1], corners[2], 0.75, &right);
            Vector3D_Lerp(&left, &right, 0.5, sim->control_points + ipanel);
        }
    }
}

void Simulation_ComputeSurfaceAreas(Simulation *sim) {
    size_t ipanel;
    Vector3D *corners[4];

    double a1, a2, b1, b2;

    for (int i = 0; i < sim->wing.num_chordwise_panels; i++) {
        for (int j = 0; j < sim->wing.num_spanwise_panels; j++) {
            ipanel = Sub2Ind(i, j, sim->wing.num_spanwise_panels);
            Simulation_GetCorners(sim, SURFACE_POINTS, i, j, corners);

            a1 = Vector3D_Distance(corners[0], corners[1]);
            b1 = Vector3D_Distance(corners[1], corners[2]);
            a2 = Vector3D_Distance(corners[2], corners[3]);
            b2 = Vector3D_Distance(corners[3], corners[0]);

            sim->surface_areas[ipanel] = 0.5 * (a1 * b1 + a2 * b2);
        }
    }
}

void Simulation_ComputeBoundRingPoints(Simulation *sim) {
    size_t ipoint;

    Vector3D *next;
    Vector3D *point;

    int num_rows = sim->wing.num_chordwise_panels + 1;
    int num_cols = sim->wing.num_spanwise_panels + 1;

    for (int j = 0; j < num_cols; j++) {
        for (int i = 0; i < num_rows; i++) {
            ipoint = Sub2Ind(i, j, num_cols);
            point = sim->surface_points + ipoint;

            if (i == num_rows - 1) {
                sim->bound_ring_points[ipoint].x = point->x + sim->starting_vortex_offset;
                sim->bound_ring_points[ipoint].y = point->y;
                sim->bound_ring_points[ipoint].z = point->z;
            } else {
                next = sim->surface_points + Sub2Ind(i + 1, j, num_cols);
                Vector3D_Lerp(point, next, 0.25, sim->bound_ring_points + ipoint);
            }
        }
    }
}

void Simulation_Process(Simulation *sim) {
    sim->iteration++;

    printf("Solving step %d...", sim->iteration);

    Simulation_ComputeKinematicVelocities(sim);

    /*
    * Shed Wake
    * Solve System of Equations
    * Roll up Wake if 0 < Iteration < Steps - 1
    */

    printf("done\n");

    if (sim->iteration == sim->num_time_steps - 1) {
        sim->is_complete = true;

        return;
    } else {
        sim->wing.last_position = sim->wing.position;
        sim->wing.last_rotation = sim->wing.rotation;
    }
}

void Simulation_ComputeKinematicVelocities(Simulation *sim) {
    size_t num_control_points = Simulation_GetNumPoints(sim, CONTROL_POINTS);
    
    if (sim->iteration) {
        Vector3D zeros = {0.0, 0.0, 0.0};

        for (size_t i = 0; i < num_control_points; i++) {
            sim->kinematic_velocities[i] = zeros;
            sim->last_control_points[i] = sim->control_points[i];
        }
    } else {
        Vector3D current;
        Vector3D previous;
        Vector3D inverse_rotation = sim->wing.rotation;

        Vector3D *velocity;

        Vector3D_Multiply(&inverse_rotation, -1.0);
        
        double rotation_matrix[9];
        double last_rotation_matrix[9];
        double inverse_rotation_matrix[9];

        FillRotationMatrix(&sim->wing.rotation, rotation_matrix);
        FillRotationMatrix(&sim->wing.last_rotation, last_rotation_matrix);
        FillRotationMatrix(&inverse_rotation, inverse_rotation_matrix);

        for (size_t i = 0; i < num_control_points; i++) {
            current = sim->control_points[i];
            previous = sim->last_control_points[i];

            velocity = sim->kinematic_velocities + i;

            Vector3D_Rotate(&current, rotation_matrix);
            Vector3D_Rotate(&previous, last_rotation_matrix);
            Vector3D_Add(&current, &sim->wing.position, &current);
            Vector3D_Add(&previous, &sim->wing.last_position, &previous);
            Vector3D_Subtract(&current, &previous, velocity);
            Vector3D_Divide(velocity, sim->delta_time);

            sim->last_control_points[i] = sim->control_points[i];
        }
    }
}