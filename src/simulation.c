#include <stdbool.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "wing.h"
#include "dgesv.h"
#include "sub2ind.h"
#include "vector3d.h"
#include "constants.h"
#include "simulation.h"
#include "allocate_doubles.h"
#include "fill_rotation_matrix.h"
#include "induce_unit_velocity.h"

Simulation *Simulation_Init(Wing *wing, int num_time_steps, double delta_time,
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

    sim->pivot_vector = (int *) malloc(sizeof(int) * num_control_points);

    if (sim->pivot_vector == NULL) {
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

    sim->wing = wing;

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
    sim->spanwise_velocity_buffer = Vector3D_Allocate(wing->num_spanwise_panels);
    sim->wake_point_displacements = calloc(num_wake_points_max, sizeof(Vector3D));

    if (sim->wake_point_displacements == NULL) {
        fprintf(stderr, "Simulation_Init: calloc returned NULL");

        return NULL;
    }

    Simulation_ComputeSurfacePoints(sim);
    Simulation_ComputeSurfaceAreas(sim);
    Simulation_ComputeSurfaceVectors(sim);
    Simulation_ComputeBoundRingPoints(sim);
    Simulation_ComputeControlPoints(sim);
    Simulation_ComputeCoefficients(sim);
    
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
    free(sim->spanwise_velocity_buffer);
    free(sim->wake_point_displacements);

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
    double p = sim->wing->naca_p / 10.0;
    double m = sim->wing->naca_m / 100.0;
    double leading_tangent = tan((90.0 - sim->wing->leading_sweep_angle) * PI / 180.0);
    double trailing_tangent = tan((90.0 - sim->wing->trailing_sweep_angle) * PI / 180.0);

    Vector3D rotation = {0.0, sim->wing->angle_of_attack * PI / 180.0, 0.0};
    Vector3D *point;

    FillRotationMatrix(&rotation, rotation_matrix);

    for (int j = 0; j < sim->wing->num_spanwise_panels + 1; j++) {
        for (int i = 0; i < sim->wing->num_chordwise_panels + 1; i++) {
            point = sim->surface_points + Sub2Ind(i, j, sim->wing->num_spanwise_panels + 1);

            point->y = sim->wing->semi_span * j / sim->wing->num_spanwise_panels;
            normalized_x = ((double) i) / sim->wing->num_chordwise_panels;

            constant = 2.0 * p * normalized_x - normalized_x * normalized_x;

            if (normalized_x >= p || !sim->wing->naca_p) {
                normalized_z = (m / ((1.0 - p) * (1.0 - p))) * (1.0 - 2.0 * p + constant);
            } else {
                normalized_z = (m / (p * p)) * constant;
            }

            leading_offset = point->y * leading_tangent;
            trailing_offset = point->y * trailing_tangent;

            local_chord = sim->wing->root_chord + trailing_offset - leading_offset;

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

    for (int i = 0; i < sim->wing->num_chordwise_panels; i++) {
        for (int j = 0; j < sim->wing->num_spanwise_panels; j++) {
            ipanel = Sub2Ind(i, j, sim->wing->num_spanwise_panels);
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
            return sim->wing->num_chordwise_panels;

        case SURFACE_POINTS:
        case BOUND_RING_POINTS:
            return sim->wing->num_chordwise_panels + 1;

        case WAKE_RING_POINTS:
            return sim->iteration + 1;
    }
}

int Simulation_GetNumColumns(const Simulation* sim, Geometry geometry) {
    switch (geometry) {
        case SURFACE_POINTS:
        case BOUND_RING_POINTS:
            return sim->wing->num_spanwise_panels + 1;

        case CONTROL_POINTS:
            return sim->wing->num_spanwise_panels;

        case WAKE_RING_POINTS:
            return sim->iteration < 0 ? 0 : sim->wing->num_spanwise_panels + 1;
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

    for (int i = 0; i < sim->wing->num_chordwise_panels; i++) {
        for (int j = 0; j < sim->wing->num_spanwise_panels; j++) {
            ipanel = Sub2Ind(i, j, sim->wing->num_spanwise_panels);

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

    for (int i = 0; i < sim->wing->num_chordwise_panels; i++) {
        for (int j = 0; j < sim->wing->num_spanwise_panels; j++) {
            ipanel = Sub2Ind(i, j, sim->wing->num_spanwise_panels);
            Simulation_GetCorners(sim, SURFACE_POINTS, i, j, corners);

            // Compute the area of two triangles and add together
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

    int num_rows = sim->wing->num_chordwise_panels + 1;
    int num_cols = sim->wing->num_spanwise_panels + 1;

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
    Simulation_ShedWake(sim);

    if (sim->iteration) {
        Simulation_ComputeWakeInducedVelocities(sim);
        Simulation_Solve(sim);
    }
    
    if (sim->iteration > 0 && sim->iteration < sim->num_time_steps - 1) {
        Simulation_RollupWake(sim);
    }

    printf("done\n");

    if (sim->iteration == sim->num_time_steps - 1) {
        sim->is_complete = true;

        return;
    } else {
        sim->wing->last_position = sim->wing->position;
        sim->wing->last_rotation = sim->wing->rotation;
    }
}

void Simulation_ComputeKinematicVelocities(Simulation *sim) {
    size_t num_control_points = Simulation_GetNumPoints(sim, CONTROL_POINTS);

    if (!sim->iteration) {
        Vector3D zeros = {0.0, 0.0, 0.0};

        for (size_t i = 0; i < num_control_points; i++) {
            sim->kinematic_velocities[i] = zeros;
            sim->last_control_points[i] = sim->control_points[i];
        }

        return;
    }
    
    Vector3D current;
    Vector3D previous;
    Vector3D inverse_rotation = sim->wing->rotation;

    Vector3D *velocity;

    Vector3D_Multiply(&inverse_rotation, -1.0);
    
    double rotation_matrix[9];
    double last_rotation_matrix[9];
    double inverse_rotation_matrix[9];

    FillRotationMatrix(&sim->wing->rotation, rotation_matrix);
    FillRotationMatrix(&sim->wing->last_rotation, last_rotation_matrix);
    FillRotationMatrix(&inverse_rotation, inverse_rotation_matrix);

    for (size_t i = 0; i < num_control_points; i++) {
        current = sim->control_points[i];
        previous = sim->last_control_points[i];

        velocity = sim->kinematic_velocities + i;

        Vector3D_Rotate(&current, rotation_matrix);
        Vector3D_Rotate(&previous, last_rotation_matrix);
        Vector3D_Add(&current, &sim->wing->position, &current);
        Vector3D_Add(&previous, &sim->wing->last_position, &previous);
        Vector3D_Subtract(&current, &previous, velocity);
        Vector3D_Divide(velocity, sim->delta_time);

        sim->last_control_points[i] = sim->control_points[i];
    }
}

void Simulation_ShedWake(Simulation *sim) {
    Vector3D *point;
    Vector3D *position = &sim->wing->position;
    Vector3D *last_position = &sim->wing->last_position;

    double rotation_matrix[9];
    
    size_t num_points = Simulation_GetNumPoints(sim, WAKE_RING_POINTS);

    // Convert the wake points from local to global coordinates

    if (sim->iteration) {
        size_t last_num_points = num_points - Simulation_GetNumColumns(sim, WAKE_RING_POINTS);

        FillRotationMatrix(&sim->wing->last_rotation, rotation_matrix);

        for (size_t i = 0; i < last_num_points; i++) {
            point = sim->wake_ring_points + i;

            Vector3D_Rotate(point, rotation_matrix);
            Vector3D_Add(point, last_position, point);
        }
    }

    size_t icurr;
    size_t iprev;

    Vector3D point_copy;

    int num_rows = Simulation_GetNumRows(sim, BOUND_RING_POINTS);
    int num_cols = Simulation_GetNumColumns(sim, BOUND_RING_POINTS);

    FillRotationMatrix(&sim->wing->rotation, rotation_matrix);

    for (int j = 0; j < num_cols; j++) {
        // Shift the old wake points down one row

        if (sim->iteration) {
            for (int i = sim->iteration; i > 0; i--) {
                icurr = Sub2Ind(i, j, num_cols);
                iprev = Sub2Ind(i - 1, j, num_cols);

                sim->wake_ring_points[icurr] = sim->wake_ring_points[iprev];
            }
        }

        // Insert the new wake point into the first row, column j
        point_copy = sim->bound_ring_points[Sub2Ind(num_rows - 1, j, num_cols)];

        Vector3D_Rotate(&point_copy, rotation_matrix);
        Vector3D_Add(&point_copy, position, &point_copy);

        sim->wake_ring_points[Sub2Ind(0, j, num_cols)] = point_copy;
    }

    Vector3D opposite_rotation = sim->wing->rotation;

    Vector3D_Multiply(&opposite_rotation, -1.0);
    FillRotationMatrix(&opposite_rotation, rotation_matrix);

    // Convert wake points back from global to local coordinates
    for (size_t i = 0; i < num_points; i++) {
        point = sim->wake_ring_points + i;

        Vector3D_Subtract(point, position, point);
        Vector3D_Rotate(point, rotation_matrix);
    }

    // Assign trailing wing vorticity strengths to the first row of the wake

    num_rows = Simulation_GetNumRows(sim, CONTROL_POINTS);
    num_cols = Simulation_GetNumColumns(sim, CONTROL_POINTS);
    
    int num_wake_point_rows = Simulation_GetNumRows(sim, WAKE_RING_POINTS);

    size_t ibound;
    size_t iwake;

    if (num_wake_point_rows > 1) {
        for (int j = 0; j < num_cols; j++) {
            if (num_wake_point_rows > 2) {
                // Shift the old wake vortex strengths down one row

                for (int i = num_wake_point_rows - 2; i > 0; i--) {
                    icurr = Sub2Ind(i, j, num_cols);
                    iprev = Sub2Ind(i - 1, j, num_cols);

                    sim->wake_vortex_strengths[icurr] = sim->wake_vortex_strengths[iprev];
                }
            }

            iwake = Sub2Ind(0, j, num_cols);
            ibound = Sub2Ind(num_rows - 1, j, num_cols);

            // Insert the new vortex strength into the first row, column j
            sim->wake_vortex_strengths[ibound] = sim->bound_vortex_strengths[iwake];
        }
    }
}

void Simulation_InduceUnitVelocities(Simulation *sim, Geometry geometry, int i, int j, const Vector3D *point) {
    Vector3D **corners = sim->corner_buffer;
    Vector3D *velocities = sim->unit_velocity_buffer;

    Simulation_GetCorners(sim, geometry, i, j, corners);

    if (i > 0) {
        velocities[0].x = -sim->spanwise_velocity_buffer[j].x;
        velocities[0].y = -sim->spanwise_velocity_buffer[j].y;
        velocities[0].z = -sim->spanwise_velocity_buffer[j].z;
    } else {
        InduceUnitVelocity(point, corners[0], corners[1], velocities, sim->cutoff_radius);
    }

    InduceUnitVelocity(point, corners[1], corners[2], velocities + 1, sim->cutoff_radius);
    InduceUnitVelocity(point, corners[2], corners[3], velocities + 2, sim->cutoff_radius);

    if (j > 0) {
        velocities[3].x = -sim->chordwise_velocity_buffer.x;
        velocities[3].y = -sim->chordwise_velocity_buffer.y;
        velocities[3].z = -sim->chordwise_velocity_buffer.z;
    } else {
        InduceUnitVelocity(point, corners[3], corners[0], velocities + 3, sim->cutoff_radius);
    }

    sim->chordwise_velocity_buffer = velocities[1];
    sim->spanwise_velocity_buffer[j] = velocities[2];
}

void Simulation_ComputeCoefficients(Simulation *sim) {
    Vector3D w_induced;
    Vector3D v_induced;

    Vector3D point_copy;
    Vector3D normal_copy;

    Vector3D *unit_velocities = sim->unit_velocity_buffer;

    double a;
    double b;

    int num_rows = Simulation_GetNumRows(sim, BOUND_RING_POINTS) - 1;
    int num_cols = Simulation_GetNumColumns(sim, BOUND_RING_POINTS) - 1;

    size_t iring;
    size_t imatrix;
    size_t num_control_points = Simulation_GetNumPoints(sim, CONTROL_POINTS);
    size_t num_rings = Simulation_GetNumQuads(sim, BOUND_RING_POINTS);

    for (int mirror = 0; mirror < 2; mirror++) {
        for (size_t ipoint = 0; ipoint < num_control_points; ipoint++) {
            normal_copy = sim->normals[ipoint];
            point_copy = sim->control_points[ipoint];

            if (mirror) {
                point_copy.y = -point_copy.y;
                normal_copy.y = -normal_copy.y;
            }

            iring = 0;

            for (int i = 0; i < num_rows; i++) {
                for (int j = 0; j < num_cols; j++) {
                    Simulation_InduceUnitVelocities(sim, BOUND_RING_POINTS, i, j, &point_copy);
                                      
                    imatrix = ipoint * num_rings + iring;

                    Vector3D_Add(unit_velocities + 1, unit_velocities + 3, &w_induced);
                    Vector3D_Add(unit_velocities, unit_velocities + 2, &v_induced);
                    Vector3D_Add(&v_induced, &w_induced, &v_induced);

                    a = Vector3D_Dot(&v_induced, &normal_copy);
                    b = Vector3D_Dot(&w_induced, &normal_copy);

                    if (mirror || iring) {
                        sim->a_wing_on_wing[imatrix] += a;
                        sim->b_wing_on_wing[imatrix] += b;
                    } else {
                        sim->a_wing_on_wing[imatrix] = a;
                        sim->b_wing_on_wing[imatrix] = b;
                    }

                    iring++;
                }
            }
        }
    }
}

void Simulation_ComputeWakeInducedVelocities(Simulation *sim) {
    Vector3D induced_velocity;
    Vector3D point_copy;

    Vector3D *vels = sim->unit_velocity_buffer;

    int num_rows = Simulation_GetNumRows(sim, WAKE_RING_POINTS) - 1;
    int num_cols = Simulation_GetNumColumns(sim, WAKE_RING_POINTS) - 1;

    size_t iring;
    size_t imatrix;
    size_t num_control_points = Simulation_GetNumPoints(sim, CONTROL_POINTS);
    size_t num_rings = Simulation_GetNumQuads(sim, WAKE_RING_POINTS);

    for (int mirror = 0; mirror < 2; mirror++) {
        for (size_t ipoint = 0; ipoint < num_control_points; ipoint++) {
            point_copy = sim->control_points[ipoint];

            if (mirror) {
                point_copy.y = -point_copy.y;
            }

            iring = 0;

            for (int i = 0; i < num_rows; i++) {
                for (int j = 0; j < num_cols; j++) {
                    Simulation_InduceUnitVelocities(sim, WAKE_RING_POINTS, i, j, &point_copy);
                                      
                    imatrix = ipoint * num_rings + iring;

                    induced_velocity.x = vels[0].x + vels[1].x + vels[2].x + vels[3].x;
                    induced_velocity.y = vels[0].y + vels[1].y + vels[2].y + vels[3].y;
                    induced_velocity.z = vels[0].z + vels[1].z + vels[2].z + vels[3].z;

                    Vector3D_Multiply(&induced_velocity, sim->wake_vortex_strengths[iring]);

                    if (mirror || iring) {
                        Vector3D_Add(&induced_velocity, sim->wake_induced_velocities + imatrix, 
                                                        sim->wake_induced_velocities + imatrix);
                    } else {
                        sim->wake_induced_velocities[imatrix] = induced_velocity;
                    }

                    iring++;
                }
            }
        }
    }
}

void Simulation_Solve(Simulation *sim) {
    double rhs;

    Vector3D *normal;

    int info;
    int num_right_hand_sides = 1;
    int column_size = sim->wing->num_chordwise_panels * sim->wing->num_spanwise_panels;

    size_t imatrix;
    size_t num_rings = Simulation_GetNumQuads(sim, WAKE_RING_POINTS);
    size_t num_points = Simulation_GetNumPoints(sim, CONTROL_POINTS);

    for (size_t ipoint = 0; ipoint < num_points; ipoint++) {
        normal = sim->normals + ipoint;
        rhs = Vector3D_Dot(sim->kinematic_velocities + ipoint, normal);

        for (size_t iring = 0; iring < num_rings; iring++) {
            imatrix = ipoint * num_rings + iring;

            rhs -= Vector3D_Dot(sim->wake_induced_velocities + imatrix, normal);
        }

        sim->right_hand_side[ipoint] = rhs;
    }

    dgesv_(&column_size, &num_right_hand_sides, sim->a_wing_on_wing, &column_size,
           sim->pivot_vector, sim->right_hand_side, &column_size, &info);

    if (info < 0) {
        fprintf(stderr, "Simulation_Solve: argument %d provided to dgesv is invalid or illegal", -info);
    } else if (info > 0) {
        fprintf(stderr, "Simulation_Solve: dgesv encountered a singular matrix");
    } else {
        for (size_t ipoint = 0; ipoint < num_points; ipoint++) {
            sim->bound_vortex_strengths[ipoint] = sim->right_hand_side[ipoint];
        }
    }
}

void Simulation_RollupWake(Simulation *sim) {
    int num_rows;
    int num_cols;

    size_t iring;
    size_t num_points = Simulation_GetNumPoints(sim, WAKE_RING_POINTS);

    bool rollup_y;

    Vector3D point_copy;
    Vector3D induced_velocity;
    Vector3D *vels = sim->unit_velocity_buffer;

    Vector3D *points = sim->wake_ring_points;
    Vector3D *displacements = sim->wake_point_displacements;

    Geometry geometry;
    Geometry geometries[2] = {BOUND_RING_POINTS, WAKE_RING_POINTS};

    double gamma;
    double *gammas;

    for (int igeometry = 0; igeometry < 2; igeometry++) {
        geometry = geometries[igeometry];

        num_rows = Simulation_GetNumRows(sim, geometry) - 1;
        num_cols = Simulation_GetNumColumns(sim, geometry) - 1;

        if (geometry == BOUND_RING_POINTS) {
            gammas = sim->bound_vortex_strengths;
        } else {
            gammas = sim->wake_vortex_strengths;
        }

        for (int mirror = 0; mirror < 2; mirror++) {
            for (size_t ipoint = 0; ipoint < num_points; ipoint++) {
                rollup_y = ipoint % sim->wing->num_spanwise_panels;

                point_copy = points[ipoint];

                if (mirror) {
                    point_copy.y = -point_copy.y;
                }

                iring = 0;

                for (int i = 0; i < num_rows; i++) {
                    for (int j = 0; j < num_cols; j++) {
                        gamma = gammas[iring];

                        Simulation_InduceUnitVelocities(sim, geometry, i, j, &point_copy);

                        induced_velocity.x = (vels[0].x + vels[1].x + vels[2].x + vels[3].x) * gamma;
                        induced_velocity.z = (vels[0].z + vels[1].z + vels[2].z + vels[3].z) * gamma;

                        displacements[ipoint].z += induced_velocity.x * sim->delta_time;
                        displacements[ipoint].x += induced_velocity.z * sim->delta_time;

                        if (rollup_y) {
                            induced_velocity.y = (vels[0].y + vels[1].y + vels[2].y + vels[3].y) * gamma;

                            displacements[ipoint].y += induced_velocity.y * sim->delta_time;
                        }

                        iring++;
                    }
                }
            }
        }
    }

    for (size_t ipoint = 0; ipoint < num_points; ipoint++) {
        Vector3D_Add(points + ipoint, displacements + ipoint, points + ipoint);
        Vector3D_Multiply(displacements + ipoint, 0.0);
    }
}