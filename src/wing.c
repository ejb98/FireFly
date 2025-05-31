#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#include "wing.h"
#include "sub2ind.h"
#include "vector3d.h"
#include "constants.h"
#include "allocate_doubles.h"

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
                     double trailing_edge_sweep_angle) {

    Wing *wing = (Wing *) malloc(sizeof(Wing));

    if (wing == NULL) {
        fprintf(stderr, "Wing_construct: memory allocation for struct Wing returned NULL");
    }

    Vector3D zeros = {0.0, 0.0, 0.0};

    size_t num_panels = num_spanwise_panels * num_chordwise_panels;
    size_t num_points = (num_spanwise_panels + 1) * (num_chordwise_panels + 1);
    size_t max_num_wake_rings = (num_time_steps - 1) * num_spanwise_panels;
    size_t max_num_wake_points = num_time_steps * (num_spanwise_panels + 1);
    
    wing->naca_m = naca_m;
    wing->naca_p = naca_p;
    wing->iteration = 0;
    wing->num_time_steps = num_time_steps;
    wing->num_spanwise_panels = num_spanwise_panels;
    wing->num_chordwise_panels = num_chordwise_panels;

    wing->pivot_vector = (int *) malloc(sizeof(int) * num_panels);

    if (wing->pivot_vector == NULL) {
        fprintf(stderr, "Wing_construct: memory allocation for pivot_vector returned NULL");
    }

    wing->semi_span = semi_span;
    wing->root_chord = root_chord;
    wing->cutoff_radius = cutoff_radius;
    wing->angle_of_attack = angle_of_attack;
    wing->starting_vortex_offset = starting_vortex_offset;
    wing->leading_edge_sweep_angle = leading_edge_sweep_angle;
    wing->trailing_edge_sweep_angle = trailing_edge_sweep_angle;

    wing->pressures = allocate_doubles(num_panels);
    wing->a_wing_on_wing = allocate_doubles(num_panels * num_panels);
    wing->b_wing_on_wing = allocate_doubles(num_panels * num_panels);
    wing->right_hand_side = allocate_doubles(num_panels);
    wing->wake_vortex_strengths = allocate_doubles(max_num_wake_rings);
    wing->bound_vortex_strengths = allocate_doubles(num_panels);
    wing->previous_bound_vortex_strengths = allocate_doubles(num_panels);

    wing->position = zeros;
    wing->rotation = zeros;
    wing->previous_position = zeros;
    wing->previous_rotation = zeros;
    wing->vertical_velocity_buffer = zeros;

    wing->normal_vectors = Vector3D_malloc(num_panels);
    wing->surface_points = Vector3D_malloc(num_points);
    wing->control_points = Vector3D_malloc(num_panels);
    wing->wake_ring_points = Vector3D_malloc(max_num_wake_points);
    wing->bound_ring_points = Vector3D_malloc(num_points);
    wing->kinematic_velocities = Vector3D_malloc(num_panels);
    wing->wake_induced_velocities = Vector3D_malloc(num_panels * max_num_wake_rings);
    wing->wake_point_displacements = Vector3D_malloc(max_num_wake_points);
    wing->spanwise_tangent_vectors = Vector3D_malloc(num_panels);
    wing->chordwise_tangent_vectors = Vector3D_malloc(num_panels);
    wing->horizontal_velocity_buffer = Vector3D_malloc(wing->num_spanwise_panels);

    Wing_compute_surface_points(wing);

    return wing;
}

void Wing_destruct(Wing *wing) {    
    free(wing->pivot_vector);

    free(wing->pressures);
    free(wing->a_wing_on_wing);
    free(wing->b_wing_on_wing);
    free(wing->right_hand_side);
    free(wing->wake_vortex_strengths);
    free(wing->bound_vortex_strengths);
    free(wing->previous_bound_vortex_strengths);

    free(wing->normal_vectors);
    free(wing->surface_points);
    free(wing->control_points);
    free(wing->wake_ring_points);
    free(wing->bound_ring_points);
    free(wing->kinematic_velocities);
    free(wing->wake_induced_velocities);
    free(wing->wake_point_displacements);
    free(wing->spanwise_tangent_vectors);
    free(wing->chordwise_tangent_vectors);
    free(wing->horizontal_velocity_buffer);

    free(wing);
}

void Wing_compute_surface_points(Wing *wing) {
    double constant;
    double local_chord;
    double normalized_x;
    double normalized_z;
    double leading_offset;
    double trailing_offset;
    double rotation_matrix[3][3];
    double p = wing->naca_p / 10.0;
    double m = wing->naca_m / 100.0;
    double leading_tangent = tan((90.0 - wing->leading_edge_sweep_angle) * PI / 180.0);
    double trailing_tangent = tan((90.0 - wing->trailing_edge_sweep_angle) * PI / 180.0);

    Vector3D rotation = {0.0, wing->angle_of_attack * PI / 180.0, 0.0};
    Vector3D *point;

    Vector3D_fill_rotation_matrix(&rotation, rotation_matrix);

    for (int j = 0; j < wing->num_spanwise_panels + 1; j++) {
        for (int i = 0; i < wing->num_chordwise_panels + 1; i++) {
            point = wing->surface_points + sub2ind(i, j, wing->num_spanwise_panels + 1);

            point->y = wing->semi_span * j / wing->num_spanwise_panels;
            normalized_x = ((double) i) / wing->num_chordwise_panels;

            constant = 2.0 * p * normalized_x - normalized_x * normalized_x;

            if (normalized_x >= p || !wing->naca_p) {
                normalized_z = (m / ((1.0 - p) * (1.0 - p))) * (1.0 - 2.0 * p + constant);
            } else {
                normalized_z = (m / (p * p)) * constant;
            }

            leading_offset = point->y * leading_tangent;
            trailing_offset = point->y * trailing_tangent;

            local_chord = wing->root_chord + trailing_offset - leading_offset;

            point->x = leading_offset + normalized_x * local_chord;
            point->z = local_chord * normalized_z;

            Vector3D_rotate(point, rotation_matrix);
        }
    }
}

void Wing_write_points_to_vtk(Wing *wing, Geometry geometry, const char *file_path) {
    size_t length = strlen(file_path);

    char back_slash = '\\';
    char forward_slash = '/';
    char last_character = file_path[length - 1];

    if (last_character != forward_slash && last_character != back_slash) {
        fprintf(stderr, "Wing_write_points_to_vtk: last character of file path must be either '\\' or '/'");

        return;
    }

    int num_rows;
    int num_cols;

    char full_path[175];
    char description[20];

    Vector3D *points;

    switch (geometry) {
        case CONTROL_POINTS:
            strcpy(description, "control_points");
            points = wing->control_points;
            num_rows = wing->num_chordwise_panels;
            num_cols = wing->num_spanwise_panels;

            break;
        case SURFACE_POINTS:
            strcpy(description, "surface_points");
            points = wing->surface_points;
            num_rows = wing->num_chordwise_panels + 1;
            num_cols = wing->num_spanwise_panels + 1;

            break;
        case WAKE_RING_POINTS:
            strcpy(description, "wake_ring_points");
            points = wing->wake_ring_points;
            num_rows = wing->iteration;
            num_cols = wing->num_spanwise_panels + 1;

            break;
        case BOUND_RING_POINTS:
            strcpy(description, "bound_ring_points");
            points = wing->bound_ring_points;
            num_rows = wing->num_chordwise_panels + 1;
            num_cols = wing->num_spanwise_panels + 1;

            break;
        default:
            fprintf(stderr, "Wing_write_points_to_vtk: inapplicable geometry cannot be written to file");

            return;
    }
    
    if (num_rows < 2) {
        fprintf(stderr, "Wing_write_points_to_vtk: points must have at least two rows");

        return;
    }

    if (num_cols < 2) {
        fprintf(stderr, "Wing_write_points_to_vtk: points must have at least two columns");

        return;
    }

    snprintf(full_path, sizeof(full_path), "%s%s.vtk.%d", file_path, description, wing->iteration);

    FILE *file = fopen(full_path, "w");

    if (file == NULL) {
        fprintf(stderr, "Wing_write_points_to_vtk: failed to open %s", full_path);

        return;
    }

    size_t i0, i1, i2, i3;
    size_t num_quads = (num_rows - 1) * (num_cols - 1);
    size_t num_points = num_rows * num_cols;

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
            i0 = sub2ind(i, j, num_cols);
            i1 = sub2ind(i, j + 1, num_cols);
            i2 = sub2ind(i + 1, j + 1, num_cols);
            i3 = sub2ind(i + 1, j, num_cols);

            fprintf(file, "4 %zu %zu %zu %zu\n", i0, i1, i2, i3);
        }
    }

    fprintf(file, "CELL_TYPES %zu\n", num_quads);
    for (size_t i = 0; i < num_quads; i++) {
        fprintf(file, "9\n");
    }

    fclose(file);
}