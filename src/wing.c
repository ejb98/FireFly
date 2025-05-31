#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#include "wing.h"
#include "sub2ind.h"
#include "vector3d.h"
#include "constants.h"
#include "allocate_doubles.h"
#include "fill_rotation_matrix.h"

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
    wing->surface_area = 0.0;
    wing->cutoff_radius = cutoff_radius;
    wing->angle_of_attack = angle_of_attack;
    wing->starting_vortex_offset = starting_vortex_offset;
    wing->leading_edge_sweep_angle = leading_edge_sweep_angle;
    wing->trailing_edge_sweep_angle = trailing_edge_sweep_angle;

    wing->pressures = allocate_doubles(num_panels);
    wing->a_wing_on_wing = allocate_doubles(num_panels * num_panels);
    wing->b_wing_on_wing = allocate_doubles(num_panels * num_panels);
    wing->surface_areas = allocate_doubles(num_panels);
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
    Wing_compute_surface_vectors(wing);
    Wing_compute_surface_areas(wing);
    Wing_compute_control_points(wing);
    Wing_compute_bound_ring_points(wing);

    return wing;
}

void Wing_free(Wing *wing) {    
    free(wing->pivot_vector);

    free(wing->pressures);
    free(wing->a_wing_on_wing);
    free(wing->b_wing_on_wing);
    free(wing->surface_areas);
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

void Wing_print_attributes(const Wing *wing) {
    printf("Semi Span: %.2f m\n", wing->semi_span);
    printf("Root Chord: %.2f m\n", wing->root_chord);
    printf("Leading Edge Sweep Angle: %.2f deg\n", wing->leading_edge_sweep_angle);
    printf("Trailing Edge Sweep Angle: %.2f deg\n", wing->trailing_edge_sweep_angle);
    printf("Surface Area: %.2f sq. m\n", wing->surface_area);
    printf("Angle of Attack: %.2f deg\n", wing->angle_of_attack);
    printf("Airfoil: NACA %d%dXX\n", wing->naca_m, wing->naca_p);
    printf("Spanwise Panels: %d\n", wing->num_spanwise_panels);
    printf("Chordwise Panels: %d\n", wing->num_chordwise_panels);
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

    fill_rotation_matrix(&rotation, rotation_matrix);

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

Vector3D *Wing_get_points(const Wing *wing, Geometry geometry, int *num_rows, int *num_cols) {
    Vector3D *points;

    if (geometry == CONTROL_POINTS) {
        *num_cols = wing->num_spanwise_panels;
    } else {
        *num_cols = wing->num_spanwise_panels + 1;
    }

    switch (geometry) {
        case CONTROL_POINTS:
            points = wing->control_points;
            *num_rows = wing->num_chordwise_panels;

            break;
        case SURFACE_POINTS:
            points = wing->surface_points;
            *num_rows = wing->num_chordwise_panels + 1;

            break;
        case WAKE_RING_POINTS:
            points = wing->wake_ring_points;
            *num_rows = wing->iteration;

            break;
        case BOUND_RING_POINTS:
            points = wing->bound_ring_points;
            *num_rows = wing->num_chordwise_panels + 1;
    }

    return points;
}

void Wing_get_corners(const Wing *wing, Geometry geometry, int i, int j, Vector3D **corners) {
    int num_rows;
    int num_cols;

    Vector3D *points = Wing_get_points(wing, geometry, &num_rows, &num_cols);

    size_t indices[4];
    size_t num_points = ((size_t) num_rows) * num_cols;

    indices[0] = sub2ind(i, j, num_cols);
    indices[1] = sub2ind(i, j + 1, num_cols);
    indices[2] = sub2ind(i + 1, j + 1, num_cols);
    indices[3] = sub2ind(i + 1, j, num_cols);

    if (indices[2] > num_points - 1) {
        fprintf(stderr, "Wing_get_corners: index %zu is out of bounds for array of size %zu", indices[2], num_points);
    }

    for (int i = 0; i < 4; i++) {
        corners[i] = points + indices[i];
    }
}

void Wing_write_points_to_vtk(const Wing *wing, Geometry geometry, const char *file_path) {
    int num_rows;
    int num_cols;

    char full_path[175];
    char description[20];

    Vector3D *points = Wing_get_points(wing, geometry, &num_rows, &num_cols);

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

void Wing_compute_surface_vectors(Wing *wing) {
    size_t ipanel;
    Vector3D left;
    Vector3D back;
    Vector3D front;
    Vector3D right;
    Vector3D vectora;
    Vector3D vectorb;
    Vector3D *normal;
    Vector3D *corners[4];

    for (int i = 0; i < wing->num_chordwise_panels; i++) {
        for (int j = 0; j < wing->num_spanwise_panels; j++) {
            ipanel = sub2ind(i, j, wing->num_spanwise_panels);
            normal = wing->normal_vectors + ipanel;
            Wing_get_corners(wing, SURFACE_POINTS, i, j, corners);
            Vector3D_lerp(corners[0], corners[3], 0.5, &left);
            Vector3D_lerp(corners[2], corners[3], 0.5, &back);
            Vector3D_lerp(corners[0], corners[1], 0.5, &front);
            Vector3D_lerp(corners[1], corners[2], 0.5, &right);
            Vector3D_subtract(corners[2], corners[0], &vectora);
            Vector3D_subtract(corners[1], corners[3], &vectorb);
            Vector3D_cross(&vectora, &vectorb, normal);
            Vector3D_normalize(normal);
            Vector3D_direction(&front, &back, wing->chordwise_tangent_vectors + ipanel);
            Vector3D_direction(&left, &right, wing->spanwise_tangent_vectors + ipanel);
        }
    }
}

void Wing_compute_control_points(Wing *wing) {
    size_t ipanel;
    Vector3D left;
    Vector3D right;
    Vector3D *corners[4];

    for (int i = 0; i < wing->num_chordwise_panels; i++) {
        for (int j = 0; j < wing->num_spanwise_panels; j++) {
            ipanel = sub2ind(i, j, wing->num_spanwise_panels);
            Wing_get_corners(wing, SURFACE_POINTS, i, j, corners);
            Vector3D_lerp(corners[0], corners[3], 0.75, &left);
            Vector3D_lerp(corners[1], corners[2], 0.75, &right);
            Vector3D_lerp(&left, &right, 0.5, wing->control_points + ipanel);
        }
    }
}

void Wing_compute_surface_areas(Wing *wing) {
    size_t ipanel;
    Vector3D *corners[4];

    double a1, a2, b1, b2;

    wing->surface_area = 0.0;

    for (int i = 0; i < wing->num_chordwise_panels; i++) {
        for (int j = 0; j < wing->num_spanwise_panels; j++) {
            ipanel = sub2ind(i, j, wing->num_spanwise_panels);
            Wing_get_corners(wing, SURFACE_POINTS, i, j, corners);

            a1 = Vector3D_distance(corners[0], corners[1]);
            b1 = Vector3D_distance(corners[1], corners[2]);
            a2 = Vector3D_distance(corners[2], corners[3]);
            b2 = Vector3D_distance(corners[3], corners[0]);

            wing->surface_areas[ipanel] = 0.5 * (a1 * b1 + a2 * b2);
            wing->surface_area += wing->surface_areas[ipanel];
        }
    }
}

void Wing_compute_bound_ring_points(Wing *wing) {
    size_t ipoint;

    Vector3D *next;
    Vector3D *point;

    int num_rows = wing->num_chordwise_panels + 1;
    int num_cols = wing->num_spanwise_panels + 1;

    for (int j = 0; j < num_cols; j++) {
        for (int i = 0; i < num_rows; i++) {
            ipoint = sub2ind(i, j, num_cols);
            point = wing->surface_points + ipoint;

            if (i == num_rows - 1) {
                wing->bound_ring_points[ipoint].x = point->x + wing->starting_vortex_offset;
                wing->bound_ring_points[ipoint].y = point->y;
                wing->bound_ring_points[ipoint].z = point->z;
            } else {
                next = wing->surface_points + sub2ind(i + 1, j, num_cols);
                Vector3D_lerp(point, next, 0.25, wing->bound_ring_points + ipoint);
            }
        }
    }
}