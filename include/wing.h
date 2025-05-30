#ifndef WING_H
#define WING_H

#include <stddef.h>

#include "arena.h"
#include "vector3d.h"
#include "constants.h"

typedef struct Wing {
    int naca_m;
    int naca_p;
    int iteration;
    int num_time_steps;
    int num_wake_rings;
    int num_wake_points;
    int num_control_points;
    int num_surface_points;
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

    DoubleArena double_arena;
    Vector3DArena vector3d_arena;
} Wing;

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

    Vector3D rotation = {0.0, wing->angle_of_attack * pi / 180.0, 0.0};
    Vector3D *point;

    Vector3D_fill_rotation_matrix(&rotation, rotation_matrix);

    for (int j = 0; j < wing->num_chordwise_panels + 1; j++) {
        for (int i = 0; i < wing->num_spanwise_panels + 1; i++) {
            point = wing->surface_points + sub2ind(i, j, wing->num_spanwise_panels + 1);

            point->y = wing->semi_span * j / wing->num_spanwise_panels;
            normalized_x = ((double) i) / wing->num_chordwise_panels;

            constant = 2.0 * p * normalized_x - normalized_x * normalized_x;

            if (normalized_x >= p || !wing->naca_p) {
                normalized_z = (m / ((1.0 - p) * (1.0 - p))) * (1.0 - 2.0 * p + constant);
            } else {
                normalized_z = (m / (p * p)) * constant;
            }

            leading_offset = point->y * tan_leading;
            trailing_offset = point->y * tan_trailing;

            local_chord = wing->root_chord + trailing_offset - leading_offset;

            point->x = leading_offset + normalized_x * local_chord;
            point->z = local_chord * normalized_z;

            Vector3D_rotate(point, rotation_matrix);
        }
    }
}

void Wing_init_memory_arenas(Wing *wing) {
    size_t num_control_points = wing->num_control_points;
    size_t num_normal_vectors = num_control_points;
    size_t num_surface_points = wing->num_surface_points;
    size_t num_wake_ring_points = wing->num_time_steps * (wing->num_spanwise_panels + 1);
    size_t num_bound_ring_points = num_surface_points;
    size_t num_kinematic_velocities = num_control_points;
    size_t num_wake_induced_velocities = num_control_points * (wing->num_time_steps - 1) * (wing->num_spanwise_panels);
    size_t num_wake_point_displacements = num_wake_ring_points;
    size_t num_spanwise_tangent_vectors = num_control_points;
    size_t num_chordwise_tangent_vectors = num_control_points;
    size_t num_horizontal_velocity_buffer_elements = wing->num_spanwise_panels;

    size_t num_pressures = num_control_points;
    size_t num_a_wing_on_wing_elements = num_control_points * num_control_points;
    size_t num_b_wing_on_wing_elements = num_a_wing_on_wing_elements;
    size_t num_right_hand_side_elements = num_control_points;
    size_t num_wake_vortex_strengths = (wing->num_time_steps - 1) * (wing->num_spanwise_panels);
    size_t num_bound_vortex_strengths = num_control_points;
    size_t num_previous_bound_vortex_strengths = num_control_points;

    size_t num_vector3d_elements = num_control_points +
                                   num_normal_vectors +
                                   num_surface_points +
                                   num_wake_ring_points +
                                   num_bound_ring_points +
                                   num_kinematic_velocities +
                                   num_wake_induced_velocities +
                                   num_wake_point_displacements +
                                   num_spanwise_tangent_vectors +
                                   num_chordwise_tangent_vectors +
                                   num_horizontal_velocity_buffer_elements;

    size_t num_double_elements = num_pressures +
                                 num_a_wing_on_wing_elements +
                                 num_b_wing_on_wing_elements +
                                 num_right_hand_side_elements +
                                 num_wake_vortex_strengths +
                                 num_bound_vortex_strengths +
                                 num_previous_bound_vortex_strengths;
    
    wing->vector3d_arena = Vector3DArena_construct(num_vector3d_elements);
    wing->double_arena = DoubleArena_allocate(num_double_elements);
}

Wing Wing_construct(int naca_m,
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

    Wing wing_obj = {.naca_m = naca_m,
                     .naca_p = naca_p,
                     .iteration = 0,
                     .num_wake_rings = 0,
                     .num_wake_points = 0,
                     .num_time_steps = num_time_steps,
                     .num_spanwise_panels = num_spanwise_panels,
                     .num_chordwise_panels = num_chordwise_panels,
                     .semi_span = semi_span,
                     .root_chord = root_chord,
                     .cutoff_radius = cutoff_radius,
                     .angle_of_attack = angle_of_attack,
                     .starting_vortex_offset = starting_vortex_offset,
                     .leading_edge_sweep_angle = leading_edge_sweep_angle,
                     .trailing_edge_sweep_angle = trailing_edge_sweep_angle,
                     .num_control_points = num_chordwise_panels * num_spanwise_panels,
                     .num_surface_points = (num_chordwise_panels + 1) * (num_spanwise_panels + 1)};

    Wing *wing = &wing;
    
    Wing_init_memory_arenas(wing);
    Wing_compute_surface_points(wing);
}

#endif