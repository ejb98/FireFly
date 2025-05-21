#include <stdio.h>
#include <stdlib.h>

#include "wing.h"
#include "vector.h"
#include "isempty.h"
#include "check_null.h"
#include "construct_mesh.h"
#include "init_bound_rings.h"
#include "init_control_points.h"
#include "init_surface_panels.h"

int init_wing(Wing *wing) {
    wing->iteration = -1;
    
    int num_rows_wake = wing->num_wake_point_rows;
    int num_rows = wing->num_chordwise_panels;
    int num_cols = wing->num_spanwise_panels;
    int num_points = num_rows * num_cols;
    int num_wake_vortices = (num_rows_wake - 1) * (num_cols);

    wing->surface_panels = construct_mesh(num_rows + 1, num_cols + 1);
    if (isempty(&wing->surface_panels)) {
        fprintf(stderr, "init_wing: failed to allocate memory for field surface_panels");
        return 1;
    }

    wing->control_points = construct_mesh(num_rows, num_cols);
    if (isempty(&wing->control_points)) return 1;

    wing->normal_vectors = construct_mesh(num_rows, num_cols);
    if (isempty(&wing->normal_vectors)) return 1;

    wing->bound_rings = construct_mesh(num_rows + 1, num_cols + 1);
    if (isempty(&wing->bound_rings)) return 1;

    wing->wake_rings = construct_mesh(num_rows_wake, num_cols + 1);
    if (isempty(&wing->wake_rings)) return 1;
    
    wing->wake_vorticity = (double *) calloc(num_wake_vortices, sizeof(double));
    if (check_null("init_wing", "wake_vorticity", wing->wake_vorticity)) return 1;

    wing->bound_vorticity = (double *) calloc(num_points, sizeof(double));
    if (check_null("init_wing", "bound_vorticity", wing->bound_vorticity)) return 1;

    wing->normal_velocities = (double *) calloc(num_points, sizeof(double));
    if (check_null("init_wing", "normal_velocities", wing->normal_velocities)) return 1;

    wing->right_hand_side = (double *) malloc(sizeof(double) * num_points);
    if (check_null("init_wing", "right_hand_side", wing->right_hand_side)) return 1;

    wing->a_wing_on_wing = (double *) malloc(sizeof(double) * num_points * num_points);
    if (check_null("init_wing", "a_wing_on_wing", wing->a_wing_on_wing)) return 1;

    wing->b_wing_on_wing = (double *) malloc(sizeof(double) * num_points * num_points);
    if (check_null("init_wing", "b_wing_on_wing", wing->b_wing_on_wing)) return 1;

    wing->a_wake_on_wing = (double *) malloc(sizeof(double) * num_points * num_wake_vortices);
    if (check_null("init_wing", "a_wake_on_wing", wing->a_wake_on_wing)) return 1;

    wing->b_wake_on_wing = (double *) malloc(sizeof(double) * num_points * num_wake_vortices);
    if (check_null("init_wing", "b_wake_on_wing", wing->b_wake_on_wing)) return 1;

    wing->pivot_vector = (int *) malloc(sizeof(int) * num_points);
    if (check_null("init_wing", "pivot_vector", wing->pivot_vector)) return 1;

    wing->buffer = (Vector *) malloc(sizeof(Vector) * num_cols * 4);
    if (check_null("init_wing", "buffer", wing->buffer)) return 1;

    init_surface_panels(wing);
    init_control_points(wing);
    init_bound_rings(wing);

    return 0;
}