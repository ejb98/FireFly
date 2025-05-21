#include <stdio.h>
#include <stdlib.h>

#include "wing.h"
#include "vector.h"
#include "isempty.h"
#include "construct_mesh.h"
#include "init_bound_rings.h"
#include "init_control_points.h"
#include "init_surface_panels.h"
#include "compute_coefficients.h"

int init_wing(Wing *wing) {
    int num_rows_wake = wing->num_wake_point_rows;
    int num_rows = wing->num_chordwise_panels;
    int num_cols = wing->num_spanwise_panels;
    int num_points = num_rows * num_cols;
    int mirror;
    int append;

    wing->surface_panels = construct_mesh(num_rows + 1, num_cols + 1);
    if (isempty(&wing->surface_panels)) {
        perror("init_wing: failed to allocate memory for surface_panels");
        return 1;
    }

    wing->control_points = construct_mesh(num_rows, num_cols);
    if (isempty(&wing->control_points)) {
        perror("init_wing: failed to allocate memory for control_points");
        return 1;
    }

    wing->normal_vectors = construct_mesh(num_rows, num_cols);
    if (isempty(&wing->normal_vectors)) {
        perror("init_wing: failed to allocate memory for normal_vectors");
        return 1;
    }

    
    wing->bound_rings = construct_mesh(num_rows + 1, num_cols + 1);
    if (isempty(&wing->bound_rings)) {
        perror("init_wing: failed to allocate memory for bound_rings");
        return 1;
    }

    wing->wake_rings = construct_mesh(num_rows_wake, num_cols + 1);
    if (isempty(&wing->wake_rings)) {
        perror("init_wing: failed to allocate memory for wake_rings");
        return 1;
    }
    
    wing->vorticity_strengths = (double *) calloc(num_points, sizeof(double));
    if (!wing->vorticity_strengths) {
        perror("init_wing: failed to allocate memory for vorticity_strengths");
        return 1;
    }

    wing->normal_velocities = (double *) calloc(num_points, sizeof(double));
    if (!wing->normal_velocities) {
        perror("init_wing: failed to allocate memory for normal_velocities");
        return 1;
    }

    wing->a_wing_on_wing = (double *) malloc(sizeof(double) * num_points * num_points);
    if (!wing->a_wing_on_wing) {
        perror("init_wing: failed to allocate memory for a_wing_on_wing");
        return 1;
    }

    wing->b_wing_on_wing = (double *) malloc(sizeof(double) * num_points * num_points);
    if (!wing->b_wing_on_wing) {
        perror("init_wing: failed to allocate memory for b_wing_on_wing");
        return 1;
    }

    wing->pivot_vector = (int *) malloc(sizeof(int) * num_points);
    if (!wing->pivot_vector) {
        perror("init_wing: failed to allocate memory for pivot_vector");
        return 1;
    }

    wing->buffer = (Vector *) malloc(sizeof(Vector) * num_cols * 4);
    if (!wing->buffer) {
        perror("init_wing: failed to allocate memory for buffer");
        return 1;
    }

    init_surface_panels(wing);
    init_control_points(wing);
    init_bound_rings(wing);

    for (int i = 0; i < 2; i++) {
        mirror = i;
        append = i;

        compute_coefficients(&wing->control_points,
                             &wing->bound_rings,
                             &wing->normal_vectors,
                             wing->a_wing_on_wing,
                             wing->b_wing_on_wing,
                             wing->buffer,
                             wing->cutoff, mirror, append);
    }

    return 0;
}