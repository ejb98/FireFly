#include <stdio.h>
#include <stdlib.h>

#include "wing.h"
#include "arena.h"
#include "vector.h"
#include "isempty.h"
#include "check_null.h"
#include "construct_mesh.h"
#include "arena_allocate.h"
#include "init_bound_rings.h"
#include "init_control_points.h"
#include "init_surface_panels.h"

int init_wing(Wing *wing) {
    wing->iteration = -1;
    
    int num_rows = wing->num_chordwise_panels;
    int num_cols = wing->num_spanwise_panels;
    int num_step = wing->num_wake_point_rows;

    size_t nr = num_rows;
    size_t nrp1 = nr + 1;
    size_t nc = num_cols;
    size_t ncp1 = nc + 1;
    size_t nt = num_step;
    size_t ntm1 = nt - 1;
    size_t nrnc = nr * nc;
    size_t np = 6 * (nrnc + ncp1 * (nt + nrp1)) + ntm1 * nc + nrnc * (3 + 2 * nc * (nr + ntm1));

    wing->memory.num_elements = np;
    wing->memory.next_free_index = 0;
    wing->memory.elements = (double *) calloc(sizeof(double), np);

    Arena *arena = &wing->memory;

    if (wing->memory.elements == NULL) {
        fprintf(stderr, "init_wing: failed to allocate memory arena for %zu points", np);

        return 1;
    }

    wing->surface_panels = construct_mesh(num_rows + 1, num_cols + 1, arena);
    wing->control_points = construct_mesh(num_rows, num_cols, arena);
    wing->normal_vectors = construct_mesh(num_rows, num_cols, arena);
    wing->bound_rings = construct_mesh(num_rows + 1, num_cols + 1, arena);
    wing->wake_rings = construct_mesh(num_step, num_cols + 1, arena);
    wing->wake_displacements = construct_mesh(num_step, num_cols + 1, arena);    
    wing->wake_vorticity = arena_allocate((nt - 1) * nc, arena);
    wing->bound_vorticity = arena_allocate(nr * nc, arena);
    wing->normal_velocities = arena_allocate(nr * nc, arena);
    wing->right_hand_side = arena_allocate(nr * nc, arena);
    wing->a_wing_on_wing = arena_allocate(nr * nc * nr * nc, arena);
    wing->b_wing_on_wing = arena_allocate(nr * nc * nr * nc, arena);
    wing->a_wake_on_wing = arena_allocate(nr * nc * (nt - 1) * nc, arena);
    wing->b_wake_on_wing = arena_allocate(nr * nc * (nt - 1) * nc, arena);
    wing->pivot_vector = (int *) calloc(nr * nc, sizeof(int));
    
    if (wing->pivot_vector == NULL) {
        fprintf(stderr, "init_wing: failed to allocate pivot vector of length %zu", nr * nc);

        return 1;
    }

    init_surface_panels(wing);
    init_control_points(wing);
    init_bound_rings(wing);

    return 0;
}