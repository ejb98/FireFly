#include <stdlib.h>

#include "wing.h"
#include "dgesv.h"
#include "sub2ind.h"
#include "get_size.h"
#include "geometry.h"
#include "compute_coefficients.h"

void solve(Wing *wing) {
    if (wing->iteration) {
        if (wing->iteration > 1) {
            size_t num_points = (size_t) get_size(&wing->control_points);
            size_t num_rows = (size_t) wing->wake_rings.num_rows - 1;
            size_t num_cols = (size_t) wing->wake_rings.num_cols - 1;
            size_t num_rings = num_rows * num_cols;
            size_t num_coef = num_points * num_rings;

            for (size_t i = 0; i < num_coef; i++) {
                wing->a_wake_on_wing[i] = 0.0;
                wing->b_wake_on_wing[i] = 0.0;
            }

            for (size_t i = 0; i < get_size(&wing->wake_displacements); i++) {
                wing->wake_displacements.x[i] = 0.0;
                wing->wake_displacements.y[i] = 0.0;
                wing->wake_displacements.z[i] = 0.0;
            }
        }

        compute_coefficients(wing, WAKE_RINGS);
    } else {
        compute_coefficients(wing, BOUND_RINGS);
    }

    if (wing->iteration) {
        int info;
        int num_right_hand_sides = 1;
        int column_size = wing->num_chordwise_panels * wing->num_spanwise_panels;

        size_t imatrix;
        size_t num_rings = ((size_t) (wing->wake_rings.num_rows - 1)) * (wing->wake_rings.num_cols - 1);
        size_t num_points = get_size(&wing->control_points);

        for (size_t i = 0; i < num_points; i++) {
            wing->right_hand_side[i] = wing->normal_velocities[i];

            for (size_t j = 0; j < num_rings; j++) {
                imatrix = sub2ind(i, j, num_rings);

                wing->right_hand_side[i] -= wing->a_wake_on_wing[imatrix] * wing->wake_vorticity[j];
            }
        }

        dgesv_(&column_size, &num_right_hand_sides, wing->a_wing_on_wing,
               &column_size, wing->pivot_vector, wing->right_hand_side, &column_size, &info);

        for (size_t i = 0; i < num_points; i++) {
            wing->bound_vorticity[i] = wing->right_hand_side[i];
        }
    }
}