#include "wing.h"
#include "dgesv.h"
#include "sub2ind.h"
#include "compute_coefficients.h"

void solve(Wing *wing) {
    int mirror;
    int append;

    for (int i = 0; i < 2; i++) {
        mirror = i;
        append = i;

        if (wing->iteration) {
            compute_coefficients(&wing->control_points,
                                 &wing->wake_rings,
                                 &wing->normal_vectors,
                                 wing->a_wake_on_wing,
                                 wing->b_wake_on_wing,
                                 wing->buffer,
                                 wing->cutoff, mirror, append);
        } else {
            compute_coefficients(&wing->control_points,
                                 &wing->bound_rings,
                                 &wing->normal_vectors,
                                 wing->a_wing_on_wing,
                                 wing->b_wing_on_wing,
                                 wing->buffer,
                                 wing->cutoff, mirror, append);
        }
    }

    if (wing->iteration) {
        int info;
        int imatrix;
        int num_points = wing->num_chordwise_panels * wing->num_spanwise_panels;
        int num_right_hand_sides = 1;

        for (int i = 0; i < num_points; i++) {
            wing->right_hand_side[i] = wing->normal_velocities[i];

            for (int j = 0; j < wing->iteration; j++) {
                imatrix = sub2ind(i, j, wing->iteration);

                wing->right_hand_side[i] -= wing->a_wake_on_wing[imatrix] * wing->wake_vorticity[j];
            }
        }

        dgesv_(&num_points, &num_right_hand_sides, wing->a_wing_on_wing,
               &num_points, wing->pivot_vector, wing->right_hand_side, &num_points, &info);

        for (int i = 0; i < num_points; i++) {
            wing->bound_vorticity[i] = wing->right_hand_side[i];
        }
    }
}