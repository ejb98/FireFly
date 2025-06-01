#include <stdlib.h>

#include "add.h"
#include "wing.h"
#include "sub2ind.h"
#include "get_size.h"
#include "subtract.h"
#include "vector_to_mesh.h"
#include "mesh_to_vector.h"
#include "apply_rotation.h"
#include "assign_rotation.h"

void shed_wake(Wing *wing) {
    Vector3D point;

    double rot_mat[3][3];

    size_t iwing;
    size_t iring;
    size_t iwake;
    size_t icurr;
    size_t iprev;
    size_t num_points;

    if (wing->iteration) {
        num_points = get_size(&wing->wake_rings);

        assign_rotation(rot_mat, &wing->rotation_prev);

        for (size_t i = 0; i < num_points; i++) {
            mesh_to_vector(&wing->wake_rings, i, &point);
            apply_rotation(rot_mat, &point);
            add(&point, &wing->position_prev, &point);
            vector_to_mesh(&point, &wing->wake_rings, i);
        }
    }

    assign_rotation(rot_mat, &wing->rotation);

    for (int j = 0; j < wing->bound_rings.num_cols; j++) {
        iring = Sub2Ind(wing->bound_rings.num_rows - 1, j, wing->bound_rings.num_cols);
        iwake = Sub2Ind(0, j, wing->wake_rings.num_cols);

        mesh_to_vector(&wing->bound_rings, iring, &point);
        apply_rotation(rot_mat, &point);

        if (wing->iteration) {
            for (int i = wing->iteration; i > 0; i--) {
                icurr = Sub2Ind(i, j, wing->wake_rings.num_cols);
                iprev = Sub2Ind(i - 1, j, wing->wake_rings.num_cols);

                wing->wake_rings.x[icurr] = wing->wake_rings.x[iprev];
                wing->wake_rings.y[icurr] = wing->wake_rings.y[iprev];
                wing->wake_rings.z[icurr] = wing->wake_rings.z[iprev];
            }
        }

        add(&point, &wing->position, &point);
        vector_to_mesh(&point, &wing->wake_rings, iwake);
    }

    wing->wake_rings.num_rows = wing->iteration + 1;

    if (wing->wake_rings.num_rows < wing->num_wake_deforming_rows) {
        wing->wake_displacements.num_rows = wing->wake_rings.num_rows;
    } else {
        wing->wake_displacements.num_rows = wing->num_wake_deforming_rows;
    }

    num_points = get_size(&wing->wake_rings);

    Vector3D rotation_negative = {-wing->rotation.x, -wing->rotation.y, -wing->rotation.z};

    assign_rotation(rot_mat, &rotation_negative);

    for (size_t i = 0; i < num_points; i++) {
        mesh_to_vector(&wing->wake_rings, i, &point);
        subtract(&point, &wing->position, &point);
        apply_rotation(rot_mat, &point);
        vector_to_mesh(&point, &wing->wake_rings, i);
    }

    if (wing->wake_rings.num_rows > 1) {
        for (int j = 0; j < wing->num_spanwise_panels; j++) {
            if (wing->wake_rings.num_rows > 2) {
                for (int i = wing->wake_rings.num_rows - 2; i > 0; i--) {
                    icurr = Sub2Ind(i, j, wing->num_spanwise_panels);
                    iprev = Sub2Ind(i - 1, j, wing->num_spanwise_panels);

                    wing->wake_vorticity[icurr] = wing->wake_vorticity[iprev];
                }
            }

            iwing = Sub2Ind(wing->num_chordwise_panels - 1, j, wing->num_spanwise_panels);
            iwake = Sub2Ind(0, j, wing->num_spanwise_panels);

            wing->wake_vorticity[iwake] = wing->bound_vorticity[iwing];
        }
    }
}