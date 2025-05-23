#include "wing.h"
#include "sub2ind.h"
#include "get_size.h"
#include "apply_rotation.h"
#include "assign_rotation.h"

void shed_wake(Wing *wing) {
    Vector point;

    double rot_mat[3][3];

    int iwing;
    int iring;
    int iwake;
    int icurr;
    int iprev;
    int num_points;

    if (wing->iteration) {
        num_points = get_size(&wing->wake_rings);

        assign_rotation(rot_mat, &wing->rotation_prev);

        for (int i = 0; i < num_points; i++) {
            point.x = wing->wake_rings.x[i];
            point.y = wing->wake_rings.y[i];
            point.z = wing->wake_rings.z[i];

            apply_rotation(rot_mat, &point);

            wing->wake_rings.x[i] = point.x + wing->position_prev.x;
            wing->wake_rings.y[i] = point.y + wing->position_prev.y;
            wing->wake_rings.z[i] = point.z + wing->position_prev.z;
        }
    }

    assign_rotation(rot_mat, &wing->rotation);

    for (int j = 0; j < wing->bound_rings.num_cols; j++) {
        iring = sub2ind(wing->bound_rings.num_rows - 1, j, wing->bound_rings.num_cols);
        iwake = sub2ind(0, j, wing->wake_rings.num_cols);

        point.x = wing->bound_rings.x[iring];
        point.y = wing->bound_rings.y[iring];
        point.z = wing->bound_rings.z[iring];

        apply_rotation(rot_mat, &point);

        if (wing->iteration) {
            for (int i = wing->iteration; i > 0; i--) {
                icurr = sub2ind(i, j, wing->wake_rings.num_cols);
                iprev = sub2ind(i - 1, j, wing->wake_rings.num_cols);

                wing->wake_rings.x[icurr] = wing->wake_rings.x[iprev];
                wing->wake_rings.y[icurr] = wing->wake_rings.y[iprev];
                wing->wake_rings.z[icurr] = wing->wake_rings.z[iprev];
            }
        }

        wing->wake_rings.x[iwake] = point.x + wing->position.x;
        wing->wake_rings.y[iwake] = point.y + wing->position.y;
        wing->wake_rings.z[iwake] = point.z + wing->position.z;
    }

    wing->wake_rings.num_rows = wing->iteration + 1;

    if (wing->wake_rings.num_rows < wing->num_wake_deforming_rows) {
        wing->wake_displacements.num_rows = wing->wake_rings.num_rows;
    } else {
        wing->wake_displacements.num_rows = wing->num_wake_deforming_rows;
    }

    num_points = get_size(&wing->wake_rings);

    Vector rotation_negative = {-wing->rotation.x, -wing->rotation.y, -wing->rotation.z};

    assign_rotation(rot_mat, &rotation_negative);

    for (int i = 0; i < num_points; i++) {
        point.x = wing->wake_rings.x[i] - wing->position.x;
        point.y = wing->wake_rings.y[i] - wing->position.y;
        point.z = wing->wake_rings.z[i] - wing->position.z;

        apply_rotation(rot_mat, &point);

        wing->wake_rings.x[i] = point.x;
        wing->wake_rings.y[i] = point.y;
        wing->wake_rings.z[i] = point.z;
    }

    if (wing->wake_rings.num_rows > 1) {
        for (int j = 0; j < wing->num_spanwise_panels; j++) {
            if (wing->wake_rings.num_rows > 2) {
                for (int i = wing->wake_rings.num_rows - 2; i > 0; i--) {
                    icurr = sub2ind(i, j, wing->num_spanwise_panels);
                    iprev = sub2ind(i - 1, j, wing->num_spanwise_panels);

                    wing->wake_vorticity[icurr] = wing->wake_vorticity[iprev];
                }
            }

            iwing = sub2ind(wing->num_chordwise_panels - 1, j, wing->num_spanwise_panels);
            iwake = sub2ind(0, j, wing->num_spanwise_panels);

            wing->wake_vorticity[iwake] = wing->bound_vorticity[iwing];
        }
    }
}