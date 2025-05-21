#include "wing.h"
#include "sub2ind.h"
#include "get_size.h"
#include "apply_rotation.h"
#include "assign_rotation.h"

void shed_wake(Wing *wing) {
    double x, y, z;
    double rot_mat[3][3];

    int iring;
    int iwake;
    int icurr;
    int iprev;
    int num_points;

    if (wing->iteration) {
        num_points = get_size(&wing->wake_rings);

        assign_rotation(rot_mat, wing->pitch_prev, wing->roll_prev, wing->yaw_prev);

        for (int i = 0; i < num_points; i++) {
            apply_rotation(rot_mat, wing->wake_rings.x + i,
                            wing->wake_rings.y + i,
                            wing->wake_rings.z + i);

            wing->wake_rings.x[i] += wing->x_pos_prev;
            wing->wake_rings.y[i] += wing->y_pos_prev;
            wing->wake_rings.z[i] += wing->z_pos_prev;
        }
    }

    assign_rotation(rot_mat, wing->pitch, wing->roll, wing->yaw);

    for (int j = 0; j < wing->bound_rings.num_cols; j++) {
        iring = sub2ind(wing->bound_rings.num_rows - 1, j, wing->bound_rings.num_cols);
        iwake = sub2ind(0, j, wing->wake_rings.num_cols);

        x = wing->bound_rings.x[iring];
        y = wing->bound_rings.y[iring];
        z = wing->bound_rings.z[iring];

        apply_rotation(rot_mat, &x, &y, &z);

        if (wing->iteration) {
            for (int i = wing->iteration; i > 0; i--) {
                icurr = sub2ind(i, j, wing->wake_rings.num_cols);
                iprev = sub2ind(i - 1, j, wing->wake_rings.num_cols);

                wing->wake_rings.x[icurr] = wing->wake_rings.x[iprev];
                wing->wake_rings.y[icurr] = wing->wake_rings.y[iprev];
                wing->wake_rings.z[icurr] = wing->wake_rings.z[iprev];
            }
        }

        wing->wake_rings.x[iwake] = x + wing->x_pos;
        wing->wake_rings.y[iwake] = y + wing->y_pos;
        wing->wake_rings.z[iwake] = z + wing->z_pos;
    }

    wing->wake_rings.num_rows = wing->iteration + 1;
    num_points = get_size(&wing->wake_rings);

    assign_rotation(rot_mat, -wing->pitch, -wing->roll, -wing->yaw);

    for (int i = 0; i < num_points; i++) {
        wing->wake_rings.x[i] -= wing->x_pos;
        wing->wake_rings.y[i] -= wing->y_pos;
        wing->wake_rings.z[i] -= wing->z_pos;

        apply_rotation(rot_mat, wing->wake_rings.x + i,
                        wing->wake_rings.y + i,
                        wing->wake_rings.z + i);
    }
}