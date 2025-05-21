#include "dot.h"
#include "wing.h"
#include "vector.h"
#include "sub2ind.h"
#include "shed_wake.h"
#include "mesh_to_vector.h"
#include "apply_rotation.h"
#include "assign_rotation.h"

void process(Wing *wing, double delta_time) {
    int ipoint;

    double x_prev, y_prev, z_prev;
    double x_curr, y_curr, z_curr;
    double rot_mat[3][3];
    double rot_mat_prev[3][3];

    wing->iteration++;

    Vector normal;
    Vector velocity;

    if (wing->iteration) {
        assign_rotation(rot_mat, wing->pitch, wing->roll, wing->yaw);
        assign_rotation(rot_mat_prev, wing->pitch_prev, wing->roll_prev, wing->yaw_prev);
    }

    for (int j = 0; j < wing->num_spanwise_panels; j++) {
        for (int i = 0; i < wing->num_chordwise_panels; i++) {
            ipoint = sub2ind(i, j, wing->num_spanwise_panels);

            if (!wing->iteration) {
                wing->normal_velocities[ipoint] = 0.0;
            } else {
                x_curr = wing->control_points.x[ipoint];
                y_curr = wing->control_points.y[ipoint];
                z_curr = wing->control_points.z[ipoint];

                x_prev = x_curr;
                y_prev = y_curr;
                z_prev = z_curr;

                apply_rotation(rot_mat, &x_curr, &y_curr, &z_curr);
                apply_rotation(rot_mat_prev, &x_prev, &y_prev, &z_prev);

                x_curr += wing->x_pos;
                y_curr += wing->y_pos;
                z_curr += wing->z_pos;

                x_prev += wing->x_pos_prev;
                y_prev += wing->y_pos_prev;
                z_prev += wing->z_pos_prev;
                
                mesh_to_vector(&wing->normal_vectors, ipoint, &normal);

                velocity.x = (x_curr - x_prev) / delta_time;
                velocity.y = (y_curr - y_prev) / delta_time;
                velocity.z = (z_curr - z_prev) / delta_time;

                wing->normal_velocities[ipoint] = dot(&velocity, &normal);
            }
        }
    }

    shed_wake(wing);

    wing->x_pos_prev = wing->x_pos;
    wing->y_pos_prev = wing->y_pos;
    wing->z_pos_prev = wing->z_pos;

    wing->yaw_prev = wing->yaw;
    wing->roll_prev = wing->roll;
    wing->pitch_prev = wing->pitch;
}