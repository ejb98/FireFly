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

    double rot_mat[3][3];
    double rot_mat_prev[3][3];

    wing->iteration++;

    Vector curr;
    Vector prev;
    Vector normal;
    Vector velocity;

    if (wing->iteration) {
        assign_rotation(rot_mat, &wing->rotation);
        assign_rotation(rot_mat_prev, &wing->rotation_prev);
    }

    for (int j = 0; j < wing->num_spanwise_panels; j++) {
        for (int i = 0; i < wing->num_chordwise_panels; i++) {
            ipoint = sub2ind(i, j, wing->num_spanwise_panels);

            if (!wing->iteration) {
                wing->normal_velocities[ipoint] = 0.0;
            } else {
                curr.x = wing->control_points.x[ipoint];
                curr.y = wing->control_points.y[ipoint];
                curr.z = wing->control_points.z[ipoint];

                prev.x = curr.x;
                prev.y = curr.y;
                prev.z = curr.z;

                apply_rotation(rot_mat, &curr);
                apply_rotation(rot_mat_prev, &prev);

                curr.x += wing->position.x;
                curr.y += wing->position.y;
                curr.z += wing->position.z;

                prev.x += wing->position_prev.x;
                prev.y += wing->position_prev.y;
                prev.z += wing->position_prev.z;
                
                mesh_to_vector(&wing->normal_vectors, ipoint, &normal);

                velocity.x = (curr.x - prev.x) / delta_time;
                velocity.y = (curr.y - prev.y) / delta_time;
                velocity.z = (curr.z - prev.z) / delta_time;

                wing->normal_velocities[ipoint] = dot(&velocity, &normal);
            }
        }
    }

    shed_wake(wing);

    wing->position_prev = wing->position;
    wing->rotation_prev = wing->rotation;
}