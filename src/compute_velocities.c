#include <stdlib.h>

#include "add.h"
#include "dot.h"
#include "wing.h"
#include "vector.h"
#include "divide.h"
#include "sub2ind.h"
#include "subtract.h"
#include "mesh_to_vector.h"
#include "vector_to_mesh.h"
#include "apply_rotation.h"
#include "assign_rotation.h"

void compute_velocities(Simulation *wing, double delta_time) {
    size_t ipoint;

    double rot_mat[3][3];
    double rot_mat_prev[3][3];

    Vector3D curr;
    Vector3D prev;
    Vector3D normal;
    Vector3D velocity;

    assign_rotation(rot_mat, &wing->rotation);
    assign_rotation(rot_mat_prev, &wing->rotation_prev);

    for (int j = 0; j < wing->nspanwise_panel; j++) {
        for (int i = 0; i < wing->nchordwise_panels; i++) {
            ipoint = Sub2Ind(i, j, wing->nspanwise_panel);

            mesh_to_vector(&wing->control_points, ipoint, &curr);

            prev.x = curr.x;
            prev.y = curr.y;
            prev.z = curr.z;

            apply_rotation(rot_mat, &curr);
            apply_rotation(rot_mat_prev, &prev);
            add(&curr, &wing->position, &curr);
            add(&prev, &wing->position_prev, &prev);
            
            mesh_to_vector(&wing->normals, ipoint, &normal);

            subtract(&curr, &prev, &velocity);
            divide(&velocity, delta_time);

            vector_to_mesh(&velocity, &wing->kinematic_velocities, ipoint);

            wing->normal_velocities[ipoint] = dot(&velocity, &normal);
        }
    }
}