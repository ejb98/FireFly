#include <stdlib.h>

#include "add.h"
#include "dot.h"
#include "wing.h"
#include "vector.h"
#include "divide.h"
#include "sub2ind.h"
#include "subtract.h"
#include "mesh_to_vector.h"
#include "apply_rotation.h"
#include "assign_rotation.h"

void compute_velocities(Wing *wing, double delta_time) {
    size_t ipoint;

    double rot_mat[3][3];
    double rot_mat_prev[3][3];

    Vector curr;
    Vector prev;
    Vector normal;
    Vector tangent_chordwise;
    Vector tangent_spanwise;
    Vector velocity;

    assign_rotation(rot_mat, &wing->rotation);
    assign_rotation(rot_mat_prev, &wing->rotation_prev);

    for (int j = 0; j < wing->num_spanwise_panels; j++) {
        for (int i = 0; i < wing->num_chordwise_panels; i++) {
            ipoint = sub2ind(i, j, wing->num_spanwise_panels);

            mesh_to_vector(&wing->control_points, ipoint, &curr);

            prev.x = curr.x;
            prev.y = curr.y;
            prev.z = curr.z;

            apply_rotation(rot_mat, &curr);
            apply_rotation(rot_mat_prev, &prev);
            add(&curr, &wing->position, &curr);
            add(&prev, &wing->position_prev, &prev);
            
            mesh_to_vector(&wing->normal_vectors, ipoint, &normal);
            mesh_to_vector(&wing->tangent_vectors_chordwise, ipoint, &tangent_chordwise);
            mesh_to_vector(&wing->tangent_vectors_spanwise, ipoint, &tangent_spanwise);

            subtract(&curr, &prev, &velocity);
            divide(&velocity, delta_time);

            wing->normal_velocities[ipoint] = dot(&velocity, &normal);
            wing->chordwise_velocities[ipoint] = dot(&velocity, &tangent_chordwise);
            wing->spanwise_velocities[ipoint] = dot(&velocity, &tangent_spanwise);
        }
    }
}