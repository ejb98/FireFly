#include <stdlib.h>
#include <math.h>

#include "wing.h"
#include "cross.h"
#include "divide.h"
#include "vector.h"
#include "sub2ind.h"
#include "subtract.h"
#include "assign_corners.h"
#include "vector_to_mesh.h"
#include "mesh_to_vector.h"
#include "compute_normal.h"
#include "compute_between.h"
#include "compute_direction.h"
#include "compute_magnitude.h"

void init_control_points(Wing *wing) {
    Vector normal;
    Vector control;
    Vector tangent_spanwise;
    Vector tangent_chordwise;
    Vector corners[4];
    Vector middle_back;
    Vector middle_front;
    Vector middle_left;
    Vector middle_right;
    Vector control_left;
    Vector control_right;
    
    size_t ipoint;

    for (int i = 0; i < wing->control_points.num_rows; i++) {
        for (int j = 0; j < wing->control_points.num_cols; j++) {
            ipoint = sub2ind(i, j, wing->control_points.num_cols);

            assign_corners(&wing->surface_panels, i, j, corners);

            compute_between(corners, corners + 3, 0.75, &control_left);
            compute_between(corners, corners + 3, 0.5, &middle_left);
            compute_between(corners + 1, corners + 2, 0.75, &control_right);
            compute_between(corners + 1, corners + 2, 0.5, &middle_right);
            compute_between(corners + 3, corners + 2, 0.5, &middle_back);
            compute_between(corners, corners + 1, 0.5, &middle_front);
            compute_between(&control_left, &control_right, 0.5, &control);
            compute_direction(&middle_front, &middle_back, &tangent_chordwise);
            compute_direction(&middle_left, &middle_right, &tangent_spanwise);
            compute_normal(corners, &normal);

            vector_to_mesh(&control, &wing->control_points, ipoint);
            vector_to_mesh(&normal, &wing->normal_vectors, ipoint);
            vector_to_mesh(&tangent_chordwise, &wing->tangent_vectors_chordwise, ipoint);
            vector_to_mesh(&tangent_spanwise, &wing->tangent_vectors_spanwise, ipoint);
        }
    }
}