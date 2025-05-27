#include <stdlib.h>
#include <stdio.h>

#include "dot.h"
#include "add.h"
#include "wing.h"
#include "mesh.h"
#include "vector.h"
#include "geometry.h"
#include "get_size.h"
#include "mesh_to_vector.h"
#include "induce_velocities.h"

void compute_coefficients(Wing *wing, Geometry inducing_rings) {
    Vector point;
    Vector normal;
    Vector w_induced;
    Vector v_induced;
    Vector vertical_buffer;
    Vector vel_norm[4];

    Mesh *rings;

    double *a_matrix;
    double *b_matrix;

    if (inducing_rings == BOUND_RINGS) {
        rings = &wing->bound_rings;
        a_matrix = wing->a_wing_on_wing;
        b_matrix = wing->b_wing_on_wing;
    } else if (inducing_rings == WAKE_RINGS) {
        rings = &wing->wake_rings;
        a_matrix = wing->a_wake_on_wing;
        b_matrix = wing->b_wake_on_wing;
    } else {
        fprintf(stderr, "compute_coefficients: invalid inducing geometry type");

        return;
    }

    int num_rows = rings->num_rows - 1;
    int num_cols = rings->num_cols - 1;

    size_t iring;
    size_t imatrix;
    size_t num_control_points = get_size(&wing->control_points);
    size_t num_rings = ((size_t) num_rows) * num_cols;

    for (int mirror = 0; mirror < 2; mirror++) {
        for (size_t ipoint = 0; ipoint < num_control_points; ipoint++) {
            mesh_to_vector(&wing->normal_vectors, ipoint, &normal);
            mesh_to_vector(&wing->control_points, ipoint, &point);

            if (mirror) {
                point.y = -point.y;
                normal.y = -normal.y;
            }

            iring = 0;

            for (int i = 0; i < num_rows; i++) {
                for (int j = 0; j < num_cols; j++) {
                    induce_velocities(&point, rings, i, j, &vertical_buffer, 
                                      wing->horizontal_buffer, vel_norm, wing->cutoff);
                                      
                    imatrix = ipoint * num_rings + iring;

                    add(vel_norm + 1, vel_norm + 3, &w_induced);
                    add(vel_norm, vel_norm + 2, &v_induced);
                    add(&v_induced, &w_induced, &v_induced);

                    a_matrix[imatrix] += dot(&v_induced, &normal);
                    b_matrix[imatrix] += dot(&w_induced, &normal);

                    iring++;
                }
            }
        }
    }
}