#include <stdlib.h>
#include <stdio.h>

#include "wing.h"
#include "mesh.h"
#include "vector.h"
#include "sub2ind.h"
#include "geometry.h"
#include "get_size.h"
#include "mesh_to_vector.h"
#include "induce_by_ring.h"

void compute_coefficients(Wing *wing, Geometry inducing_rings) {
    Vector u_induced;
    Vector v_induced;
    Vector point;
    Vector normal;
    Vector corners[4];

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

    size_t iring;
    size_t imatrix;
    size_t i0, i1, i2, i3;
    size_t num_control_points = get_size(&wing->control_points);
    size_t num_rings = ((size_t) (rings->num_rows - 1)) * (rings->num_cols - 1);

    for (int mirror = 0; mirror < 2; mirror++) {
        for (size_t ipoint = 0; ipoint < num_control_points; ipoint++) {
            mesh_to_vector(&wing->normal_vectors, ipoint, &normal);
            mesh_to_vector(&wing->control_points, ipoint, &point);

            if (mirror) {
                point.y = -point.y;
                normal.y = -normal.y;
            }

            iring = 0;
            for (int i = 0; i < rings->num_rows - 1; i++) {
                for (int j = 0; j < rings->num_cols - 1; j++) {
                    i0 = sub2ind(i, j, rings->num_cols);
                    i1 = sub2ind(i, j + 1, rings->num_cols);
                    i2 = sub2ind(i + 1, j + 1, rings->num_cols);
                    i3 = sub2ind(i + 1, j, rings->num_cols);

                    imatrix = sub2ind(ipoint, iring, num_rings);

                    mesh_to_vector(rings, i0, corners);
                    mesh_to_vector(rings, i1, corners + 1);
                    mesh_to_vector(rings, i2, corners + 2);
                    mesh_to_vector(rings, i3, corners + 3);

                    u_induced = induce_by_ring(&point, corners, &v_induced, 1.0, wing->cutoff);

                    a_matrix[imatrix] += v_induced.x * normal.x + 
                                         v_induced.y * normal.y + 
                                         v_induced.z * normal.z;

                    b_matrix[imatrix] += u_induced.x * normal.x + 
                                         u_induced.y * normal.y + 
                                         u_induced.z * normal.z;

                    iring++;
                }
            }
        }
    }
}