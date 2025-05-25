#include <stdlib.h>
#include <stdio.h>

#include "wing.h"
#include "mesh.h"
#include "vector.h"
#include "sub2ind.h"
#include "geometry.h"
#include "get_size.h"
#include "mesh_to_vector.h"
#include "induce_by_segment.h"

void compute_coefficients(Wing *wing, Geometry inducing_rings) {
    Vector point;
    Vector normal;
    Vector u_induced;
    Vector v_induced;
    Vector vertical_buffer;
    Vector corners[4];
    Vector vel_norm[4];

    Mesh *rings;

    double cutoff = wing->cutoff;

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
    size_t i0, i1, i2, i3;
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
                    i0 = sub2ind(i, j, rings->num_cols);
                    i1 = sub2ind(i, j + 1, rings->num_cols);
                    i2 = sub2ind(i + 1, j + 1, rings->num_cols);
                    i3 = sub2ind(i + 1, j, rings->num_cols);

                    imatrix = sub2ind(ipoint, iring, num_rings);

                    mesh_to_vector(rings, i0, corners);
                    mesh_to_vector(rings, i1, corners + 1);
                    mesh_to_vector(rings, i2, corners + 2);
                    mesh_to_vector(rings, i3, corners + 3);

                    if (i > 0) {
                        vel_norm[0].x = -wing->horizontal_buffer[j].x;
                        vel_norm[0].y = -wing->horizontal_buffer[j].y;
                        vel_norm[0].z = -wing->horizontal_buffer[j].z;
                    } else {
                        induce_by_segment(&point, corners, corners + 1, vel_norm, cutoff);
                    }

                    induce_by_segment(&point, corners + 1, corners + 2, vel_norm + 1, cutoff);
                    induce_by_segment(&point, corners + 2, corners + 3, vel_norm + 2, cutoff);

                    if (j > 0) {
                        vel_norm[3].x = -vertical_buffer.x;
                        vel_norm[3].y = -vertical_buffer.y;
                        vel_norm[3].z = -vertical_buffer.z;
                    } else {
                        induce_by_segment(&point, corners + 3, corners, vel_norm + 3, cutoff);
                    }

                    if (j < num_cols - 1) {
                        vertical_buffer = vel_norm[1];
                    }

                    if (i < num_rows - 1) {
                        wing->horizontal_buffer[j] = vel_norm[2];
                    }

                    u_induced.x = vel_norm[1].x + vel_norm[3].x;
                    u_induced.y = vel_norm[1].y + vel_norm[3].y;
                    u_induced.z = vel_norm[1].z + vel_norm[3].z;

                    v_induced.x = vel_norm[0].x + vel_norm[2].x + u_induced.x;
                    v_induced.y = vel_norm[0].y + vel_norm[2].y + u_induced.y;
                    v_induced.z = vel_norm[0].z + vel_norm[2].z + u_induced.z;

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