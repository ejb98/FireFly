#include <stdlib.h>

#include "wing.h"
#include "mesh.h"
#include "vector.h"
#include "sub2ind.h"
#include "get_size.h"
#include "mesh_to_vector.h"
#include "induce_by_segment.h"

void rollup_wake(Wing *wing, double delta_time) {
    int num_rows;
    int num_cols;

    size_t iring;
    size_t i0, i1, i2, i3;
    size_t num_points = get_size(&wing->wake_rings);
    size_t num_points_trunc = ((size_t) wing->num_wake_deforming_rows) * wing->wake_rings.num_cols;

    if (num_points > num_points_trunc) {
        num_points = num_points_trunc;
    }

    int rollup_y;

    Vector point;
    Vector vel_induced;
    Vector vertical_buffer;
    Vector corners[4];
    Vector vel_norm[4];

    Mesh *rings;
    Mesh *points = &wing->wake_rings;
    Mesh *displacements = &wing->wake_displacements;
    Mesh *ring_meshes[2] = {&wing->bound_rings, &wing->wake_rings};

    double gamma;
    double cutoff = wing->cutoff;
    double *gammas = wing->wake_vorticity;

    for (int imesh = 0; imesh < 2; imesh++) {
        rings = ring_meshes[imesh];

        if (imesh) {
            num_rows = rings->num_rows - 2;
        } else {
            num_rows = rings->num_rows - 1;
        }
        
        num_cols = rings->num_cols - 1;

        for (int mirror = 0; mirror < 2; mirror++) {
            for (size_t ipoint = 0; ipoint < num_points; ipoint++) {
                rollup_y = ipoint % points->num_cols;

                mesh_to_vector(points, ipoint, &point);

                if (mirror) {
                    point.y = -point.y;
                }

                iring = 0;
                for (int i = 0; i < num_rows; i++) {
                    for (int j = 0; j < num_cols; j++) {
                        gamma = gammas[iring];

                        i0 = sub2ind(i, j, rings->num_cols);
                        i1 = sub2ind(i, j + 1, rings->num_cols);
                        i2 = sub2ind(i + 1, j + 1, rings->num_cols);
                        i3 = sub2ind(i + 1, j, rings->num_cols);

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

                        vel_induced.x = (vel_norm[0].x + vel_norm[1].x + vel_norm[2].x + vel_norm[3].x) * gamma;
                        vel_induced.z = (vel_norm[0].z + vel_norm[1].z + vel_norm[2].z + vel_norm[3].z) * gamma;

                        displacements->x[ipoint] += vel_induced.x * delta_time;
                        displacements->z[ipoint] += vel_induced.z * delta_time;

                        if (rollup_y) {
                            vel_induced.y = (vel_norm[0].y + vel_norm[1].y + vel_norm[2].y + vel_norm[3].y) * gamma;

                            displacements->y[ipoint] += vel_induced.y * delta_time;
                        }

                        iring++;
                    }
                }
            }
        }
    }

    for (size_t ipoint = 0; ipoint < num_points; ipoint++) {
        points->x[ipoint] += displacements->x[ipoint];
        points->y[ipoint] += displacements->y[ipoint];
        points->z[ipoint] += displacements->z[ipoint];
    }
}