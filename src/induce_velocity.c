#include <stdlib.h>

#include "dot.h"
#include "mesh.h"
#include "vector.h"
#include "sub2ind.h"
#include "get_size.h"
#include "mesh_to_vector.h"
#include "induce_by_segment.h"

void induce_velocity(Mesh *points,
                     Mesh *rings,
                     Mesh *normals_or_velocities,
                     double *a,
                     double *b_or_gammas,
                     Vector *velocity_buffer,
                     double cutoff,
                     int mirror,
                     int append) {

    double gamma;

    double *b;
    double *gammas;

    int compute_coefficients = a != NULL;
    
    int num_ring_rows = rings->num_rows - 1;
    int num_ring_cols = rings->num_cols - 1;
    int num_rings = num_ring_rows * num_ring_cols;
    int num_points = compute_coefficients ? get_size(points)
                                          : get_size(normals_or_velocities);


    int i0, i1, i2, i3;
    int imatrix;
    int iring;
    int ibuff;

    Vector point;
    Vector normal;
    Vector velocity;
    Vector downwash;
    Vector p0, p1, p2, p3;

    Mesh *normals;
    Mesh *velocities;

    for (int ipoint = 0; ipoint < num_points; ipoint++) {
        mesh_to_vector(points, ipoint, &point);

        if (mirror) {
            point.y = -point.y;
        }

        if (compute_coefficients) {
            b = b_or_gammas;
            normals = normals_or_velocities;

            mesh_to_vector(normals, ipoint, &normal);

            if (mirror) {
                normal.y = -normal.y;
            }
        } else {
            gammas = b_or_gammas;
            velocities = normals_or_velocities;
        }

        iring = 0;
        for (int i = 0; i < num_ring_rows; i++) {

            ibuff = 0;
            for (int j = 0; j < num_ring_cols; j++) {
                gamma = compute_coefficients ? 1.0 : gammas[iring];

                i0 = sub2ind(i, j, rings->num_cols);
                i1 = sub2ind(i, j + 1, rings->num_cols);
                i2 = sub2ind(i + 1, j + 1, rings->num_cols);
                i3 = sub2ind(i + 1, j, rings->num_cols);

                mesh_to_vector(rings, i0, &p0);
                mesh_to_vector(rings, i1, &p1);
                mesh_to_vector(rings, i2, &p2);
                mesh_to_vector(rings, i3, &p3);

                if (i > 0) {
                    velocity_buffer[ibuff].x = -velocity_buffer[4 * j + 2].x;
                    velocity_buffer[ibuff].y = -velocity_buffer[4 * j + 2].y;
                    velocity_buffer[ibuff].z = -velocity_buffer[4 * j + 2].z;

                    ibuff++;
                } else {
                    induce_by_segment(&point, &p0, &p1, velocity_buffer + ibuff++, gamma, cutoff);
                }

                induce_by_segment(&point, &p1, &p2, velocity_buffer + ibuff++, gamma, cutoff);
                induce_by_segment(&point, &p2, &p3, velocity_buffer + ibuff++, gamma, cutoff);

                if (j > 0) {
                    velocity_buffer[ibuff].x = -velocity_buffer[ibuff - 6].x;
                    velocity_buffer[ibuff].y = -velocity_buffer[ibuff - 6].y;
                    velocity_buffer[ibuff].z = -velocity_buffer[ibuff - 6].z;

                    ibuff++;
                } else {
                    induce_by_segment(&point, &p3, &p0, velocity_buffer + ibuff++, gamma, cutoff);
                }

                downwash.x = velocity_buffer[ibuff - 3].x + velocity_buffer[ibuff - 1].x;
                downwash.y = velocity_buffer[ibuff - 3].y + velocity_buffer[ibuff - 1].y;
                downwash.z = velocity_buffer[ibuff - 3].z + velocity_buffer[ibuff - 1].z;

                velocity.x = downwash.x + velocity_buffer[ibuff - 4].x + velocity_buffer[ibuff - 2].x;
                velocity.y = downwash.y + velocity_buffer[ibuff - 4].y + velocity_buffer[ibuff - 2].y;
                velocity.z = downwash.z + velocity_buffer[ibuff - 4].z + velocity_buffer[ibuff - 2].z;

                imatrix = sub2ind(ipoint, iring, num_rings);

                if (compute_coefficients) {
                    if (append) {
                        a[imatrix] += dot(&velocity, &normal);
                    } else {
                        a[imatrix] = dot(&velocity, &normal);
                    }

                    if (b != NULL) {
                        if (append) {
                            b[imatrix] += dot(&downwash, &normal);
                        } else {
                            b[imatrix] = dot(&downwash, &normal);
                        }
                    }
                } else {
                    if (append) {
                        velocities->x[ipoint] += velocity.x;
                        velocities->y[ipoint] += velocity.y;
                        velocities->z[ipoint] += velocity.z;
                    } else {
                        velocities->x[ipoint] = velocity.x;
                        velocities->y[ipoint] = velocity.y;
                        velocities->z[ipoint] = velocity.z;
                    }
                }

                iring++;
            }
        }
    }
}