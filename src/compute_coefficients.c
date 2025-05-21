#include <stdlib.h>

#include "dot.h"
#include "mesh.h"
#include "vector.h"
#include "sub2ind.h"
#include "get_size.h"
#include "mesh_to_vector.h"
#include "induce_by_segment.h"

void compute_coefficients(Mesh *points, Mesh *rings, Mesh *normals,
                          double *a, double *b, Vector *buffer,
                          double cutoff, int mirror, int append) {
    
    int num_points = get_size(points);
    int num_ring_rows = rings->num_rows - 1;
    int num_ring_cols = rings->num_cols - 1;
    int num_rings = num_ring_rows * num_ring_cols;
    int i1, i2, i3, i4;
    int imatrix;
    int iring;
    int ibuff;

    Vector point;
    Vector normal;
    Vector velocity;
    Vector downwash;
    Vector p1, p2, p3, p4;

    for (int ipoint = 0; ipoint < num_points; ipoint++) {
        mesh_to_vector(points, ipoint, &point);
        mesh_to_vector(normals, ipoint, &normal);

        if (mirror) {
            point.y = -point.y;
            normal.y = -normal.y;
        }

        iring = 0;
        for (int i = 0; i < num_ring_rows; i++) {

            ibuff = 0;
            for (int j = 0; j < num_ring_cols; j++) {
                i1 = sub2ind(i, j, rings->num_cols);
                i2 = sub2ind(i, j + 1, rings->num_cols);
                i3 = sub2ind(i + 1, j + 1, rings->num_cols);
                i4 = sub2ind(i + 1, j, rings->num_cols);

                mesh_to_vector(rings, i1, &p1);
                mesh_to_vector(rings, i2, &p2);
                mesh_to_vector(rings, i3, &p3);
                mesh_to_vector(rings, i4, &p4);

                if (i > 0) {
                    buffer[ibuff].x = -buffer[4 * j + 2].x;
                    buffer[ibuff].y = -buffer[4 * j + 2].y;
                    buffer[ibuff].z = -buffer[4 * j + 2].z;

                    ibuff++;
                } else {
                    induce_by_segment(&point, &p1, &p2, buffer + ibuff++, 1.0, cutoff);
                }

                induce_by_segment(&point, &p2, &p3, buffer + ibuff++, 1.0, cutoff);
                induce_by_segment(&point, &p3, &p4, buffer + ibuff++, 1.0, cutoff);

                if (j > 0) {
                    buffer[ibuff].x = -buffer[ibuff - 6].x;
                    buffer[ibuff].y = -buffer[ibuff - 6].y;
                    buffer[ibuff].z = -buffer[ibuff - 6].z;

                    ibuff++;
                } else {
                    induce_by_segment(&point, &p4, &p1, buffer + ibuff++, 1.0, cutoff);
                }

                downwash.x = buffer[ibuff - 3].x + buffer[ibuff - 1].x;
                downwash.y = buffer[ibuff - 3].y + buffer[ibuff - 1].y;
                downwash.z = buffer[ibuff - 3].z + buffer[ibuff - 1].z;

                velocity.x = downwash.x + buffer[ibuff - 4].x + buffer[ibuff - 2].x;
                velocity.y = downwash.y + buffer[ibuff - 4].y + buffer[ibuff - 2].y;
                velocity.z = downwash.z + buffer[ibuff - 4].z + buffer[ibuff - 2].z;

                imatrix = sub2ind(ipoint, iring, num_rings);

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

                iring++;
            }
        }
    }
}