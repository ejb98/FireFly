#include "dot.h"
#include "mesh.h"
#include "vector.h"
#include "sub2ind.h"
#include "get_size.h"
#include "mesh_to_vector.h"
#include "induce_by_ring.h"

void compute_coefficients(Mesh *points, Mesh *normals, Mesh *rings,
                          double *a, double *b, double cutoff) {
    int iring;
    int imatrix;
    int i0, i1, i2, i3;
    int num_rings = (rings->num_rows - 1) * (rings->num_cols - 1);
    int num_control_points = get_size(points);

    Vector u_induced;
    Vector v_induced;
    Vector point;
    Vector normal;
    Vector corners[4];

    for (int mirror = 0; mirror < 2; mirror++) {
        for (int ipoint = 0; ipoint < num_control_points; ipoint++) {
            mesh_to_vector(normals, ipoint, &normal);
            mesh_to_vector(points, ipoint, &point);

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

                    u_induced = induce_by_ring(&point, corners, &v_induced, 1.0, cutoff);

                    a[imatrix] += dot(&v_induced, &normal);
                    b[imatrix] += dot(&u_induced, &normal);

                    iring++;
                }
            }
        }
    }
}