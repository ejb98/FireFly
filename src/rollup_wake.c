#include <stdlib.h>

#include "wing.h"
#include "mesh.h"
#include "vector.h"
#include "sub2ind.h"
#include "get_size.h"
#include "mesh_to_vector.h"
#include "induce_by_ring.h"

void rollup_wake(Wing *wing, double delta_time) {
    size_t iring;
    size_t i0, i1, i2, i3;
    size_t num_points = get_size(&wing->wake_rings);
    // size_t num_points_trunc = ((size_t) wing->num_wake_deforming_rows) * wing->wake_rings.num_cols;

    // if (num_points > num_points_trunc) {
    //     num_points = num_points_trunc;
    // }

    int rollup_y;

    Vector point;
    Vector v_induced;
    Vector corners[4];

    Mesh *rings;
    Mesh *points = &wing->wake_rings;
    Mesh *displacements = &wing->wake_displacements;
    Mesh *ring_meshes[2] = {&wing->bound_rings, &wing->wake_rings};

    double cutoff = wing->cutoff;
    double *gammas = wing->wake_vorticity;

    for (int imesh = 0; imesh < 2; imesh++) {
        rings = ring_meshes[imesh];

        for (int mirror = 0; mirror < 2; mirror++) {
            for (size_t ipoint = 0; ipoint < num_points; ipoint++) {
                rollup_y = ipoint % points->num_cols;

                mesh_to_vector(points, ipoint, &point);

                if (mirror) {
                    point.y = -point.y;
                }

                iring = 0;
                for (int i = 0; i < rings->num_rows - 1; i++) {
                    for (int j = 0; j < rings->num_cols - 1; j++) {
                        i0 = sub2ind(i, j, rings->num_cols);
                        i1 = sub2ind(i, j + 1, rings->num_cols);
                        i2 = sub2ind(i + 1, j + 1, rings->num_cols);
                        i3 = sub2ind(i + 1, j, rings->num_cols);

                        mesh_to_vector(rings, i0, corners);
                        mesh_to_vector(rings, i1, corners + 1);
                        mesh_to_vector(rings, i2, corners + 2);
                        mesh_to_vector(rings, i3, corners + 3);

                        induce_by_ring(&point, corners, &v_induced, gammas[iring], cutoff);

                        displacements->x[ipoint] += v_induced.x * delta_time;
                        displacements->z[ipoint] += v_induced.z * delta_time;

                        if (rollup_y) {
                            displacements->y[ipoint] += v_induced.y * delta_time;
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