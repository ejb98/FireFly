#include <stdlib.h>
#include <math.h>

#include "wing.h"
#include "vector.h"
#include "sub2ind.h"
#include "mesh_to_vector.h"

void init_control_points(Wing *wing) {
    Vector a, b;
    Vector left;
    Vector back;
    Vector right;
    Vector front;
    Vector normal;
    Vector tangent;
    Vector corners[4];

    Mesh *points = &wing->surface_panels;

    int num_cols = points->num_cols;
    double magnitude;
    
    size_t ipoint;

    for (int i = 0; i < wing->control_points.num_rows; i++) {
        for (int j = 0; j < wing->control_points.num_cols; j++) {
            ipoint = sub2ind(i, j, wing->control_points.num_cols);

            mesh_to_vector(points, sub2ind(i, j, num_cols), corners);
            mesh_to_vector(points, sub2ind(i, j + 1, num_cols), corners + 1);
            mesh_to_vector(points, sub2ind(i + 1, j + 1, num_cols), corners + 2);
            mesh_to_vector(points, sub2ind(i + 1, j, num_cols), corners + 3);

            left.y = corners[0].y;
            left.x = corners[0].x + 0.75 * (corners[3].x - corners[0].x);
            left.z = corners[0].z + 0.75 * (corners[3].z - corners[0].z);

            right.y = corners[1].y;
            right.x = corners[1].x + 0.75 * (corners[2].x - corners[1].x);
            right.z = corners[1].z + 0.75 * (corners[2].z - corners[1].z);

            back.x = (corners[3].x + corners[2].x) / 2.0;
            back.y = (corners[3].y + corners[2].y) / 2.0;
            back.z = (corners[3].z + corners[2].z) / 2.0;

            front.x = (corners[0].x + corners[1].x) / 2.0;
            front.y = (corners[0].y + corners[1].y) / 2.0;
            front.z = (corners[0].z + corners[1].z) / 2.0;

            wing->control_points.x[ipoint] = (left.x + right.x) / 2.0;
            wing->control_points.y[ipoint] = (left.y + right.y) / 2.0;
            wing->control_points.z[ipoint] = (left.z + right.z) / 2.0;

            a.x = corners[2].x - corners[0].x;
            a.y = corners[2].y - corners[0].y;
            a.z = corners[2].z - corners[0].z;

            b.x = corners[1].x - corners[3].x;
            b.y = corners[1].y - corners[3].y;
            b.z = corners[1].z - corners[3].z;

            normal.x = a.y * b.z - a.z * b.y;
            normal.y = a.z * b.x - a.x * b.z;
            normal.z = a.x * b.y - a.y * b.x;

            magnitude = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

            wing->normal_vectors.x[ipoint] = normal.x / magnitude;
            wing->normal_vectors.y[ipoint] = normal.y / magnitude;
            wing->normal_vectors.z[ipoint] = normal.z / magnitude;

            tangent.x = front.x - back.x;
            tangent.y = front.y - back.y;
            tangent.z = front.z - back.z;

            magnitude = sqrt(tangent.x * tangent.x + tangent.y * tangent.y + tangent.z * tangent.z);

            wing->tangent_vectors.x[ipoint] = tangent.x / magnitude;
            wing->tangent_vectors.y[ipoint] = tangent.y / magnitude;
            wing->tangent_vectors.z[ipoint] = tangent.z / magnitude;
        }
    }
}