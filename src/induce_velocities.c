#include "mesh.h"
#include "vector.h"
#include "sub2ind.h"
#include "mesh_to_vector.h"
#include "induce_by_segment.h"

void induce_velocities(Vector *point, Mesh *rings, int i, int j, Vector *vertical_buffer,
                       Vector *horizontal_buffer, Vector *velocities, double cutoff) {
    Vector corners[4];

    mesh_to_vector(rings, sub2ind(i, j, rings->num_cols), corners);
    mesh_to_vector(rings, sub2ind(i, j + 1, rings->num_cols), corners + 1);
    mesh_to_vector(rings, sub2ind(i + 1, j + 1, rings->num_cols), corners + 2);
    mesh_to_vector(rings, sub2ind(i + 1, j, rings->num_cols), corners + 3);

    if (i > 0) {
        velocities[0].x = -horizontal_buffer[j].x;
        velocities[0].y = -horizontal_buffer[j].y;
        velocities[0].z = -horizontal_buffer[j].z;
    } else {
        induce_by_segment(point, corners, corners + 1, velocities, cutoff);
    }

    induce_by_segment(point, corners + 1, corners + 2, velocities + 1, cutoff);
    induce_by_segment(point, corners + 2, corners + 3, velocities + 2, cutoff);

    if (j > 0) {
        velocities[3].x = -vertical_buffer->x;
        velocities[3].y = -vertical_buffer->y;
        velocities[3].z = -vertical_buffer->z;
    } else {
        induce_by_segment(point, corners + 3, corners, velocities + 3, cutoff);
    }

    *vertical_buffer = velocities[1];
    horizontal_buffer[j] = velocities[2];
}
