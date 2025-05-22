#include <stdlib.h>

#include "dot.h"
#include "wing.h"
#include "mesh.h"
#include "sub2ind.h"
#include "induce_velocity.h"

void rollup_wake(Wing *wing, double delta_time) {
    int mirror;
    int append;

    Mesh *rings;
    Mesh *points = &wing->wake_rings;

    double *gammas;

    for (int i = 0; i < 4; i++) {
        mirror = i % 2;
        append = i;

        if (i < 2) {
            rings = &wing->bound_rings; 
            gammas = wing->bound_vorticity;
        } else {
            rings = &wing->wake_rings;
            gammas = wing->wake_vorticity;
        }

        induce_velocity(points, rings, &wing->wake_velocities, NULL, gammas,
                        wing->velocity_buffer, wing->cutoff, mirror, append);
    }

    int index;

    for (int i = 0; i < wing->wake_rings.num_rows; i++) {
        for (int j = 0; j < wing->wake_rings.num_cols; j++) {
            index = sub2ind(i, j, wing->wake_rings.num_cols);

            wing->wake_rings.x[index] += wing->wake_velocities.x[index] * delta_time;
            wing->wake_rings.z[index] += wing->wake_velocities.z[index] * delta_time;

            if (j > 0) {
                wing->wake_rings.y[index] += wing->wake_velocities.y[index] * delta_time;
            }
        }
    }
}