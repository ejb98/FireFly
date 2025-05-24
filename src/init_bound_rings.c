#include <stdlib.h>

#include "wing.h"
#include "sub2ind.h"

void init_bound_rings(Wing *wing) {
    double dx;
    double dz;

    size_t index_curr;
    size_t index_next;

    for (int j = 0; j < wing->bound_rings.num_cols; j++) {
        for (int i = 0; i < wing->bound_rings.num_rows - 1; i++) {
            index_curr = sub2ind(i, j, wing->bound_rings.num_cols);
            index_next = sub2ind(i + 1, j, wing->bound_rings.num_cols);

            dx = wing->surface_panels.x[index_next] - wing->surface_panels.x[index_curr];
            dz = wing->surface_panels.z[index_next] - wing->surface_panels.z[index_curr];

            wing->bound_rings.y[index_curr] = wing->surface_panels.y[index_curr];
            wing->bound_rings.z[index_curr] = wing->surface_panels.z[index_curr] + dz / 4.0;
            wing->bound_rings.x[index_curr] = wing->surface_panels.x[index_curr] + dx / 4.0;

            if (i == wing->bound_rings.num_rows - 2) {
                wing->bound_rings.y[index_next] = wing->surface_panels.y[index_next];
                wing->bound_rings.z[index_next] = wing->surface_panels.z[index_next];
                wing->bound_rings.x[index_next] = wing->surface_panels.x[index_next] + wing->wake_offset;
            }
        }
    }
}