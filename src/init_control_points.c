#include <stdlib.h>
#include <math.h>

#include "wing.h"
#include "vector.h"
#include "sub2ind.h"

void init_control_points(Wing *wing) {
    Vector left;
    Vector right;
    Vector a, b, c;

    double dx_left;
    double dz_left;
    double dx_right;
    double dz_right;
    double magnitude;

    size_t i0, i1, i2, i3, ip;

    for (int i = 0; i < wing->control_points.num_rows; i++) {
        for (int j = 0; j < wing->control_points.num_cols; j++) {
            ip = sub2ind(i, j, wing->control_points.num_cols);
            i0 = sub2ind(i, j, wing->surface_panels.num_cols);
            i1 = sub2ind(i, j + 1, wing->surface_panels.num_cols);
            i2 = sub2ind(i + 1, j + 1, wing->surface_panels.num_cols);
            i3 = sub2ind(i + 1, j, wing->surface_panels.num_cols);

            dx_left = wing->surface_panels.x[i3] - wing->surface_panels.x[i0];
            dz_left = wing->surface_panels.z[i3] - wing->surface_panels.z[i0];

            dx_right = wing->surface_panels.x[i2] - wing->surface_panels.x[i1];
            dz_right = wing->surface_panels.z[i2] - wing->surface_panels.z[i1];

            left.y = wing->surface_panels.y[i0];
            left.x = wing->surface_panels.x[i0] + 0.75 * dx_left;
            left.z = wing->surface_panels.z[i0] + 0.75 * dz_left;

            right.y = wing->surface_panels.y[i1];
            right.x = wing->surface_panels.x[i1] + 0.75 * dx_right;
            right.z = wing->surface_panels.z[i1] + 0.75 * dz_right;

            wing->control_points.x[ip] = (left.x + right.x) / 2.0;
            wing->control_points.y[ip] = (left.y + right.y) / 2.0;
            wing->control_points.z[ip] = (left.z + right.z) / 2.0;

            a.x = wing->surface_panels.x[i2] - wing->surface_panels.x[i0];
            a.y = wing->surface_panels.y[i2] - wing->surface_panels.y[i0];
            a.z = wing->surface_panels.z[i2] - wing->surface_panels.z[i0];

            b.x = wing->surface_panels.x[i1] - wing->surface_panels.x[i3];
            b.y = wing->surface_panels.y[i1] - wing->surface_panels.y[i3];
            b.z = wing->surface_panels.z[i1] - wing->surface_panels.z[i3];

            c.x = a.y * b.z - a.z * b.y;
            c.y = a.z * b.x - a.x * b.z;
            c.z = a.x * b.y - a.y * b.x;

            magnitude = sqrt(c.x * c.x + c.y * c.y + c.z * c.z);

            wing->normal_vectors.x[ip] = c.x / magnitude;
            wing->normal_vectors.y[ip] = c.y / magnitude;
            wing->normal_vectors.z[ip] = c.z / magnitude;
        }
    }
}