#include <math.h>

#include "dot.h"
#include "wing.h"
#include "cross.h"
#include "vector.h"
#include "sub2ind.h"

void init_control_points(Wing *wing) {
    double x_left;
    double y_left;
    double z_left;
    double x_right;
    double y_right;
    double z_right;
    double dx_left;
    double dz_left;
    double dx_right;
    double dz_right;
    double magnitude;

    Vector a, b, c;
    size_t i1, i2, i3, i4, ip;

    for (int i = 0; i < wing->control_points.num_rows; i++) {
        for (int j = 0; j < wing->control_points.num_cols; j++) {
            ip = sub2ind(i, j, wing->control_points.num_cols);
            i1 = sub2ind(i, j, wing->surface_panels.num_cols);
            i2 = sub2ind(i, j + 1, wing->surface_panels.num_cols);
            i3 = sub2ind(i + 1, j + 1, wing->surface_panels.num_cols);
            i4 = sub2ind(i + 1, j, wing->surface_panels.num_cols);

            dx_left = wing->surface_panels.x[i4] - wing->surface_panels.x[i1];
            dz_left = wing->surface_panels.z[i4] - wing->surface_panels.z[i1];

            dx_right = wing->surface_panels.x[i3] - wing->surface_panels.x[i2];
            dz_right = wing->surface_panels.z[i3] - wing->surface_panels.z[i2];

            y_left = wing->surface_panels.y[i1];
            x_left = wing->surface_panels.x[i1] + 0.75 * dx_left;
            z_left = wing->surface_panels.z[i1] + 0.75 * dz_left;

            y_right = wing->surface_panels.y[i2];
            x_right = wing->surface_panels.x[i2] + 0.75 * dx_right;
            z_right = wing->surface_panels.z[i2] + 0.75 * dz_right;

            wing->control_points.x[ip] = (x_left + x_right) / 2.0;
            wing->control_points.y[ip] = (y_left + y_right) / 2.0;
            wing->control_points.z[ip] = (z_left + z_right) / 2.0;

            a.x = wing->surface_panels.x[i3] - wing->surface_panels.x[i1];
            a.y = wing->surface_panels.y[i3] - wing->surface_panels.y[i1];
            a.z = wing->surface_panels.z[i3] - wing->surface_panels.z[i1];

            b.x = wing->surface_panels.x[i2] - wing->surface_panels.x[i4];
            b.y = wing->surface_panels.y[i2] - wing->surface_panels.y[i4];
            b.z = wing->surface_panels.z[i2] - wing->surface_panels.z[i4];

            cross(&a, &b, &c);

            magnitude = sqrt(dot(&c, &c));

            wing->normal_vectors.x[ip] = c.x / magnitude;
            wing->normal_vectors.y[ip] = c.y / magnitude;
            wing->normal_vectors.z[ip] = c.z / magnitude;
        }
    }
}