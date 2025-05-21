#include <math.h>

#include "wing.h"
#include "sub2ind.h"
#include "apply_rotation.h"
#include "assign_rotation.h"

void init_surface_panels(Wing *wing) {
    int index;

    const double pi = 3.141592654;

    double x_norm;
    double z_norm;
    double local_chord;
    double offset_leading;
    double offset_trailing;
    double x_tip_leading;
    double x_tip_trailing;
    double constant;
    double rotation_matrix[3][3];
    double m = wing->naca_m / 100.0;
    double p = wing->naca_p / 10.0;
    double lambda_leading = wing->sweep_angle_leading * pi / 180.0;
    double lambda_trailing = wing->sweep_angle_trailing * pi / 180.0;
    double tan_leading = tan(pi / 2.0 - lambda_leading);
    double tan_trailing = tan(pi / 2.0 - lambda_trailing);

    Vector point;
    Vector rotation_aoa = {0.0, wing->angle_of_attack * pi / 180.0, 0.0};

    assign_rotation(rotation_matrix, &rotation_aoa);

    for (int j = 0; j < wing->surface_panels.num_cols; j++) {
        point.y = wing->semi_span * j / wing->num_spanwise_panels;

        for (int i = 0; i < wing->surface_panels.num_rows; i++) {
            x_norm = ((double) i) / wing->num_chordwise_panels;

            constant = 2 * p * x_norm - x_norm * x_norm;

            if (x_norm < p && wing->naca_p != 0) {
                z_norm = (m / (p * p)) * constant;
            } else {
                z_norm = (m / ((1 - p) * (1 - p))) * (1 - 2 * p + constant);
            }

            offset_leading = point.y * tan_leading;
            offset_trailing = point.y * tan_trailing;

            local_chord = wing->root_chord + offset_trailing - offset_leading;

            point.x = offset_leading + x_norm * local_chord;
            point.z = local_chord * z_norm;

            if (j == wing->surface_panels.num_cols - 1) {
                if (i == 0) {
                    x_tip_leading = point.x;
                }

                if (i == wing->surface_panels.num_rows - 1) {
                    x_tip_trailing = point.x;
                }
            }

            apply_rotation(rotation_matrix, &point);

            index = sub2ind(i, j, wing->surface_panels.num_cols);

            wing->surface_panels.x[index] = point.x;
            wing->surface_panels.y[index] = point.y;
            wing->surface_panels.z[index] = point.z;
        }
    }

    double tip_chord = x_tip_trailing - x_tip_leading;

    wing->surface_area = wing->semi_span * (wing->root_chord + tip_chord) / 2.0;
    wing->aspect_ratio = 2.0 * wing->semi_span * wing->semi_span / wing->surface_area;
}