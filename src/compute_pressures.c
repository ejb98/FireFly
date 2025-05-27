#include <stdlib.h>

#include "dot.h"
#include "wing.h"
#include "vector.h"
#include "sub2ind.h"
#include "subtract.h"
#include "compute_area.h"
#include "mesh_to_vector.h"
#include "assign_corners.h"
#include "compute_between.h"
#include "compute_magnitude.h"

void compute_pressures(Wing *wing, double delta_time, double rho) {
    if (!wing->iteration) {
        wing->lift = 0.0;

        return;
    }

    double integral;
    double gammasum;
    double gammaij;
    double gammapi;
    double gammapj;
    double normz;
    double gmpj;
    double gmpi;
    double lift;
    double area;
    double dgdt;
    double vi;
    double vj;
    double dx;
    double dy;
    double a;
    double b;
    double h;

    size_t ivortex;
    size_t ivortex_first;
    size_t ivortex_previ;
    size_t ivortex_prevj;

    Vector point;
    Vector front;
    Vector leading;
    Vector corners[4];

    lift = 0.0;
    for (int j = 0; j < wing->num_spanwise_panels; j++) {
        ivortex_first = sub2ind(0, j, wing->num_spanwise_panels);

        gammasum = 0.0;
        for (int i = 0; i < wing->num_chordwise_panels; i++) {
            ivortex = sub2ind(i, j, wing->num_spanwise_panels);

            mesh_to_vector(&wing->control_points, ivortex, &point);
            assign_corners(&wing->surface_panels, i, j, corners);
            subtract(corners + 1, corners, &front);

            normz = wing->normal_vectors.z[ivortex];
            area = compute_area(corners);
            
            dy = compute_magnitude(&front);
            dx = 0.5 * ((corners[3].x - corners[0].x) + (corners[2].x - corners[1].x));
            vi = wing->chordwise_velocities[ivortex];
            vj = wing->spanwise_velocities[ivortex];

            gammaij = wing->bound_vorticity[ivortex];

            if (i) {
                ivortex_previ = sub2ind(i - 1, j, wing->num_spanwise_panels);
                gammapi = wing->bound_vorticity[ivortex_previ];
                gmpi = gammaij - gammapi;

                b = point.x - leading.x;
                h = (b - a) / i;
            } else {
                compute_between(corners, corners + 1, 0.5, &leading);
                gmpi = gammaij;

                a = point.x - leading.x;
            }

            if (j) {
                ivortex_prevj = sub2ind(i, j - 1, wing->num_spanwise_panels);
                gammapj = wing->bound_vorticity[ivortex_prevj];
                gmpj = gammaij - gammapj;
            } else {
                gmpj = gammaij;
            }

            gammasum += 2.0 * gammaij;

            if (i == 0) {
                integral = 0.0;
            } else if (i == 1) {
                integral = 0.25 * h * (gammapi + gammaij);
            } else {
                integral = 0.25 * h * (gammasum - wing->bound_vorticity[ivortex_first] - gammaij);
            }

            dgdt = (integral - wing->vorticity_integral_buffer[ivortex]) / delta_time;

            wing->vorticity_integral_buffer[ivortex] = integral;
            wing->pressures[ivortex] = rho * (vi * gmpi / dx + vj * gmpj / dy + dgdt);

            lift -= wing->pressures[ivortex] * area * normz;
        }
    }

    wing->lift = lift;
}