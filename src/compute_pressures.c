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

    double gamma_previ;
    double gamma_prevj;
    double integral;
    double partial;
    double gamma;
    double normz;
    double lift;
    double area;
    double vi;
    double vj;
    double dx;
    double dy;

    size_t ivortex;

    Vector point;
    Vector front;
    Vector corners[4];

    lift = 0.0;
    for (int j = 0; j < wing->num_spanwise_panels; j++) {

        integral = 0.0;
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

            gamma = wing->bound_vorticity[ivortex];

            if (i > 1) {
                gamma_previ = wing->bound_vorticity[sub2ind(i - 1, j, wing->num_spanwise_panels)];
            } else {
                gamma_previ = 0.0;
            }

            if (j) {
                gamma_prevj = wing->bound_vorticity[sub2ind(i, j - 1, wing->num_spanwise_panels)];
            } else {
                gamma_prevj = 0.0;
            }

            integral += gamma * dx;
            partial = (integral - wing->vorticity_integral_buffer[ivortex]) / delta_time;

            wing->vorticity_integral_buffer[ivortex] = integral;
            wing->pressures[ivortex] = rho * (vi * (gamma - gamma_previ) / dx + 
                                              vj * (gamma - gamma_prevj) / dy + partial);

            lift -= wing->pressures[ivortex] * area * normz;
        }
    }

    wing->lift = lift;
}