#include <stdlib.h>

#include "dot.h"
#include "wing.h"
#include "vector.h"
#include "sub2ind.h"
#include "subtract.h"
#include "compute_area.h"
#include "mesh_to_vector.h"
#include "assign_corners.h"
#include "compute_magnitude.h"

void compute_pressures(Wing *wing, double delta_time, double rho) {
    double sigma;
    double sigma1;
    double gamaij;
    double normz;
    double width;
    double dfdt;
    double lift;
    double vinf;
    double area;
    double dx;

    size_t ivortex;
    size_t ivortex_previ;

    Vector front;
    Vector corners[4];

    lift = 0.0;

    if (!wing->iteration) {
        wing->lift = lift;

        return;
    }

    for (int j = 0; j < wing->num_spanwise_panels; j++) {
        sigma = 0.0;
        sigma1 = 0.0;

        for (int i = 0; i < wing->num_chordwise_panels; i++) {
            ivortex = sub2ind(i, j, wing->num_spanwise_panels);

            vinf = -wing->chordwise_velocities[ivortex];

            assign_corners(&wing->surface_panels, i, j, corners);
            subtract(corners + 1, corners, &front);

            area = compute_area(corners);
            width = compute_magnitude(&front);
            normz = wing->normal_vectors.z[ivortex];

            dx = 0.5 * ((corners[3].x - corners[0].x) + (corners[2].x - corners[1].x));

            if (i) {
                ivortex_previ = sub2ind(i - 1, j, wing->num_spanwise_panels);
                gamaij = wing->bound_vorticity[ivortex] -
                         wing->bound_vorticity[ivortex_previ];
            } else {
                gamaij = wing->bound_vorticity[ivortex];
            }

            sigma1 = (0.5 * gamaij + sigma) * dx;
            sigma = wing->bound_vorticity[ivortex];

            dfdt = (sigma1 - wing->vorticity_integral_buffer[ivortex]) / delta_time;

            wing->vorticity_integral_buffer[ivortex] = sigma1;
            wing->pressures[ivortex] = rho * (vinf * gamaij + dfdt) * width;

            lift += wing->pressures[ivortex] * area * normz;
        }
    }

    wing->lift = lift;
}