#include <stdlib.h>

#include "add.h"
#include "dot.h"
#include "wing.h"
#include "vector.h"
#include "sub2ind.h"
#include "subtract.h"
#include "mesh_to_vector.h"
#include "assign_corners.h"

void compute_pressures(Wing *wing, double delta_time, double rho) {
    if (!wing->iteration) {
        wing->lift = 0.0;

        return;
    }

    double gamma_previ;
    double gamma_prevj;
    double derivative;
    double gamma;
    double normz;
    double lift;
    double dx;
    double dy;

    size_t ivortex;

    Vector front;
    Vector corners[4];
    Vector velocity;
    Vector tangent_spanwise;
    Vector tangent_chordwise;
    Vector kinematic_velocity;
    Vector wake_induced_velocity;

    lift = 0.0;
    for (int j = 0; j < wing->num_spanwise_panels; j++) {

        for (int i = 0; i < wing->num_chordwise_panels; i++) {
            ivortex = sub2ind(i, j, wing->num_spanwise_panels);

            assign_corners(&wing->surface_panels, i, j, corners);
            subtract(corners + 1, corners, &front);

            normz = wing->normal_vectors.z[ivortex];
            
            dy = corners[1].y - corners[0].y;
            dx = 0.5 * ((corners[3].x - corners[0].x) + (corners[2].x - corners[1].x));
            
            gamma = wing->bound_vorticity[ivortex];

            if (i) {
                gamma_previ = wing->bound_vorticity[sub2ind(i - 1, j, wing->num_spanwise_panels)];
            } else {
                gamma_previ = 0.0;
            }

            if (j) {
                gamma_prevj = wing->bound_vorticity[sub2ind(i, j - 1, wing->num_spanwise_panels)];
            } else {
                gamma_prevj = 0.0;
            }

            derivative = (gamma - wing->vorticity_prev[ivortex]) / delta_time;

            wing->vorticity_prev[ivortex] = gamma;

            mesh_to_vector(&wing->tangent_vectors_chordwise, ivortex, &tangent_chordwise);
            mesh_to_vector(&wing->tangent_vectors_spanwise, ivortex, &tangent_spanwise);
            mesh_to_vector(&wing->wake_induced_velocities, ivortex, &wake_induced_velocity);
            mesh_to_vector(&wing->kinematic_velocities, ivortex, &kinematic_velocity);

            add(&kinematic_velocity, &wake_induced_velocity, &velocity);

            wing->pressures[ivortex] = rho * (dot(&velocity, &tangent_chordwise) * (gamma - gamma_previ) / dx +
                                              dot(&velocity, &tangent_spanwise) * (gamma - gamma_prevj) / dy +
                                              derivative);

            lift -= dx * dy * wing->pressures[ivortex] * normz;
        }
    }

    wing->lift = lift;
}