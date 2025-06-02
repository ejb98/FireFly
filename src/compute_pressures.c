#include <stdlib.h>
#include <stdio.h>

#include "add.h"
#include "dot.h"
#include "wing.h"
#include "vector.h"
#include "sub2ind.h"
#include "subtract.h"
#include "mesh_to_vector.h"
#include "assign_corners.h"
#include "compute_magnitude.h"

void compute_pressures(Simulation *wing, double delta_time, double rho) {
    if (!wing->iteration) {
        wing->total_lift = 0.0;

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

    Vector3D front;
    Vector3D left;
    Vector3D right;
    Vector3D corners[4];
    Vector3D velocity;
    Vector3D tangent_spanwise;
    Vector3D tangent_chordwise;
    Vector3D kinematic_velocity;
    Vector3D wake_induced_velocity;

    lift = 0.0;
    for (int j = 0; j < wing->nspanwise_panel; j++) {

        for (int i = 0; i < wing->nchordwise_panels; i++) {
            ivortex = Sub2Ind(i, j, wing->nspanwise_panel);

            assign_corners(&wing->surface_panels, i, j, corners);
            subtract(corners + 1, corners, &front);
            subtract(corners + 3, corners, &left);
            subtract(corners + 2, corners + 1, &right);

            normz = wing->normals.z[ivortex];
            
            dx = 0.5 * (compute_magnitude(&left) + compute_magnitude(&right));
            dy = compute_magnitude(&front);

            if (i == 0 && j == 0) {
                printf("...(%f, %f)...", dx, dy);
            }

            gamma = wing->bound_vorticity[ivortex];

            if (i) {
                gamma_previ = wing->bound_vorticity[Sub2Ind(i - 1, j, wing->nspanwise_panel)];
            } else {
                gamma_previ = 0.0;
            }

            if (j) {
                gamma_prevj = wing->bound_vorticity[Sub2Ind(i, j - 1, wing->nspanwise_panel)];
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

    wing->total_lift = lift;
}