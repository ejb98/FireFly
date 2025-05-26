#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "wing.h"
#include "subtract.h"
#include "compute_magnitude.h"
#include "vector.h"
#include "process.h"
#include "sub2ind.h"
#include "get_size.h"
#include "init_wing.h"
#include "mesh_to_vector.h"
#include "write_vtk_file.h"
#include "print_firefly.h"
#include "print_attributes.h"
#include "compute_mean_chord.h"
#include "compute_between.h"
#include "compute_area.h"
#include "assign_corners.h"

#define NEWTON_TO_POUND 0.224809
#define PASCAL_TO_PSI 0.000145038
#define SAVE_VTK_FILES 1
#define NUM_TIME_STEPS 160
#define SWEEP_ANGLE_LEADING 70.0
#define SWEEP_ANGLE_TRAILING 80.0
#define ANGLE_OF_ATTACK 5.0
#define PITCHING_FREQUENCY 2.0
#define PITCHING_AMPLITUDE 0.0
#define HEAVING_FREQUENCY 2.0
#define HEAVING_AMPLITUDE 0.0
#define FAR_FIELD_VELOCITY 10.0
#define NUM_CHORDWISE_PANELS 4
#define NUM_SPANWISE_PANELS 13
#define NUM_WAKE_DEFORMING_ROWS 160
#define AIR_DENSITY 1.0
#define ROOT_CHORD 1.0
#define SEMI_SPAN 4.0
#define CUTOFF 1e-6
#define NACA_M 6
#define NACA_P 4
#define PI 3.141592654

int main(int argc, char **argv) {
    print_firefly();

    clock_t last;
    clock_t start;
    clock_t current;

    start = clock();

    double t;
    double dx = ROOT_CHORD / NUM_CHORDWISE_PANELS;
    double dt = dx / FAR_FIELD_VELOCITY / 4.0;
    double dx_wake = 0.3 * FAR_FIELD_VELOCITY * dt;
    double lift;

    char file_name[50];

    Wing wing_obj = {.num_chordwise_panels = NUM_CHORDWISE_PANELS,
                     .num_spanwise_panels = NUM_SPANWISE_PANELS,
                     .num_wake_point_rows = NUM_TIME_STEPS,
                     .num_wake_deforming_rows = NUM_WAKE_DEFORMING_ROWS,
                     .naca_m = NACA_M,
                     .naca_p = NACA_P,
                     .semi_span = SEMI_SPAN,
                     .root_chord = ROOT_CHORD,
                     .sweep_angle_leading = SWEEP_ANGLE_LEADING,
                     .sweep_angle_trailing = SWEEP_ANGLE_TRAILING,
                     .angle_of_attack = ANGLE_OF_ATTACK,
                     .cutoff = CUTOFF,
                     .wake_offset = dx_wake};

    if (init_wing(&wing_obj)) {
        fprintf(stderr, "main: init_wing failed to allocate memory");

        return 1;
    }

    Wing *wing = &wing_obj;

    if (SAVE_VTK_FILES) {
        write_vtk_file(&wing_obj.surface_panels, "surface_panels.vtk");
        write_vtk_file(&wing_obj.bound_rings, "bound_rings.vtk");
        write_vtk_file(&wing_obj.control_points, "control_points.vtk");
    }

    last = clock();

    for (int istep = 0; istep < NUM_TIME_STEPS; istep++) {
        printf("Solving Step %d...", istep);

        t = istep * dt;

        wing_obj.position.x = -FAR_FIELD_VELOCITY * t;
        wing_obj.position.y = 0.0;
        wing_obj.position.z = HEAVING_AMPLITUDE * sin(2.0 * PI * HEAVING_FREQUENCY * t);
        wing_obj.rotation.x = 0.0;
        wing_obj.rotation.y = PITCHING_AMPLITUDE * sin(2.0 * PI * PITCHING_FREQUENCY * t) * PI / 180.0;
        wing_obj.rotation.z = 0.0;

        process(wing, dt);

        lift = 0.0;

        if (istep) {
            Mesh *mesh = &wing->surface_panels;
            Vector corners[4];
            Vector corners_local[4];
            Vector normal;
            Vector front;
            double mean_chord;
            double width;
            double area;
            double gamma;
            double gamma_chordwise;
            double gamma_spanwise;
            double dgammadt;
            double dgammadx;
            double dx_left;
            double dx_right;
            double dx_mean;
            double delta_pressure;
            double *gammas = wing->bound_vorticity;
            size_t ipanel;

            for (int j = 0; j < wing_obj.num_spanwise_panels; j++) {
                for (int i = 0; i < wing_obj.num_chordwise_panels; i++) {
                    ipanel = sub2ind(i, j, wing_obj.num_spanwise_panels);
                    assign_corners(mesh, i, j, corners);
                    area = compute_area(corners);
                    mean_chord = compute_mean_chord(corners);
                    subtract(corners + 1, corners, &front);
                    width = compute_magnitude(&front);
                    mesh_to_vector(&wing->normal_vectors, ipanel, &normal);

                    dgammadx = 0.0;
                    for (int ichord = 0; ichord <= i; ichord++) {

                        if (!ichord) {
                            gamma = gammas[sub2ind(ichord, j, mesh->num_cols)];
                        } else {
                            gamma = gammas[sub2ind(ichord, j, mesh->num_cols)] - 
                                    gammas[sub2ind(ichord - 1, j, mesh->num_cols)];
                        }
                        
                        assign_corners(mesh, ichord, j, corners_local);

                        dx_left = corners_local[3].x - corners_local[0].x;
                        dx_right = corners_local[2].x - corners_local[1].x;

                        dx_mean = (dx_left + dx_right) / 2.0;

                        if (!ichord) {
                            dgammadx += 0.5 * gamma * dx_mean;
                        } else {
                            dgammadx += (0.5 * gamma + gammas[sub2ind(ichord - 1, j, mesh->num_cols)]) * dx_mean;
                        }
                    }

                    dgammadt = (dgammadx - wing->bound_vorticity_prev[ipanel]) / dt;

                    if (!i) {
                        gamma_chordwise = gammas[ipanel];
                    } else {
                        gamma_chordwise = gammas[ipanel] - gammas[sub2ind(i - 1, j, mesh->num_cols)];
                    }

                    if (!j) {
                        gamma_spanwise = gammas[ipanel];
                    } else {
                        gamma_chordwise = gammas[ipanel] - gammas[sub2ind(i, j - 1, mesh->num_cols)];
                    }

                    delta_pressure = AIR_DENSITY * (wing->freestream_velocities[ipanel] * gamma_chordwise / mean_chord +
                                                    wing->spanwise_velocities[ipanel] * gamma_spanwise / width + dgammadt);
                    
                    lift += delta_pressure * area * normal.x;

                    wing->bound_vorticity_prev[ipanel] = dgammadx;
                }
            }
        }

        if (istep && SAVE_VTK_FILES) {
            snprintf(file_name, sizeof(file_name), "wake_rings.vtk.%d", istep);
            write_vtk_file(&wing_obj.wake_rings, file_name);
        }

        current = clock();

        printf("completed in %.0f msec...", ((double) (current - last)) * 1000.0 / CLOCKS_PER_SEC);
        printf("Lift = %f lbf\n", lift);

        last = current;
    }

    putchar('\n');
    print_attributes(wing);
    printf("Elapsed Time: %.2f sec\n", ((double) (current - start)) / CLOCKS_PER_SEC);

    free(wing->horizontal_buffer);
    free(wing->pivot_vector);
    free(wing->memory.elements);

    return 0;
}