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
#include "compute_pressures.h"
#include "assign_corners.h"

#define SAVE_VTK_FILES 1
#define NUM_TIME_STEPS 160
#define SWEEP_ANGLE_LEADING 90.0
#define SWEEP_ANGLE_TRAILING 90.0
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

    FILE *fp = fopen("lift_results.csv", "w");

    if (fp == NULL) {
        fprintf(stderr, "main: failed to open lift_results.csv");

        return 1;
    }

    fprintf(fp, "Time (sec),Lift Coefficient\n");

    double t;
    double dx = ROOT_CHORD / NUM_CHORDWISE_PANELS;
    double dt = dx / FAR_FIELD_VELOCITY / 4.0;
    double dx_wake = 0.3 * FAR_FIELD_VELOCITY * dt;

    const double vel_sq = FAR_FIELD_VELOCITY * FAR_FIELD_VELOCITY;
    const double rho = AIR_DENSITY;

    char file_name[50];

    Wing wing_obj = {.num_chordwise_panels = NUM_CHORDWISE_PANELS,
                     .num_spanwise_panels = NUM_SPANWISE_PANELS,
                     .num_time_steps = NUM_TIME_STEPS,
                     .num_wake_deforming_rows = NUM_WAKE_DEFORMING_ROWS,
                     .naca_m = NACA_M,
                     .naca_p = NACA_P,
                     .semi_span = SEMI_SPAN,
                     .root_chord = ROOT_CHORD,
                     .sweep_angle_leading = SWEEP_ANGLE_LEADING,
                     .sweep_angle_trailing = SWEEP_ANGLE_TRAILING,
                     .angle_of_attack = ANGLE_OF_ATTACK,
                     .cutoff_radius = CUTOFF,
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
        compute_pressures(wing, dt, AIR_DENSITY);

        if (istep && SAVE_VTK_FILES) {
            snprintf(file_name, sizeof(file_name), "wake_rings.vtk.%d", istep);
            write_vtk_file(&wing_obj.wake_rings, file_name);
        }

        current = clock();

        printf("completed in %.0f msec...\n", ((double) (current - last)) * 1000.0 / CLOCKS_PER_SEC);
        fprintf(fp, "%f,%f\n", t, wing->total_lift / (0.5 * rho * vel_sq));

        last = current;
    }

    putchar('\n');
    print_attributes(wing);
    printf("Elapsed Time: %.2f sec\n", ((double) (current - start)) / CLOCKS_PER_SEC);

    fclose(fp);
    free(wing->horizontal_buffer);
    free(wing->pivot_vector);
    free(wing->memory.elements);

    return 0;
}