#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "wing.h"
#include "process.h"
#include "init_wing.h"

#define NUM_TIME_STEPS 160
#define SWEEP_ANGLE_LEADING 70.0
#define SWEEP_ANGLE_TRAILING 80.0
#define ANGLE_OF_ATTACK 5.0
#define PITCHING_FREQUENCY 2.0
#define PITCHING_AMPLITUDE 5.0
#define HEAVING_FREQUENCY 2.0
#define HEAVING_AMPLITUDE 0.0
#define FAR_FIELD_VELOCITY 10.0
#define NUM_CHORDWISE_PANELS 4
#define NUM_SPANWISE_PANELS 4
#define NUM_WAKE_DEFORMING_ROWS 5
#define AIR_DENSITY 1.0
#define ROOT_CHORD 1.0
#define SEMI_SPAN 4.0
#define CUTOFF 0.0001
#define NACA_M 6
#define NACA_P 4
#define PI 3.141592654

int main(int argc, char **argv) {
    clock_t start, end;

    start = clock();

    double t;
    double dx = ROOT_CHORD / NUM_CHORDWISE_PANELS;
    double dt = dx / FAR_FIELD_VELOCITY / 4.0;
    double dx_wake = 0.3 * FAR_FIELD_VELOCITY * dt;
    double elapsed_time;

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
        perror("main: init_wing failed to allocate all memory, exiting now");
        return 1;
    }

    Wing *wing = &wing_obj;

    for (int i = 0; i < NUM_TIME_STEPS; i++) {
        t = i * dt;

        wing_obj.x_pos = -FAR_FIELD_VELOCITY * t;
        wing_obj.y_pos = 0;
        wing_obj.z_pos = HEAVING_AMPLITUDE * sin(2.0 * PI * HEAVING_FREQUENCY * t);
        wing_obj.yaw = 0;
        wing_obj.roll = 0;
        wing_obj.pitch = PITCHING_AMPLITUDE * sin(2.0 * PI * PITCHING_FREQUENCY * t) * PI / 180.0;

        process(wing, dt);

        if (i) {
            // TODO: call the solve(wing) function for i > 0 since vorticity is already zero at i == 0
        }

        // TODO: call rollup_wake(wing) function to assign vorticity and perturb the wake
    }

    end = clock();
    elapsed_time = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Elapsed Time: %.3f sec\n", elapsed_time);
    printf("Wing Size: %zu bytes\n", sizeof(Wing));

    return 0;
}