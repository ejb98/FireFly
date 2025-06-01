#include <stdio.h>
#include <math.h>
#include <time.h>

#include "wing.h"
#include "geometry.h"
#include "constants.h"
#include "print_firefly.h"

#define NACA_M 6
#define NACA_P 4
#define CUTOFF 1e-10
#define SEMI_SPAN 4.0
#define ROOT_CHORD 1.0
#define AIR_DENSITY 1.0
#define LEADING_SWEEP 90.0
#define TRAILING_SWEEP 90.0
#define NUM_TIME_STEPS 160
#define ANGLE_OF_ATTACK 5.0
#define HEAVING_FREQUENCY 0.0
#define HEAVING_AMPLITUDE 0.0
#define PITCHING_FREQUENCY 0.0
#define PITCHING_AMPLITUDE 0.0
#define FAR_FIELD_VELOCITY 10.0
#define NUM_SPANWISE_PANELS 13
#define NUM_CHORDWISE_PANELS 4

int main(int argc, char **argv) {
    PrintFireFly();

    clock_t last;
    clock_t start;
    clock_t current;

    start = clock();

    double t;
    double dx = ROOT_CHORD / NUM_CHORDWISE_PANELS;
    double dt = dx / FAR_FIELD_VELOCITY / 4.0;
    double dx_wake = 0.3 * FAR_FIELD_VELOCITY * dt;
    double elapsed_time;

    WingProperties props = {.naca_m = NACA_M, .naca_p = NACA_P,
                            .num_chordwise_panels = NUM_CHORDWISE_PANELS,
                            .num_spanwise_panels = NUM_SPANWISE_PANELS,
                            .semi_span = SEMI_SPAN, .root_chord = ROOT_CHORD,
                            .angle_of_attack = ANGLE_OF_ATTACK,
                            .leading_edge_sweep_angle = LEADING_SWEEP,
                            .trailing_edge_sweep_angle = TRAILING_SWEEP};

    Wing *wing = Wing_Construct(&props, NUM_TIME_STEPS, dx_wake, CUTOFF);
    
    Wing_WritePoints2VTK(wing, SURFACE_POINTS, "results\\");                     
    Wing_WritePoints2VTK(wing, CONTROL_POINTS, "results\\");                     
    Wing_WritePoints2VTK(wing, BOUND_RING_POINTS, "results\\");                     

    last = clock();

    for (int istep = 0; istep < NUM_TIME_STEPS; istep++) {
        printf("Solving Step %d...", istep);

        t = istep * dt;

        wing->position.x = -FAR_FIELD_VELOCITY * t;
        wing->position.z = HEAVING_AMPLITUDE * sin(2.0 * PI * HEAVING_FREQUENCY * t);
        wing->rotation.y = PITCHING_AMPLITUDE * sin(2.0 * PI * PITCHING_FREQUENCY * t) * PI / 180.0;

        Wing_Process(wing, dt);

        current = clock();
        elapsed_time = ((double) (current - last)) * 1000.0 / CLOCKS_PER_SEC;

        printf("completed in %.0f msec\n", elapsed_time);

        last = current;
    }

    putchar('\n');
    Wing_PrintAttributes(wing);

    elapsed_time = ((double) (current - start)) / CLOCKS_PER_SEC;
    printf("Elapsed Time: %.2f sec\n", elapsed_time);

    Wing_Deallocate(wing);

    return 0;
}