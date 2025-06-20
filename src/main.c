#include <stdio.h>
#include <math.h>
#include <time.h>

#include "wing.h"
#include "geometry.h"
#include "constants.h"
#include "simulation.h"
#include "print_firefly.h"

#define NACA_M 6
#define NACA_P 4
#define CUTOFF 1e-4
#define SEMI_SPAN 4.0
#define FILE_PATH "results"
#define ROOT_CHORD 1.0
#define AIR_DENSITY 1.0
#define LEADING_SWEEP 70.0
#define TRAILING_SWEEP 80.0
#define NUM_TIME_STEPS 160
#define ANGLE_OF_ATTACK 5.0
#define HEAVING_FREQUENCY 2.0
#define HEAVING_AMPLITUDE 0.5
#define PITCHING_FREQUENCY 2.0
#define PITCHING_AMPLITUDE 5.0
#define FAR_FIELD_VELOCITY 10.0
#define NUM_SPANWISE_PANELS 13
#define NUM_CHORDWISE_PANELS 4

int main(int argc, char **argv) {
    PrintFireFly();

    clock_t start = clock();

    double dx = ROOT_CHORD / NUM_CHORDWISE_PANELS;
    double dt = dx / FAR_FIELD_VELOCITY / 4.0;
    double dx_wake = 0.3 * FAR_FIELD_VELOCITY * dt;
    double t = -dt;

    char file_path[] = FILE_PATH;

    Wing wing = {.naca_m = NACA_M, .naca_p = NACA_P,
                 .num_chordwise_panels = NUM_CHORDWISE_PANELS,
                 .num_spanwise_panels = NUM_SPANWISE_PANELS,
                 .semi_span = SEMI_SPAN, .root_chord = ROOT_CHORD,
                 .angle_of_attack = ANGLE_OF_ATTACK,
                 .leading_sweep_angle = LEADING_SWEEP,
                 .trailing_sweep_angle = TRAILING_SWEEP};

    Simulation *sim = Simulation_Init(&wing, NUM_TIME_STEPS, dt, AIR_DENSITY, dx_wake, CUTOFF, file_path);

    while (!sim->is_complete) {
        t += dt;

        wing.position.x = -FAR_FIELD_VELOCITY * t;
        wing.position.z = HEAVING_AMPLITUDE * sin(2.0 * PI * HEAVING_FREQUENCY * t);
        wing.rotation.y = PITCHING_AMPLITUDE * sin(2.0 * PI * PITCHING_FREQUENCY * t) * PI / 180.0;
        
        Simulation_Process(sim);
    }

    putchar('\n');
    Wing_Print(&wing);

    double elapsed_time = ((double) (clock() - start)) / CLOCKS_PER_SEC;

    printf("Elapsed Time: %.2f sec\n", elapsed_time);
    Simulation_Deallocate(sim);

    return 0;
}