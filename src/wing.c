#include <stdio.h>

#include "wing.h"

void Wing_Print(const Wing *wing) {
    printf("Semi Span: %.2f m\n", wing->semi_span);
    printf("Root Chord: %.2f m\n", wing->root_chord);
    printf("Leading Edge Sweep Angle: %.2f deg\n", wing->leading_sweep_angle);
    printf("Trailing Edge Sweep Angle: %.2f deg\n", wing->trailing_sweep_angle);
    printf("Angle of Attack: %.2f deg\n", wing->angle_of_attack);
    printf("Airfoil: NACA %d%dXX\n", wing->naca_m, wing->naca_p);
    printf("Spanwise Panels: %d\n", wing->num_spanwise_panels);
    printf("Chordwise Panels: %d\n", wing->num_chordwise_panels);
}