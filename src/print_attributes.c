#include <stdio.h>

#include "wing.h"

void print_attributes(Wing *wing) {
    printf("Semi Span: %.2f m\n", wing->semi_span);
    printf("Root Chord: %.2f m\n", wing->root_chord);
    printf("Leading Edge Sweep Angle: %.2f deg\n", wing->sweep_angle_leading);
    printf("Trailing Edge Sweep Angle: %.2f deg\n", wing->sweep_angle_trailing);
    printf("Angle of Attack: %.2f deg\n", wing->angle_of_attack);
    printf("Surface Area: %.2f sq. m\n", wing->surface_area);
    printf("Aspect Ratio: %.2f\n", wing->aspect_ratio);
    printf("Airfoil: NACA %d%dXX\n", wing->naca_m, wing->naca_p);
    printf("Spanwise Panels: %d\n", wing->num_spanwise_panels);
    printf("Chordwise Panels: %d\n", wing->num_chordwise_panels);
    printf("Wake Deforming Rows: %d\n", wing->num_wake_deforming_rows);
}