#include "vector.h"
#include "subtract.h"
#include "compute_magnitude.h"
#include "compute_mean_chord.h"

double compute_area(Vector *corners) {
    Vector front;

    subtract(corners + 1, corners, &front);

    double mean_chord = compute_mean_chord(corners);
    double width = compute_magnitude(&front);

    return width * mean_chord;
}