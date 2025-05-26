#include "vector.h"
#include "subtract.h"
#include "compute_magnitude.h"

double compute_area(Vector *corners) {
    Vector left, right, front;

    subtract(corners + 1, corners, &front);
    subtract(corners + 3, corners, &left);
    subtract(corners + 2, corners + 1, &right);

    double mean_chord = (compute_magnitude(&left) + compute_magnitude(&right)) / 2.0;
    double width = compute_magnitude(&front);

    return width * mean_chord;
}