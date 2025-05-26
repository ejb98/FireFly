#include "vector.h"
#include "subtract.h"
#include "compute_magnitude.h"

double compute_mean_chord(Vector *corners) {
    Vector left, right;

    subtract(corners + 3, corners, &left);
    subtract(corners + 2, corners + 1, &right);

    return 0.5 * (compute_magnitude(&left) + compute_magnitude(&right));
}