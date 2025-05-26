#include "vector.h"

void divide(Vector *vector, double denominator) {
    vector->x /= denominator;
    vector->y /= denominator;
    vector->z /= denominator;
}