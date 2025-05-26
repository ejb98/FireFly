#include "vector.h"

void multiply(Vector *vector, double multiplier) {
    vector->x *= multiplier;
    vector->y *= multiplier;
    vector->z *= multiplier;
}