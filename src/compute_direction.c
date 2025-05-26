#include "vector.h"
#include "divide.h"
#include "subtract.h"
#include "compute_magnitude.h"

void compute_direction(const Vector *start, const Vector *end, Vector *direction) {
    subtract(end, start, direction);
    divide(direction, compute_magnitude(direction));
}