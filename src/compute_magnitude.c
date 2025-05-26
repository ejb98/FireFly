#include <math.h>

#include "dot.h"
#include "vector.h"

double compute_magnitude(const Vector *vector) {
    return sqrt(dot(vector, vector));
}