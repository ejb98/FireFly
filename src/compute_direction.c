#include "vector.h"
#include "divide.h"
#include "subtract.h"
#include "compute_magnitude.h"

void compute_direction(const Vector3D *start, const Vector3D *end, Vector3D *direction) {
    subtract(end, start, direction);
    divide(direction, compute_magnitude(direction));
}