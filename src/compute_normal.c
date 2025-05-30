#include "cross.h"
#include "vector.h"
#include "divide.h"
#include "subtract.h"
#include "compute_magnitude.h"

void compute_normal(Vector3D *corners, Vector3D *result) {
    Vector3D a, b;

    subtract(corners + 2, corners, &a);
    subtract(corners + 1, corners + 3, &b);
    cross(&a, &b, result);
    divide(result, compute_magnitude(result));
}