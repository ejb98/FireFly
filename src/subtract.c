#include "vector.h"

void subtract(const Vector *a, const Vector *b, Vector *result) {
    result->x = a->x - b->x;
    result->y = a->y - b->y;
    result->z = a->z - b->z;
}