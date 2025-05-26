#include "vector.h"

double dot(const Vector *a, const Vector *b) {
    return a->x * b->x + a->y * b->y + a->z * b->z;
}