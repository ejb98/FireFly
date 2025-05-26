#include "add.h"
#include "vector.h"
#include "subtract.h"
#include "multiply.h"

void compute_between(const Vector *start, const Vector *end, double fraction, Vector *result) {
    subtract(end, start, result);
    multiply(result, fraction);
    add(start, result, result);
}