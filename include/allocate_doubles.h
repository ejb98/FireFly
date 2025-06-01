#ifndef ALLOCATE_DOUBLES_H
#define ALLOCATE_DOUBLES_H

#include <stdio.h>
#include <stdlib.h>

double *AllocateDoubles(size_t num_elements) {
    double *ptr = (double *) calloc(num_elements, sizeof(double));

    if (ptr == NULL) {
        fprintf(stderr, "AllocateDoubles: memory allocation for %zu elements returned NULL", num_elements);
    }

    return ptr;
}

#endif