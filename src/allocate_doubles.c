#include <stdio.h>
#include <stdlib.h>

double *AllocateDoubles(size_t num_elements) {
    double *ptr = (double *) calloc(num_elements, sizeof(double));

    if (ptr == NULL) {
        fprintf(stderr, "AllocateDoubles: calloc returned NULL");

        return NULL;
    }

    return ptr;
}