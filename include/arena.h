#ifndef ARENA_H
#define ARENA_H

#include <stdlib.h>

typedef struct Arena {
    size_t num_elements;
    size_t next_free_index;

    double *elements;
} Arena;

#endif