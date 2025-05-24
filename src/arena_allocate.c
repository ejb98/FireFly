#include <assert.h>
#include <stdlib.h>

#include "arena.h"

double *arena_allocate(size_t num_elements, Arena *arena) {
    assert(arena->next_free_index + num_elements <= arena->num_elements);

    double *ptr = arena->elements + arena->next_free_index;

    arena->next_free_index += num_elements;

    return ptr;
}