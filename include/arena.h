#ifndef ARENA_H
#define ARENA_H

#include <assert.h>
#include <stdlib.h>

#include "vector3d.h"

typedef struct DoubleArena {
    size_t num_elements;
    size_t next_free_index;

    double *elements;
} DoubleArena;

typedef struct Vector3DArena {
    size_t num_elements;
    size_t next_free_index;

    Vector3D *elements;
} Vector3DArena;

DoubleArena DoubleArena_construct(size_t num_elements, size_t element_size) {
    DoubleArena arena;

    arena.num_elements = num_elements;
    arena.next_free_index = 0;
    arena.elements = (double *) calloc(num_elements, sizeof(double));

    assert(arena.elements != NULL);

    return arena;
}

Vector3DArena Vector3DArena_construct(size_t num_elements) {
    Vector3DArena arena;

    arena.num_elements = num_elements;
    arena.next_free_index = 0;
    arena.elements = (Vector3D *) calloc(num_elements, sizeof(Vector3D));

    assert(arena.elements != NULL);

    return arena;
}

double *DoubleArena_allocate(DoubleArena *arena, size_t num_elements) {
    assert(arena->next_free_index + num_elements <= arena->num_elements);

    double *ptr = arena->elements + arena->next_free_index;

    arena->next_free_index += num_elements;
    
    return ptr;
}

Vector3D *Vector3DArena_allocate(Vector3DArena *arena, size_t num_elements) {
    assert(arena->next_free_index + num_elements <= arena->num_elements);

    Vector3D *ptr = arena->elements + arena->next_free_index;

    arena->next_free_index += num_elements;
    
    return ptr;
}

#endif