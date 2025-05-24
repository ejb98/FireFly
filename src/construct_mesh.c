#include <stdio.h>
#include <stdlib.h>

#include "mesh.h"
#include "arena.h"
#include "arena_allocate.h"

Mesh construct_mesh(int num_rows, int num_cols, Arena *arena) {
    size_t num_elements = ((size_t) num_rows) * num_cols;

    Mesh mesh = {num_rows, num_cols};

    mesh.x = arena_allocate(num_elements, arena);
    mesh.y = arena_allocate(num_elements, arena);
    mesh.z = arena_allocate(num_elements, arena);

    return mesh;
}