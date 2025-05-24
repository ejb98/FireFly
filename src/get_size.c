#include <stdlib.h>

#include "mesh.h"

size_t get_size(const Mesh* mesh) {
    return ((size_t) mesh->num_rows) * mesh->num_cols;
}