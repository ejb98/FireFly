#include <stdlib.h>

size_t sub2ind(int i, int j, int num_cols) {
    return ((size_t) i) * num_cols + j;
}