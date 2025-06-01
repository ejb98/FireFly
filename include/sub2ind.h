#ifndef SUB2IND_H
#define SUB2IND_H

#include <stdlib.h>

size_t Sub2Ind(int i, int j, int num_cols) {
    return ((size_t) i) * num_cols + j;
}

#endif