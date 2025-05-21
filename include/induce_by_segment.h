#ifndef INDUCE_BY_SEGMENT_H
#define INDUCE_BY_SEGMENT_H

#include "vector.h"

void induce_by_segment(Vector *p, Vector *p1, Vector *p2,
                       Vector *velocity, double gamma, double cutoff);
#endif