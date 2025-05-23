#ifndef INDUCE_BY_SEGMENT_H
#define INDUCE_BY_SEGMENT_H

#include "vector.h"

void induce_by_segment(Vector *point, Vector *point1, Vector *point2,
                       Vector *velocity, double gamma, double cutoff);
#endif