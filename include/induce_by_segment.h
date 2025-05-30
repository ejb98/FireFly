#ifndef INDUCE_BY_SEGMENT_H
#define INDUCE_BY_SEGMENT_H

#include "vector.h"

void induce_by_segment(Vector3D *point, Vector3D *point1, Vector3D *point2, Vector3D *velocity_normalized, double cutoff);

#endif