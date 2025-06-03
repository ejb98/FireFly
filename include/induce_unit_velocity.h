#ifndef INDUCE_UNIT_VELOCITY_H
#define INDUCE_UNIT_VELOCITY_H

#include "vector3d.h"

void InduceUnitVelocity(const Vector3D *point, const Vector3D *point1, 
                        const Vector3D *point2, Vector3D *unit_velocity, double cutoff);

#endif