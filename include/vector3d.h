#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <stddef.h>

typedef struct Vector3D {
    double x;
    double y;
    double z;
} Vector3D;

void Vector3D_Add(const Vector3D *vector_a, const Vector3D *vector_b, Vector3D *sum);
void Vector3D_Subtract(const Vector3D *vector_a, const Vector3D *vector_b, Vector3D *difference);
void Vector3D_Divide(Vector3D *vector, double divisor);
void Vector3D_Normalize(Vector3D *vector);
void Vector3D_Multiply(Vector3D *vector, double multiplier);
void Vector3D_Rotate(Vector3D *vector, const double (*rotation_matrix)[3]);
void Vector3D_Cross(const Vector3D *vector_a, const Vector3D *vector_b, Vector3D *product);
void Vector3D_Lerp(const Vector3D *start, const Vector3D *end, double fraction, Vector3D *result);
void Vector3D_Direction(const Vector3D *start, const Vector3D *end, Vector3D *direction);
double Vector3D_Distance(const Vector3D *start, const Vector3D *end);
double Vector3D_Dot(const Vector3D *vector_a, const Vector3D *vector_b);
double Vector3D_Magnitude(const Vector3D *vector);
Vector3D *Vector3D_Allocate(size_t num_elements);

#endif