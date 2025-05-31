#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <stddef.h>

typedef struct Vector3D {
    double x;
    double y;
    double z;
} Vector3D;

void Vector3D_add(const Vector3D *vector_a, const Vector3D *vector_b, Vector3D *sum);
void Vector3D_subtract(const Vector3D *vector_a, const Vector3D *vector_b, Vector3D *difference);
void Vector3D_divide(Vector3D *vector, double divisor);
void Vector3D_multiply(Vector3D *vector, double multiplier);
void Vector3D_fill_rotation_matrix(const Vector3D *rotation, double rotation_matrix[3][3]);
void Vector3D_rotate(Vector3D *vector, double rotation_matrix[3][3]);
void Vector3D_cross(const Vector3D *vector_a, const Vector3D *vector_b, Vector3D *product);
void Vector3D_between(const Vector3D *start, const Vector3D *end, double fraction, Vector3D *result);
double Vector3D_dot(const Vector3D *vector_a, const Vector3D *vector_b);
double Vector3D_magnitude(const Vector3D *vector);
Vector3D *Vector3D_malloc(size_t num_elements);

#endif