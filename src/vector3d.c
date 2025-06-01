#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "vector3d.h"

void Vector3D_Add(const Vector3D *vector_a, const Vector3D *vector_b, Vector3D *sum) {
    sum->x = vector_a->x + vector_b->x;
    sum->y = vector_a->y + vector_b->y;
    sum->z = vector_a->z + vector_b->z;
}

void Vector3D_Subtract(const Vector3D *vector_a, const Vector3D *vector_b, Vector3D *difference) {
    difference->x = vector_a->x - vector_b->x;
    difference->y = vector_a->y - vector_b->y;
    difference->z = vector_a->z - vector_b->z;
}

void Vector3D_Multiply(Vector3D *vector, double multiplier) {
    vector->x *= multiplier;
    vector->y *= multiplier;
    vector->z *= multiplier;
}

void Vector3D_Divide(Vector3D *vector, double divisor) {
    vector->x /= divisor;
    vector->y /= divisor;
    vector->z /= divisor;
}

void Vector3D_Rotate(Vector3D *vector, const double (*rotation_matrix)[3]) {
    vector->x = rotation_matrix[0][0] * (vector->x) + 
                rotation_matrix[0][1] * (vector->y) + 
                rotation_matrix[0][2] * (vector->z);

    vector->y = rotation_matrix[1][0] * (vector->x) + 
                rotation_matrix[1][1] * (vector->y) + 
                rotation_matrix[1][2] * (vector->z);

    vector->z = rotation_matrix[2][0] * (vector->x) + 
                rotation_matrix[2][1] * (vector->y) + 
                rotation_matrix[2][2] * (vector->z);
}

void Vector3D_Cross(const Vector3D *vector_a, const Vector3D *vector_b, Vector3D *product) {
    product->x = vector_a->y * vector_b->z - vector_a->z * vector_b->y;
    product->y = vector_a->z * vector_b->x - vector_a->x * vector_b->z;
    product->z = vector_a->x * vector_b->y - vector_a->y * vector_b->x;
}

void Vector3D_Lerp(const Vector3D *start, const Vector3D *end, double fraction, Vector3D *result) {
    Vector3D_Subtract(end, start, result);
    Vector3D_Multiply(result, fraction);
    Vector3D_Add(start, result, result);
}

void Vector3D_Normalize(Vector3D *vector) {
    Vector3D_Divide(vector, Vector3D_Magnitude(vector));
}

void Vector3D_Direction(const Vector3D *start, const Vector3D *end, Vector3D *direction) {
    Vector3D_Subtract(end, start, direction);
    Vector3D_Normalize(direction);
}

double Vector3D_Dot(const Vector3D *vector_a, const Vector3D *vector_b) {
    return vector_a->x * vector_b->x + 
           vector_a->y * vector_b->y + 
           vector_a->z * vector_b->z;
}

double Vector3D_Magnitude(const Vector3D *vector) {
    return sqrt(Vector3D_Dot(vector, vector));
}

double Vector3D_Distance(const Vector3D *start, const Vector3D *end) {
    Vector3D distance;
    Vector3D_Subtract(end, start, &distance);

    return Vector3D_Magnitude(&distance);
}

Vector3D *Vector3D_Allocate(size_t num_elements) {
    Vector3D *ptr = (Vector3D *) malloc(sizeof(Vector3D) * num_elements);

    if (ptr == NULL) {
        fprintf(stderr, "Vector3D_Allocate: malloc returned NULL for %zu elements", num_elements);
    }

    return ptr;
}