#include "vector.h"

void apply_rotation(double rotation_matrix[3][3], Vector *point) {
     point->x = rotation_matrix[0][0] * (point->x) + 
                rotation_matrix[0][1] * (point->y) + 
                rotation_matrix[0][2] * (point->z);

     point->y = rotation_matrix[1][0] * (point->x) + 
                rotation_matrix[1][1] * (point->y) + 
                rotation_matrix[1][2] * (point->z);

     point->z = rotation_matrix[2][0] * (point->x) + 
                rotation_matrix[2][1] * (point->y) + 
                rotation_matrix[2][2] * (point->z);
}