void apply_rotation(double rotation_matrix[3][3],
                  double *x, double *y, double *z) {
    double x_val = *x;
    double y_val = *y;
    double z_val = *z;

    *x = rotation_matrix[0][0] * (x_val) + 
         rotation_matrix[0][1] * (y_val) + 
         rotation_matrix[0][2] * (z_val);

    *y = rotation_matrix[1][0] * (x_val) + 
         rotation_matrix[1][1] * (y_val) + 
         rotation_matrix[1][2] * (z_val);

    *z = rotation_matrix[2][0] * (x_val) + 
         rotation_matrix[2][1] * (y_val) + 
         rotation_matrix[2][2] * (z_val);
}