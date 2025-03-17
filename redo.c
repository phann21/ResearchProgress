#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DIMENSION 3

//3x3 -> 1x1
double matrix_determinant_3x3(double input_matrix[][DIMENSION]){
    double first_line = input_matrix[0][0] * (input_matrix[1][1] * input_matrix[2][2] - input_matrix[1][2] * input_matrix[2][1]); 
    double second_line = input_matrix[0][1] * (input_matrix[1][0] * input_matrix[2][2] - input_matrix[1][2] * input_matrix[2][0]);
    double third_line = input_matrix[0][3] * (input_matrix[1][0] * input_matrix[2][1] - input_matrix[1][1] * input_matrix[2][0]);

    return first_line - second_line + third_line;
}

//3x3 -> 3x3
double* matrix_invert(double input_matrix[][DIMENSION]){
    double output_matrix[DIMENSION][DIMENSION] = {  {input_matrix[1][1] * input_matrix[2][2] - input_matrix[1][2] * input_matrix[2][1], (-1) * (input_matrix[1][0] * input_matrix[2][2] - input_matrix[1][2] * input_matrix[2][0]), input_matrix[1][0] * input_matrix[2][1] - input_matrix[1][1] * input_matrix[2][0]},
                                                    {(-1) * (input_matrix[0][1] * input_matrix[2][2] - input_matrix[0][2] * input_matrix[2][1]), (input_matrix[0][0] * input_matrix[2][2] - input_matrix[0][2] * input_matrix[2][0]), (-1) * (input_matrix[0][0] * input_matrix[2][1] - input_matrix[0][1] * input_matrix[2][0])},
                                                    {input_matrix[0][1] * input_matrix[1][2] - input_matrix[0][2] * input_matrix[1][1], (-1) * (input_matrix[0][0] * input_matrix[1][2] - input_matrix[0][2] * input_matrix[1][0]), input_matrix[0][0] * input_matrix[1][1] - input_matrix[0][1] * input_matrix[1][0]}};
    
    double determinant = matrix_determinant_3x3(input_matrix);
    for (int i = 0; i < DIMENSION; i++){
        for (int j = 0; j < DIMENSION; j++){
            output_matrix[i][j] = output_matrix[i][j] / determinant;
        }
    }
    return output_matrix;
}

//3x1 * 1x3 -> 3x3
double* matrix_multiply_3x1_1x3(double matrix_1[DIMENSION], double matrix_2[DIMENSION]){
    double result[DIMENSION][DIMENSION] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i=0; i<DIMENSION; i++){
        for (int j=0; j<DIMENSION; j++){
            result[i][j] = matrix_1[i] * matrix_2[j];
        }
    }
    return result;
}

//1x3 * 3x1 -> 1x1
double matrix_multiply_1x3_3x1(double matrix_1[DIMENSION], double matrix_2[DIMENSION]){
    double sum = 0;
    for (int i=0; i<DIMENSION; i++){
        sum += matrix_1[i] * matrix_2[i];
    }
    return sum;
}

//3x1 * 1x1 -> 3x1
double* matrix_scalar_multiply_3x1(double scalar, double matrix_1[DIMENSION]){
    double result[DIMENSION] = {0,0,0};
    for (int i=0; i<DIMENSION; i++){
        result[i] = matrix_1[i] * scalar;
    }
    return result;
}

//3x1 - 3x1 -> 3x1
double* matrix_subtraction(double matrix_1[DIMENSION], double matrix_2[DIMENSION]){
    double result[DIMENSION] = {0,0,0};
    for (int i=0; i<DIMENSION; i++){
        result[i] = matrix_1[i] - matrix_2[i];
    }
    return result;
}

//3x1 + 3x1 -> 3x1
double* matrix_addition(double matrix_1[DIMENSION], double matrix_2[DIMENSION]){
    double result[DIMENSION] = {0,0,0};
    for (int i=0; i<DIMENSION; i++){
        result[i] = matrix_1[i] + matrix_2[i];
    }
    return result;
}

double solve_K(double x, double y, double z){
    double K = pow(x,(double)2) + pow(y,(double)2) + pow(z,(double)2);
}

//3x3
double* solve_A(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4){
    double A[DIMENSION][DIMENSION] = {{x2-x1, x3-x1, x4-x1},
                                    {y2-y1, y3-y1, y4-y1},
                                    {z2-z1, z3-z1, z4-z1}};
    return A;
}

//3x1
double* solve_b(double K1, double K2, double K3, double K4){
    double b[DIMENSION] = {0.5*(K2 - K1), 0.5*(K3 - K1), 0.5*(K4 - K1)};
    return b;
}

double create_ri_1(double ti, double t1){
    return ti - t1;
}

//3x1
double* solve_d(double r2_1, double r3_1, double r4_1){
    double d[DIMENSION] = {r2_1, r3_1, r4_1};
    return d;
}

//3x1
double* solve_e(double r2_1, double r3_1, double r4_1){
    double e[DIMENSION] = {pow(r2_1, (double)2), pow(r3_1, (double)2), pow(r4_1, (double)2)};
    return e;
}

// ar1^2 + br1 + c = 0
 
int main() {
    
}