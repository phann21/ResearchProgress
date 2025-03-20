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

//3x3 * 3x1 -> 3x1
double* matrix_multiply_3x3_3x1(double matrix_1[DIMENSION][DIMENSION], double matrix_2[DIMENSION]){
    double result[DIMENSION] = {0,0,0};
    for (int i=0; i<DIMENSION; i++){
        for (int j=0; j<DIMENSION; j++){
            result[i] += matrix_1[i][j] * matrix_2[j];
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

//1x1
double create_K(double x, double y, double z){
    double K = pow(x,(double)2) + pow(y,(double)2) + pow(z,(double)2);
}

//3x3
double* create_A(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4){
    double A[DIMENSION][DIMENSION] = {{x2-x1, x3-x1, x4-x1},
                                    {y2-y1, y3-y1, y4-y1},
                                    {z2-z1, z3-z1, z4-z1}};
    return A;
}

//3x1
double* create_b(double K1, double K2, double K3, double K4){
    double b[DIMENSION] = {(K2 - K1), (K3 - K1), (K4 - K1)};
    return b;
}

double create_ri_1(double ti, double t1){
    return ti - t1;
}

//3x1
double* create_d(double r2_1, double r3_1, double r4_1){
    double d[DIMENSION] = {r2_1, r3_1, r4_1};
    return d;
}

//3x1
double* create_e(double r2_1, double r3_1, double r4_1){
    double e[DIMENSION] = {pow(r2_1, (double)2), pow(r3_1, (double)2), pow(r4_1, (double)2)};
    return e;
}

// ar1^2 + br1 + c = 0
double* solve_w(double* A, double* b, double* e){
    double* inverted_A = matrix_invert(A);
    double* diff_b_e = matrix_subtraction(b,e);
    double* half_diff_b_e = matrix_scalar_multiply_3x1(0.5, diff_b_e);
    double* result = matrix_multiply_3x3_3x1(inverted_A, half_diff_b_e);
    return result;
}

double* solve_f(double* A, double* d){
    double* inverted_A = matrix_invert(A);
    double* result = matrix_multiply_3x3_3x1(inverted_A, d);
    return result;
}

double solve_alpha(double* f){
    double result = matrix_multiply_1x3_3x1(f, f);
    return result - (double)1;
}

double solve_beta(double* f, double* w){
    double result = matrix_multiply_1x3_3x1(f, w);
    return (double)-2 * result;
}

double solve_sigma(double* w){
    double result = matrix_multiply_1x3_3x1(w, w);
    return result;
}

//quadratic formula
double solve_r1(double alpha, double beta, double sigma){
    double result_1 = (-1 * beta + pow(pow(beta,2) - 4 * alpha * sigma, 0.5)) / (2 * alpha);
    double result_2 = (-1 * beta - pow(pow(beta,2) - 4 * alpha * sigma, 0.5)) / (2 * alpha);

    if (result_1 < 0){
        return result_2;
    }
    else{
        return result_1;
    }
}

// -> 3x1
double* solve_points(double* A, double* b, double* d, double* e, double r1){
    double* invert_A = matrix_invert(A);
    double* r1_d = matrix_scalar_multiply_3x1(r1, d);
    double* diff_b_e = matrix_subtraction(e,b);
    double* half_b_e = matrix_scalar_multiply_3x1(0.5, diff_b_e);
    double* add_b_d_e = matrix_addition(half_b_e, r1_d);
    double* mult_a_b_d_e = matrix_multiply_3x3_3x1(invert_A, add_b_d_e);
    double* result = matrix_scalar_multiply_3x1(-1, mult_a_b_d_e);
    return result;
}

int main(){
    //inputs
    double t1 = 0;
    double t2 = 0;
    double t3 = 0;
    double t4 = 0;

    double x1 = 0;
    double x2 = 0;
    double x3 = 0;
    double x4 = 0;

    double y1 = 0;
    double y2 = 0;
    double y3 = 0;
    double y4 = 0;

    double z1 = 0;
    double z2 = 0;
    double z3 = 0;
    double z4 = 0;

    double r2_1 = create_ri_1(t2, t1);
    double r3_1 = create_ri_1(t3, t1);
    double r4_1 = create_ri_1(t4, t1);

    double K1 = create_K(x1, y1, z1);
    double K2 = create_K(x2, y2, z2);
    double K3 = create_K(x3, y3, z3);
    double K4 = create_K(x4, y4, z4);

    double* A = create_A(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);
    double* b = create_b(K1, K2, K3, K4);
    double* d = create_d(r2_1, r3_1, r4_1);
    double* e = create_e(r2_1, r3_1, r4_1);

    double* w = solve_w(A, b, e);
    double* f = solve_f(A, d);

    double alpha = solve_alpha(f);
    double beta = solve_beta(f, w);
    double sigma = solve_sigma(w);

    double r1 = solve_r1(alpha, beta, sigma);

    double* result = solve_points(A, b, d, e, r1);
}