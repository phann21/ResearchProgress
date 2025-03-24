#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DIMENSION 3
#define C_speed 1481000 //speed of sound in mm/s

//3x3 -> 1x1
double matrix_determinant_3x3(double input_matrix[][DIMENSION]){
    double first_line = input_matrix[0][0] * (input_matrix[1][1] * input_matrix[2][2] - input_matrix[1][2] * input_matrix[2][1]); 
    double second_line = input_matrix[0][1] * (input_matrix[1][0] * input_matrix[2][2] - input_matrix[1][2] * input_matrix[2][0]);
    double third_line = input_matrix[0][3] * (input_matrix[1][0] * input_matrix[2][1] - input_matrix[1][1] * input_matrix[2][0]);

    return first_line - second_line + third_line;
}

//3x3 -> 3x3
void matrix_invert(double input_matrix[][DIMENSION], double output_matrix[][DIMENSION]){
    output_matrix[0][0] = input_matrix[1][1] * input_matrix[2][2] - input_matrix[1][2] * input_matrix[2][1];
    output_matrix[0][1] = (-1) * (input_matrix[1][0] * input_matrix[2][2] - input_matrix[1][2] * input_matrix[2][0]);
    output_matrix[0][2] = input_matrix[1][0] * input_matrix[2][1] - input_matrix[1][1] * input_matrix[2][0];
    output_matrix[1][0] = (-1) * (input_matrix[0][1] * input_matrix[2][2] - input_matrix[0][2] * input_matrix[2][1]);
    output_matrix[1][1] = (input_matrix[0][0] * input_matrix[2][2] - input_matrix[0][2] * input_matrix[2][0]);
    output_matrix[1][2] = (-1) * (input_matrix[0][0] * input_matrix[2][1] - input_matrix[0][1] * input_matrix[2][0]);
    output_matrix[2][0] = input_matrix[0][1] * input_matrix[1][2] - input_matrix[0][2] * input_matrix[1][1];
    output_matrix[2][1] = (-1) * (input_matrix[0][0] * input_matrix[1][2] - input_matrix[0][2] * input_matrix[1][0]);
    output_matrix[2][2] = input_matrix[0][0] * input_matrix[1][1] - input_matrix[0][1] * input_matrix[1][0];
    
    double determinant = matrix_determinant_3x3(input_matrix);
    for (int i = 0; i < DIMENSION; i++){
        for (int j = 0; j < DIMENSION; j++){
            output_matrix[i][j] = output_matrix[i][j] / determinant;
        }
    }
}

//3x3 * 3x1 -> 3x1
void matrix_multiply_3x3_3x1(double matrix_1[DIMENSION][DIMENSION], double matrix_2[DIMENSION], double result[DIMENSION]){
    for (int i=0; i<DIMENSION; i++){
        for (int j=0; j<DIMENSION; j++){
            result[i] += matrix_1[i][j] * matrix_2[j];
        }
    }
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
void matrix_scalar_multiply_3x1(double scalar, double matrix_1[DIMENSION], double result[DIMENSION]){
    for (int i=0; i<DIMENSION; i++){
        result[i] = matrix_1[i] * scalar;
    }
}

//3x1 - 3x1 -> 3x1
void matrix_subtraction(double matrix_1[DIMENSION], double matrix_2[DIMENSION], double result[DIMENSION]){
    for (int i=0; i<DIMENSION; i++){
        result[i] = matrix_1[i] - matrix_2[i];
    }
}

//3x1 + 3x1 -> 3x1
void matrix_addition(double matrix_1[DIMENSION], double matrix_2[DIMENSION], double result[DIMENSION]){
    for (int i=0; i<DIMENSION; i++){
        result[i] = matrix_1[i] + matrix_2[i];
    }
}

//1x1
double create_K(double x, double y, double z){
    double K = pow(x, 2) + pow(y, 2) + pow(z, 2);
    return K;
}

//3x3
// [[x2-x1; y2-y1; z2-z1;],
// [x3-x1; y3-y1; z3-z1;],
// [x4-x1; y4-y1; z4-z1;]]
void create_A(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4, double A[][DIMENSION]){
    A[0][0] = x2-x1;
    A[0][1] = y2-y1;
    A[0][2] = z2-z1;

    A[1][0] = x3-x1;
    A[1][1] = y3-y1;
    A[1][2] = z3-z1;

    A[2][0] = x4-x1;
    A[2][1] = y4-y1;
    A[2][2] = z4-z1;
}

//3x1
void create_b(double K1, double K2, double K3, double K4, double b[DIMENSION]){
    b[0] = K2 - K1;
    b[1] = K3 - K1;
    b[2] = K4 - K1;
}

double create_ri_1(double ti, double t1){
    return C_speed*(ti - t1);
}

//3x1
void create_d(double r2_1, double r3_1, double r4_1, double d[DIMENSION]){
    d[0] = r2_1;
    d[1] = r3_1;
    d[2] = r4_1;
}

//3x1
void create_e(double r2_1, double r3_1, double r4_1, double e[DIMENSION]){
    e[0] = pow(r2_1, 2);
    e[1] = pow(r3_1, 2);
    e[2] = pow(r4_1, 2);
}

// ar1^2 + br1 + c = 0
void solve_w(double A[][DIMENSION], double b[DIMENSION], double e[DIMENSION], double result[DIMENSION]){
    double inverted_A[DIMENSION][DIMENSION];
    matrix_invert(A, inverted_A);
    double diff_b_e[DIMENSION];
    matrix_subtraction(b,e, diff_b_e);
    double half_diff_b_e[DIMENSION];
    matrix_scalar_multiply_3x1(0.5, diff_b_e, half_diff_b_e);
    matrix_multiply_3x3_3x1(inverted_A, half_diff_b_e, result);
}

void solve_f(double A[][DIMENSION], double d[DIMENSION], double result[DIMENSION]){
    double inverted_A[DIMENSION][DIMENSION];
    matrix_invert(A, inverted_A);
    matrix_multiply_3x3_3x1(inverted_A, d, result);
}

double solve_alpha(double* f){
    double result = matrix_multiply_1x3_3x1(f, f);
    return result - 1;
}

double solve_beta(double* f, double* w){
    double result = matrix_multiply_1x3_3x1(f, w);
    return -2 * result;
}

double solve_gamma(double* w){
    double result = matrix_multiply_1x3_3x1(w, w);
    return result;
}

//quadratic formula
double solve_r1(double alpha, double beta, double gamma){
    double result_1 = (-1 * beta + pow(pow(beta,2) - 4 * alpha * gamma, 0.5)) / (2 * alpha);
    double result_2 = (-1 * beta - pow(pow(beta,2) - 4 * alpha * gamma, 0.5)) / (2 * alpha);
    if (result_1 < 0){
        return result_2;
    }
    else{
        return result_1;
    }
}

// -> 3x1
void solve_points(double A[][DIMENSION], double b[DIMENSION], double d[DIMENSION], double e[DIMENSION], double r1, double result[DIMENSION]){
    double invert_A[DIMENSION][DIMENSION];
    matrix_invert(A, invert_A);
    double r1_d[DIMENSION];
    matrix_scalar_multiply_3x1(r1, d, r1_d);
    double sum_e_b[DIMENSION];
    matrix_addition(e, b, sum_e_b);
    double half_b_e[DIMENSION]; 
    matrix_scalar_multiply_3x1(0.5, sum_e_b, half_b_e);
    double diff_b_e_d[DIMENSION];
    double negate_r1_d[DIMENSION];
    matrix_scalar_multiply_3x1(-1, r1_d, negate_r1_d);
    matrix_subtraction(half_b_e, negate_r1_d, diff_b_e_d);
    matrix_multiply_3x3_3x1(invert_A, diff_b_e_d, result);
}

int main(){
    //inputs us
    double t1 = 1.350438892640108e-05;
    double t2 = 1.9098293376788757e-05;
    double t3 = 1.9098293376788757e-05;
    double t4 = 1.909808997127745e-05;

    double x1 = 0;
    double y1 = 0;
    double z1 = 2;

    double x2 = -10;
    double y2 = 17.321;
    double z2 = 0;

    double x3 = -10;
    double y3 = -17.321;
    double z3 = 0;

    double x4 = 20;
    double y4 = 0;
    double z4 = 0;

    double r2_1 = create_ri_1(t2, t1);
    double r3_1 = create_ri_1(t3, t1);
    double r4_1 = create_ri_1(t4, t1);

    double K1 = create_K(x1, y1, z1);
    double K2 = create_K(x2, y2, z2);
    double K3 = create_K(x3, y3, z3);
    double K4 = create_K(x4, y4, z4);

    double A[DIMENSION][DIMENSION];
    create_A(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, A);
    double b[DIMENSION];
    create_b(K1, K2, K3, K4, b);
    double d[DIMENSION];
    create_d(r2_1, r3_1, r4_1, d);
    double e[DIMENSION];
    create_e(r2_1, r3_1, r4_1, e);

    double w[DIMENSION];
    solve_w(A, b, e, w);
    double f[DIMENSION];
    solve_f(A, d, f);

    double alpha = solve_alpha(f);
    double beta = solve_beta(f, w);
    double gamma = solve_gamma(w);

    double r1 = solve_r1(alpha, beta, gamma);

    double result[DIMENSION];
    solve_points(A, b, d, e, r1, result);

    for (int i=0; i<DIMENSION; i++){
        printf("Value %f\n", result[i]);
    }
}