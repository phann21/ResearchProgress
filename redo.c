#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DIMENSION 3
#define C_speed 1481*1000 //speed of sound in mm/s

void matrix_transpose_3x3(double input_matrix[][DIMENSION], double output_matrix[][DIMENSION]){
    for (int i=0; i<DIMENSION; i++){
        for (int j=0; j<DIMENSION; j++){
            output_matrix[i][j] = input_matrix[j][i];
        }
    }
}

//3x3 -> 1x1
double matrix_determinant_3x3(double input_matrix[][DIMENSION]){
    double first_line = input_matrix[0][0] * ((input_matrix[1][1] * input_matrix[2][2]) - (input_matrix[1][2] * input_matrix[2][1])); 
    double second_line = input_matrix[0][1] * ((input_matrix[1][0] * input_matrix[2][2]) - (input_matrix[1][2] * input_matrix[2][0]));
    double third_line = input_matrix[0][2] * ((input_matrix[1][0] * input_matrix[2][1]) - (input_matrix[1][1] * input_matrix[2][0]));

    return first_line - second_line + third_line;
}

//3x3 -> 3x3
void matrix_invert(double input_matrix[][DIMENSION], double output_matrix[][DIMENSION]){
    double adjoint_matrix[DIMENSION][DIMENSION];
    adjoint_matrix[0][0] = input_matrix[1][1] * input_matrix[2][2] - input_matrix[1][2] * input_matrix[2][1];
    adjoint_matrix[0][1] = (-1) * (input_matrix[1][0] * input_matrix[2][2] - input_matrix[1][2] * input_matrix[2][0]);
    adjoint_matrix[0][2] = input_matrix[1][0] * input_matrix[2][1] - input_matrix[1][1] * input_matrix[2][0];
    adjoint_matrix[1][0] = (-1) * (input_matrix[0][1] * input_matrix[2][2] - input_matrix[0][2] * input_matrix[2][1]);
    adjoint_matrix[1][1] = (input_matrix[0][0] * input_matrix[2][2] - input_matrix[0][2] * input_matrix[2][0]);
    adjoint_matrix[1][2] = (-1) * (input_matrix[0][0] * input_matrix[2][1] - input_matrix[0][1] * input_matrix[2][0]);
    adjoint_matrix[2][0] = input_matrix[0][1] * input_matrix[1][2] - input_matrix[0][2] * input_matrix[1][1];
    adjoint_matrix[2][1] = (-1) * (input_matrix[0][0] * input_matrix[1][2] - input_matrix[0][2] * input_matrix[1][0]);
    adjoint_matrix[2][2] = input_matrix[0][0] * input_matrix[1][1] - input_matrix[0][1] * input_matrix[1][0];
    
    double determinant = matrix_determinant_3x3(input_matrix);
    for (int i = 0; i < DIMENSION; i++){
        for (int j = 0; j < DIMENSION; j++){
            adjoint_matrix[i][j] = adjoint_matrix[i][j] / determinant;
        }
    }
    printf("determinant %lf\n", determinant);
    matrix_transpose_3x3(adjoint_matrix, output_matrix);
    for (int k=0; k<DIMENSION; k++){
        for (int l=0; l<DIMENSION; l++){
            printf("Inverted A %lf ", output_matrix[k][l]);
        }
        printf("\n");
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

double create_ri_1(double ti, double t1){
    return C_speed*(ti - t1);
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
    // for (int k=0; k<DIMENSION; k++){
    //     for (int l=0; l<DIMENSION; l++){
    //         printf("A %lf ", A[k][l]);
    //     }
    //     printf("\n");
    // }
}

//3x1
void create_b(double r2_1, double r3_1, double r4_1, double b[DIMENSION]){
    b[0] = r2_1;
    b[1] = r3_1;
    b[2] = r4_1;
    // for (int k=0; k<DIMENSION; k++){
    //     printf("B %lf\n", b[k]);
    // }
}

//3x1
void create_d(double r2_1, double r3_1, double r4_1, double K1, double K2, double K3, double K4, double d[DIMENSION]){
    d[0] = pow(r2_1, 2) - K2 + K1;
    d[1] = pow(r3_1, 2) - K3 + K1;
    d[2] = pow(r3_1, 2) - K4 + K1;
    // for (int k=0; k<DIMENSION; k++){
    //     printf("D %lf\n", d[k]);
    // }
}

// -(A)^-1 * b
void solve_w(double A[][DIMENSION], double b[DIMENSION], double result[DIMENSION]){
    double inverted_A[DIMENSION][DIMENSION];
    matrix_invert(A, inverted_A);
    double A_invert_b[DIMENSION];
    matrix_multiply_3x3_3x1(inverted_A, b, A_invert_b);
    matrix_scalar_multiply_3x1(-1, A_invert_b, result);
}

void solve_f(double A[][DIMENSION], double d[DIMENSION], double result[DIMENSION]){
    double inverted_A[DIMENSION][DIMENSION];
    matrix_invert(A, inverted_A);
    double A_invert_d[DIMENSION];
    matrix_multiply_3x3_3x1(inverted_A, d, A_invert_d);
    matrix_scalar_multiply_3x1(0.5, A_invert_d, result);
}

double solve_alpha(double* w){
    double result = matrix_multiply_1x3_3x1(w, w);
    return result - 1;
}

double solve_beta(double* w, double* f){
    double result = matrix_multiply_1x3_3x1(f, w);
    return -2 * result;
}

double solve_gamma(double* f){
    double result = matrix_multiply_1x3_3x1(f, f);
    return result;
}

//quadratic formula
double solve_r1(double alpha, double beta, double gamma){
    double result_1 = (-1 * beta + pow(pow(beta,2) - 4 * alpha * gamma, 0.5)) / (2 * alpha);
    double result_2 = (-1 * beta - pow(pow(beta,2) - 4 * alpha * gamma, 0.5)) / (2 * alpha);
    if (result_1 < 0){
        printf("result 2 chosen r1_2: %lf\n", result_2);
        printf("result 2 chosen r1_1: %lf\n", result_1);
        return result_2;
    }
    else{
        printf("result 1 chosen r1_1: %lf\n", result_1);
        printf("result 1 chosen r1_2: %lf\n", result_2);
        return result_1;
    }
}

// -> 3x1
void solve_points(double A[][DIMENSION], double b[DIMENSION], double d[DIMENSION], double r1, double result[DIMENSION]){
    double invert_A[DIMENSION][DIMENSION];
    matrix_invert(A, invert_A);
    double r1_b[DIMENSION];
    matrix_scalar_multiply_3x1(r1, b, r1_b);
    double A_invert_b_r1[DIMENSION];
    matrix_multiply_3x3_3x1(invert_A, r1_b, A_invert_b_r1);
    double negate_A_invert_b_r1[DIMENSION];
    matrix_scalar_multiply_3x1(-1, A_invert_b_r1, negate_A_invert_b_r1);
    double A_invert_d[DIMENSION];
    matrix_multiply_3x3_3x1(invert_A, d, A_invert_d);
    double half_A_invert_d[DIMENSION];
    matrix_scalar_multiply_3x1(0.5, A_invert_d, half_A_invert_d);
    matrix_subtraction(negate_A_invert_b_r1, half_A_invert_d, result);
}

int main(){
    //User inputs BEGIN//
    //inputs seconds
    double t1 = 1.350438892640108e-05;
    double t2 = 9.548939945566158e-06;
    double t3 = 9.548939945566158e-06;
    double t4 = 9.548939945566158e-06;

    double x1 = 0;
    double y1 = 0;
    double z1 = 0;

    double x2 = 0;
    double y2 = 10;
    double z2 = -10;

    double x3 = -8.66;
    double y3 = -5;
    double z3 = -10;

    double x4 = 8.66;
    double y4 = -5;
    double z4 = -10;
    //User inputs END//

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
    create_b(r2_1, r3_1, r4_1, b);
    double d[DIMENSION];
    create_d(r2_1, r3_1, r4_1, K1, K2, K3, K4, d);

    double w[DIMENSION];
    solve_w(A, b, w);
    double f[DIMENSION];
    solve_f(A, d, f);

    double alpha = solve_alpha(w);
    double beta = solve_beta(w, f);
    double gamma = solve_gamma(f);

    double r1 = solve_r1(alpha, beta, gamma);

    double result[DIMENSION];
    solve_points(A, b, d, r1, result);

    // for (int i=0; i<DIMENSION; i++){
    //     printf("Value %lf\n", result[i]);
    // }
}