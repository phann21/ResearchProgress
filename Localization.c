/* ************************************************************************************************************* */
// *****BEGIN INCLUDES AND DEFINES***************************************************************************** //
/* *********************************************************************************************************** */

#include <stdio.h>
#include <math.h>
#include <time.h>

#define DIMENSION 3
#define C_speed 1481*1000 //speed of sound in mm/s
/* *********************************************************************************************************** */
// *****END INCLUDES AND DEFINES***************************************************************************** //
/* ********************************************************************************************************* */

/* ******************************************************************************************************* */
// *****BEGIN MATRIX METHODS***************************************************************************** //
/* ***************************************************************************************************** */


// 3x3 -> 1x1
double matrix_determinant_3x3(double (*input_matrix)[DIMENSION]){
    double first_line = input_matrix[0][0] * ((input_matrix[1][1] * input_matrix[2][2]) - (input_matrix[1][2] * input_matrix[2][1])); 
    double second_line = input_matrix[0][1] * ((input_matrix[1][0] * input_matrix[2][2]) - (input_matrix[1][2] * input_matrix[2][0]));
    double third_line = input_matrix[0][2] * ((input_matrix[1][0] * input_matrix[2][1]) - (input_matrix[1][1] * input_matrix[2][0]));

    return first_line - second_line + third_line;
}

// 3x3 -> 3x3
void matrix_transpose_3x3(double (*input_matrix)[DIMENSION], double (*output_matrix)[DIMENSION]){
    for (int i=0; i<DIMENSION; i++){
        for (int j=0; j<DIMENSION; j++){
            output_matrix[i][j] = input_matrix[j][i];
        }
    }
}

// 3x3 -> 3x3
int matrix_invert(double (*input_matrix)[DIMENSION], double (*output_matrix)[DIMENSION]){
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
    if (determinant == 0){
        return 1;
    }
    for (int i = 0; i < DIMENSION; i++){
        for (int j = 0; j < DIMENSION; j++){
            adjoint_matrix[i][j] = adjoint_matrix[i][j] / determinant;
        }
    }
    // printf("determinant %lf\n", determinant);
    matrix_transpose_3x3(adjoint_matrix, output_matrix);
    return 0;
    // for (int k=0; k<DIMENSION; k++){
    //     for (int l=0; l<DIMENSION; l++){
    //         printf("Inverted A %lf ", output_matrix[k][l]);
    //     }
    //     printf("\n");
    // }
}

// 3x3 * 3x1 -> 3x1
void matrix_multiply_3x3_3x1(double (*matrix_1)[DIMENSION], double* matrix_2, double* result){
    result[0] = 0;
    result[1] = 0;
    result[2] = 0;
    for (int i=0; i<DIMENSION; i++){
        for (int j=0; j<DIMENSION; j++){
            result[i] += matrix_1[i][j] * matrix_2[j];
        }
    }
}

// 1x3 * 3x1 -> 1x1
double matrix_multiply_1x3_3x1(double* matrix_1, double* matrix_2){
    double sum = 0;
    for (int i=0; i<DIMENSION; i++){
        sum += matrix_1[i] * matrix_2[i];
    }
    return sum;
}

// 3x1 * 1x1 -> 3x1
void matrix_scalar_multiply_3x1(double scalar, double* matrix_1, double* result){
    for (int i=0; i<DIMENSION; i++){
        result[i] = matrix_1[i] * scalar;
    }
}

// 3x1 - 3x1 -> 3x1
void matrix_subtraction(double* matrix_1, double* matrix_2, double* result){
    for (int i=0; i<DIMENSION; i++){
        result[i] = matrix_1[i] - matrix_2[i];
    }
}
/* ***************************************************************************************************** */
// *****END MATRIX METHODS***************************************************************************** //
/* *************************************************************************************************** */

/* *************************************************************************************************************** */
// *****BEGIN INITIALIZING VARIABLES***************************************************************************** //
/* ************************************************************************************************************* */

// 1x1 ri_1 = C_speed(t_i - t1)
double create_ri_1(double t_i, double t1){
    return C_speed*(t_i - t1);
}

// 1x1 K_i = x_i^2 + y_i^2 + z_i^2
double create_Ki(double x, double y, double z){
    double K = pow(x, 2) + pow(y, 2) + pow(z, 2);
    return K;
}

// 3x3
// [[x2-x1; y2-y1; z2-z1;],
// [x3-x1; y3-y1; z3-z1;],
// [x4-x1; y4-y1; z4-z1;]]
void create_A(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4, double (*A)[DIMENSION]){
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

// 3x1
void create_b(double r2_1, double r3_1, double r4_1, double* b){
    b[0] = r2_1;
    b[1] = r3_1;
    b[2] = r4_1;
    // for (int k=0; k<DIMENSION; k++){
    //     printf("B %lf\n", b[k]);
    // }
}

// 3x1
void create_d(double r2_1, double r3_1, double r4_1, double K1, double K2, double K3, double K4, double* d){
    d[0] = pow(r2_1, 2) - K2 + K1;
    d[1] = pow(r3_1, 2) - K3 + K1;
    d[2] = pow(r4_1, 2) - K4 + K1;
    // for (int k=0; k<DIMENSION; k++){
    //     printf("D %lf\n", d[k]);
    // }
}

/* ************************************************************************************************************* */
// *****END INITIALIZING VARIABLES***************************************************************************** //
/* *********************************************************************************************************** */

/* *************************************************************************************************************** */
// *****BEGIN INTERMEDIATE VARIABLES***************************************************************************** //
/* ************************************************************************************************************* */

// f = -(A^-1) * b
void solve_f(double (*invert_A)[DIMENSION], double* b, double* result){
    // for (int k=0; k<DIMENSION; k++){
    //     for (int l=0; l<DIMENSION; l++){
    //         printf("inverted_A %lf ", inverted_A[k][l]);
    //     }
    //     printf("\n");
    // }

    // for (int k=0; k<DIMENSION; k++){
    //     printf("b %lf\n", b[k]);
    // }
    double A_invert_b[DIMENSION];
    matrix_multiply_3x3_3x1(invert_A, b, A_invert_b);
    // for (int k=0; k<DIMENSION; k++){
    //     printf("A_invert_b %lf\n", A_invert_b[k]);
    // }
    matrix_scalar_multiply_3x1(-1, A_invert_b, result);
    // for (int k=0; k<DIMENSION; k++){
    //     printf("f %lf\n", result[k]);
    // }
}

// g = 0.5(-(A^-1) * d)
void solve_g(double (*invert_A)[DIMENSION], double* d, double* result){
    double A_invert_d[DIMENSION];
    matrix_multiply_3x3_3x1(invert_A, d, A_invert_d);
    matrix_scalar_multiply_3x1(0.5, A_invert_d, result);
    // for (int k=0; k<DIMENSION; k++){
    //     printf("g %lf\n", result[k]);
    // }
}

// alpha = fTf - 1
double solve_alpha(double* f){
    double result = matrix_multiply_1x3_3x1(f, f);
    return result - 1;
}

// beta = -2 * gTf
double solve_beta(double* f, double* g){
    double result = matrix_multiply_1x3_3x1(f, g);
    return -2 * result;
}

// gamma = gTg
double solve_gamma(double* g){
    double result = matrix_multiply_1x3_3x1(g, g);
    return result;
}

/* ************************************************************************************************************* */
// *****END INTERMEDIATE VARIABLES***************************************************************************** // 
/* *********************************************************************************************************** */

/* ************************************************************************************************ */
// *****BEGIN SOLVERS***************************************************************************** //
/* ********************************************************************************************** */

// quadratic formula -beta +/- sqrt(b)
double solve_r1(double alpha, double beta, double gamma){
    double result_1 = (-1 * beta + sqrt(pow(beta, 2) - 4 * alpha * gamma)) / (2 * alpha);
    double result_2 = (-1 * beta - sqrt(pow(beta, 2) - 4 * alpha * gamma)) / (2 * alpha);
    if (result_1 > 0 && result_2 > 0){
        if (result_1 > result_2){
            printf("result_2 chosen r2: %lf\n", result_2);
            printf("result_2 chosen r1: %lf\n", result_1);
            return result_2;
        }
        else{
            printf("result_1 chosen r1: %lf\n", result_1);
            printf("result_1 chosen r2: %lf\n", result_2);
            return result_1;
        }
    }
    if (result_1 < 0){
        printf("result_2 chosen r2: %lf\n", result_2);
        printf("result_2 chosen r1: %lf\n", result_1);
        return result_2;
    }
    else{
        printf("result_1 chosen r1: %lf\n", result_1);
        printf("result_1 chosen r2: %lf\n", result_2);
        return result_1;
    }
}

//returns -> 3x1
void solve_points(double (*invert_A)[DIMENSION], double* b, double* d, double r1, double* result){
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

/* **************************************************************************************************** */
// *****END SOLVERS*********************************************************************************** //
/* ************************************************************************************************** */

int main(){
    clock_t begin = clock();
    /******BEGIN USER INPUTS******/
    // inputs in seconds
    // double t1 = 1.350438892640108e-05;
    // double t2 = 9.549044985638725e-06;
    // double t3 = 9.548939945566158e-06;
    // double t4 = 9.548939945566158e-06;
    // double t1 = 5.4855087134611484e-06;
    // double t2 = 1.3037280047565443e-05;
    // double t3 = 1.8108311020723282e-05;
    // double t4 = 6.2617275459120215e-06;
    double t1 = 1.1012495901620587e-05;
    double t2 = 1.3037280047565443e-05;
    double t3 = 1.8108311020723282e-05;
    double t4 = 1.181151632986894e-05;
    

    

    // double x1 = 0;
    // double y1 = 0;
    // double z1 = 0;

    // double x2 = 0;
    // double y2 = 10;
    // double z2 = -10;

    // double x3 = -8.66;
    // double y3 = -5;
    // double z3 = -10;

    // double x4 = 8.66;
    // double y4 = -5;
    // double z4 = -10;
    double x1 = 0;
    double y1 = 0;
    double z1 = 0;

    double x2 = -10;
    double y2 = 17.321;
    double z2 = -10;

    double x3 = -10;
    double y3 = -17.321;
    double z3 = -10;

    double x4 = 10;
    double y4 = 0;
    double z4 = -10;
    // *****END USER INPUTS****** //

    double r2_1 = create_ri_1(t2, t1);
    double r3_1 = create_ri_1(t3, t1);
    double r4_1 = create_ri_1(t4, t1);

    double K1 = create_Ki(x1, y1, z1);
    double K2 = create_Ki(x2, y2, z2);
    double K3 = create_Ki(x3, y3, z3);
    double K4 = create_Ki(x4, y4, z4);
    double A[DIMENSION][DIMENSION];
    create_A(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, A);
    double invert_A[DIMENSION][DIMENSION];
    if (matrix_invert(A, invert_A) == 1){
        return 1;
    }
    double b[DIMENSION];
    create_b(r2_1, r3_1, r4_1, b);
    double d[DIMENSION];
    create_d(r2_1, r3_1, r4_1, K1, K2, K3, K4, d);

    double f[DIMENSION];
    solve_f(invert_A, b, f);
    double g[DIMENSION];
    solve_g(invert_A, d, g);

    double alpha = solve_alpha(f);
    double beta = solve_beta(f, g);
    double gamma = solve_gamma(g);

    double r1 = solve_r1(alpha, beta, gamma);

    double result[DIMENSION];
    solve_points(invert_A, b, d, r1, result);
    /* here, do your time-consuming job */

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    for (int i=0; i<DIMENSION; i++){
        printf("Value %lf\n", result[i]);
    }
    printf("Time spent: %lfs\n", time_spent);
}