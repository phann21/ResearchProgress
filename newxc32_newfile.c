// /* 
//  * Author: Philip Han
//  */

// //#include "CPUConfig.h"
// //#include <xc.h>
// //#include <sys/attribs.h>
// //#include "my_delays.h"


// //MINI BOARD BUILTIN LEDs
// #define DATA_LED_IO TRISEbits.TRISE0=0
// #define STAT_LED_IO TRISEbits.TRISE1=0
// #define DATA_LED(bit) LATEbits.LATE0=bit 
// #define STAT_LED(bit) LATEbits.LATE1=bit 
// #define DIMENSION 3
// #define C 1480 // m/s^2

// //Inputs
// double t_diff[4] = {0,0,0,0};
// double t1;
// double t2;
// double t3;
// double t4;
// double input_matrix[DIMENSION][DIMENSION] = {{0,0,0},{0,0,0},{0,0,0}};
// double r_matrix[DIMENSION] = {0,0,0};
// double r1;

// //Outputs
// double output_matrix[DIMENSION];

// double ri_1_calculation(double ti, double t1){
//     return (double)(C * (ti - t1));
// }

// double k_generation(double xi, double yi, double zi){
//     return xi*xi + yi*yi + zi*zi;
// }

// //add matrices
// double* matrix_add(double input1[DIMENSION], double input2[DIMENSION]){
//     double output[DIMENSION];
//     for (int i=0; i<DIMENSION; i++){
//         output[i] = input1[i] + input2[i];
//     }
//     return output;
// }

// double* matrix_scalar_multiply(double input[DIMENSION], double scalar){
//     double output[DIMENSION];
//     for (int i=0; i<DIMENSION; i++){
//         output[i] = input[i] * scalar;
//     }
//     return output;
// }

// double* matrix_multiply(double input1[][DIMENSION], double input2[DIMENSION]){
//     double output[DIMENSION];
//     for (int i=0; i<DIMENSION; i++){
//         output[i] = 0;
//         for (int j=0; j<DIMENSION; j++){
//             output[i] += input1[i][j] * input2[j];
//         }
//     }
//     return output;
// }

// double* matrix_transpose(double input[][DIMENSION]){
//     double output[DIMENSION][DIMENSION];
//     for (int i=0; i<DIMENSION; i++){
//         for (int j=0; j<DIMENSION; j++){
//             output[i][j] = input[j][i];
//         }
//     }
//     return output;
// }

// double* populate_matrix1(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4){
//     double x2_1 = x2 - x1;
//     double x3_1 = x3 - x1;
//     double x4_1 = x4 - x1;
//     double y2_1 = y2 - y1;
//     double y3_1 = y3 - y1;
//     double y4_1 = y4 - y1;
//     double z2_1 = z2 - z1;
//     double z3_1 = z3 - z1;
//     double z4_1 = z4 - z1;
//     double matrix1[DIMENSION][DIMENSION] = {{-1 * x2_1, -1 * x3_1, -1 * x4_1}, {-1 * y2_1, -1 * y3_1, -1 * y4_1}, {-1 * z2_1, -1 * z3_1, -1 * z4_1}};
//     return matrix1;
// }

// double* populate_matrix3(double r2_1, double r3_1, double r4_1, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4){
//     double K1 = k_generation(x1, y1, z1);
//     double K2 = k_generation(x2, y2, z2);
//     double K3 = k_generation(x3, y3, z3);
//     double K4 = k_generation(x4, y4, z4);
//     double matrix3[DIMENSION];
    
//     matrix3[0] = r2_1 * r2_1 - K2 + K1;
//     matrix3[1] = r3_1 * r3_1 - K3 + K1;
//     matrix3[2] = r4_1 * r4_1 - K4 + K1;
//     return matrix3;
// }

// double calculate_r1(double matrixA[][DIMENSION], double matrixB[DIMENSION], double matrixC[][DIMENSION]){
//     double a = 1;
//     double b = matrix_scalar_multiply(matrixB, -2);
// }

// ///NOT DONE
// void calculate(double matrix2[DIMENSION], double r1, double t1, double t2, double t3, double t4, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4){
//     double output[DIMENSION];
//     double matrix1[DIMENSION][DIMENSION] = populate_matrix1(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double z1, double z2, double z3, double z4);
    
//     double r2_1 = ri_1_calculation(t2, t1);
//     double r3_1 = ri_1_calculation(t3, t1);
//     double r4_1 = ri_1_calculation(t4, t1);
    
//     matrix2 = {r2_1, r3_1, r4_1};
//     matrix2 = matrix_scalar_multiply(matrix2, r1);
    
//     double matrix3[DIMENSION] = populate_matrix3(r2_1, r3_1, r4_1, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);
//     matrix3 = matrix_scalar_multiply(matrix3, 0.5);
    
//     double combined[DIMENSION];
//     combined = matrix_add(matrix2, matrix3);
    
//     double output = matrix_multiply(matrix1, combined);
//     return output;
// }

// //https://www.aidanmocke.com/blog/2018/04/03/config/
// void set_performance_mode()
// {   
    
//     unsigned int cp0;

//     // Unlock Sequence
//     asm volatile("di"); // Disable all interrupts
//     SYSKEY = 0x00000000;
//     SYSKEY = 0xAA996655;
//     SYSKEY = 0x556699AA;  

//     // PB1DIV
//     // Peripheral Bus 1 cannot be turned off, so there's no need to turn it on
//     PB1DIVbits.PBDIV = 1; // Peripheral Bus 1 Clock Divisor Control (PBCLK1 is SYSCLK divided by 2)

//     // PB2DIV
//     PB2DIVbits.ON = 1; // Peripheral Bus 2 Output Clock Enable (Output clock is enabled)
//     PB2DIVbits.PBDIV = 1; // Peripheral Bus 2 Clock Divisor Control (PBCLK2 is SYSCLK divided by 2)

//     // PB3DIV
//     PB3DIVbits.ON = 1; // Peripheral Bus 2 Output Clock Enable (Output clock is enabled)
//     PB3DIVbits.PBDIV = 1; // Peripheral Bus 3 Clock Divisor Control (PBCLK3 is SYSCLK divided by 2)

//     // PB4DIV
//     PB4DIVbits.ON = 1; // Peripheral Bus 4 Output Clock Enable (Output clock is enabled)
//     while (!PB4DIVbits.PBDIVRDY); // Wait until it is ready to write to
//     PB4DIVbits.PBDIV = 0; // Peripheral Bus 4 Clock Divisor Control (PBCLK4 is SYSCLK divided by 1)

//     // PB5DIV
//     PB5DIVbits.ON = 1; // Peripheral Bus 5 Output Clock Enable (Output clock is enabled)
//     PB5DIVbits.PBDIV = 1; // Peripheral Bus 5 Clock Divisor Control (PBCLK5 is SYSCLK divided by 2)

//     // PB7DIV
//     PB7DIVbits.ON = 1; // Peripheral Bus 7 Output Clock Enable (Output clock is enabled)
//     PB7DIVbits.PBDIV = 0; // Peripheral Bus 7 Clock Divisor Control (PBCLK7 is SYSCLK divided by 1)

//     // PB8DIV
//     PB8DIVbits.ON = 1; // Peripheral Bus 8 Output Clock Enable (Output clock is enabled)
//     PB8DIVbits.PBDIV = 1; // Peripheral Bus 8 Clock Divisor Control (PBCLK8 is SYSCLK divided by 2)

//     // PRECON - Set up prefetch
//     PRECONbits.PFMSECEN = 0; // Flash SEC Interrupt Enable (Do not generate an interrupt when the PFMSEC bit is set)
//     PRECONbits.PREFEN = 3; // Predictive Prefetch Enable (Enable predictive prefetch for any address)
//     PRECONbits.PFMWS = 3; // PFM Access Time Defined in Terms of SYSCLK Wait States (Two wait states)
    
//     //INTCONbits.MVEC = 1;
//     // Set up caching
//     cp0 = _mfc0(16, 0);
//     cp0 &= ~0x07;
//     cp0 |= 0b011; // K0 = Cacheable, non-coherent, write-back, write allocate
//     _mtc0(16, 0, cp0);  

//     // Lock Sequence
//     SYSKEY = 0x33333333;
//     asm volatile("ei"); // Enable all interrupts
// }


// void init(void)
// { 
//     set_performance_mode();

//     // Set all of PORT pins to digital
//     ANSELB = 0;
//     ANSELE = 0;
//     ANSELG = 0;
    
//     DATA_LED_IO;
//     STAT_LED_IO;
// }

// void main(void) 
// {
//     init();
   
//     while(1)
//     {
//        DATA_LED(1); STAT_LED(0); 
//        delay_ms(500);
//        DATA_LED(0); STAT_LED(1); 
//        delay_ms(500);
//     }
// }



