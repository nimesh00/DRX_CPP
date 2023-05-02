#ifndef _DRX_HELPERS_
#define _DRX_HELPERS_
#include <math.h>
// #include "drx_grid.h"

using namespace std;

#define STATES 50
#define GRID_SIZE 200
#define ITERATIONS_MC 1000000
#define MAX_ORIENTATION 90.0

#define DIS_DEN_MEAN 8.0e11
#define DIS_DEN_STD_DEV 1.0e11
#define DIS_DEN_RESET 2.022e14
#define DIS_DEN_CRITICAL 8.08e14
#define DIS_DEN_MAX 8.08e14

#define CELL_SIZE 2.0e-6
#define CRITICAL_MISORIENTATION 15 // degrees
#define EPS_CR 0.22
#define EPS_FINAL 1.5
#define EPS_RATE 0.01
#define T 725.0 // K
#define Tm 1356.0 // K
#define R 8.3145 // J-mol^-1-K^-1
#define C 3.5e22
#define gamma_m 0.625 // J-m^-2
#define r_gamma 0.66
// #define Mo 3.6e-5 // As in Hallberg
#define Mo 1.4e-11 // As in Liu
#define Qm 126e3 // J-mol^-1
#define Qn 261e3 // J-mol^-1
#define alpha 0.5
#define b 2.56e-10 // m
#define nu 0.3 // poisson's ratio
// #define mu_o 3e10 // Pa -- As in Hallberg
#define mu_o 7.89e10 // Pa -- As in Liu
#define k1 2.78e8
#define k2 10.09
#define ff(i, s, e) for (int i = s; i < e; i++)
#define fb(i, e, s) for (int i = e; i >= s; i--)

// float cell_gamma = 0;
float nucleation_rate = C * EPS_RATE * exp(-Qn / (R * T));
float M = Mo * exp(-1 * Qm / (R * T));
float mu = mu_o * (1 - 0.5 * (T - 300.0) / (Tm));
float gamma_o = (mu * b * CRITICAL_MISORIENTATION) / (4 * M_PI * (1 - nu));
float tau = 0.5 * mu * b * b;
float p_ciritcal = DIS_DEN_MEAN;

int encoder = 0;

float gamma_l(float);

float mobility(float);

inline bool cell_on_border(int, int);

void print_array(int array[GRID_SIZE][GRID_SIZE]);

#endif
