#ifndef _DRX_HELPERS_CPP_
#define _DRX_HELPERS_CPP_

#include "helpers.h"

inline float gamma_l(float misorientation) {
    if (misorientation == 0) {
        return 0;
    }
    else if (misorientation >= CRITICAL_MISORIENTATION) {
        return gamma_o;
    } else {
        return gamma_o * (misorientation / CRITICAL_MISORIENTATION) * (1 - log(misorientation / CRITICAL_MISORIENTATION));
    }
}

inline float mobility(float misorientation) {
    if (misorientation == 0) { 
        return 0;
    } else {
        return Mo * (1 - exp(-5 * pow((misorientation / CRITICAL_MISORIENTATION), 4)));
    }
}

void print_array(int array[GRID_SIZE][GRID_SIZE]) {
    int i = 0, j = 0;
    ff(i, 0, GRID_SIZE) {
        ff(j, 0, GRID_SIZE) {
            cout << array[i][j] << " ";
        }
        cout << "\n";
    }
}

#endif