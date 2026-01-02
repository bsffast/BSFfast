#include <stdio.h>
#include <stdlib.h>
#include "BSFfast.h"

int main(void) {
    double m = 1000.;   // example mass value (GeV)
    double x = 13.;     // example x value

    // The library’s constructor (init_library) runs automatically,
    // so you don’t need to call anything to load the data.

	double result = BSFfast_sigveff_QCD_SU_plat(m, x);
    printf("BSFfast_sigveff_QCD_SU_plat(m=%e, x=%e) = %e\n", m, x, result);

    // Example of another function:
    double alpha = 1/128.9;
    double result2 = BSFfast_sigveff_dQED_S_const(alpha, m, x);

    printf("BSFfast_sigveff_QED_S_const(alpha=%e, m=%e, x=%.5f) = %e\n",
           alpha, m, x, result2);

    return 0;
}

// $ gcc -o testBSFfast BSFfast-test.c BSFfast.c -lm
// $ ./testBSFfast
