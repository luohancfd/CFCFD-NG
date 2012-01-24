/** \file inv_roberts.c
 * \ingroup nm
 * \brief Compute the required stretching to achieve a specified cell width.
 * \author Ian or Andrew (?) 
 */

#include <stdio.h>
#include <math.h>
#include "../../util/source/compiler.h"
#include "roberts.h"

int main()
{  
int i;
double eta;
double alpha;
double beta;
double lambda, etabar, etabar_act;
double l,n,h;

printf ("Enter the distane: "); fflush(stdout);
scanf ("%lf", &l);

printf ("Enter the number of cells: "); fflush(stdout);
scanf ("%lf", &n);

printf ("Enter the first height: "); fflush(stdout);
scanf ("%lf", &h);

etabar_act = 1.0 - (h/l);
eta = 1.0 - (1.0 / n);
alpha = 0.0; /* nodes will be clustered near eta=1.0 */
beta = 10.0;
for (i=1; i < 9000000; ++i)
    {
    lambda = (beta + 1.0) / (beta - 1.0);
    lambda = pow ( lambda, ((eta - alpha)/(1.0 - alpha)) );
    etabar = (beta + 2.0 * alpha) * lambda - beta + 2.0 * alpha;
    etabar = etabar / ((2.0 * alpha + 1.0) * (1.0 + lambda));
    if (etabar < etabar_act) 
        beta = 10.0 - (i * 1.0e-6);
    }
beta = beta - 1.0e-6;
printf (" The required Beta value is %lf\n",beta);
}


/*=================================================================*/
