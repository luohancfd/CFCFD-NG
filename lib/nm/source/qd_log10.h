/** \file qd_log10.h
 * \ingroup nm
 * \brief Quick_and_dirty logarithm MACRO for real arguments.
 *
 * Compute z = log10(x) approximately for values of x close to 1.
 *
 * \version 1.0     09-Apr-92
 *
 * \author PA Jacobs
 *
 */

#define QD_LOG10(x,z)						\
{								\
/* Assumed double x, z */					\
double r2, r, t, lnx, temp;					\
								\
/*								\
 * Compute natural logarithm.					\
 */								\
r = (x - 1.0) / (x + 1.0);					\
r2 = r * r;							\
temp = r * r2;							\
lnx = r + 0.333333333 * temp;					\
temp *= r2;							\
lnx += 0.2 * temp;						\
temp *= r2;							\
lnx += 0.142857143 * temp;     /* This term = r**7 / 7 */	\
lnx *= 2.0;							\
								\
/*								\
 * Scale to get logarithm to base 10.				\
 */								\
z = 0.434294482 * lnx;						\
								\
}
