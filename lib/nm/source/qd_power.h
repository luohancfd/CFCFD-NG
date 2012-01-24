/** \file qd_power.h
 * \ingroup nm
 * \brief Quick_and_dirty power MACRO for real arguments.
 *
 * Compute z = x**y approximately for moderate values of y and
 * values of x close to 1.
 * The routine produces results within 1% relative error for
 * the following ranges of x (base) and y (exponent).
 * y = 0.2:  0.12 <= x <= 8.4
 * y = 4.0:  0.75 <= x <= 1.55
 * y = -3.0: 0.55 <= x <= 1.48
 *
 * \param x : base
 * \param y : exponent
 * \param z : an approximation to x**y for x near 1 and y small.
 *
 * \version 1.0     17-Mar-91
 *
 * \author PA Jacobs, ICASE
 *
 */

#define QD_POWER(x,y,z)						\
{								\
/* Assumed double x, y, z */					\
double r, t, lnx, temp;						\
								\
/*								\
 * Compute logarithm.						\
 */								\
r = (x - 1.0) / (x + 1.0);					\
temp = r * r * r;						\
lnx = r + 0.333333333 * temp;					\
temp *= r * r;							\
lnx += 0.2 * temp;						\
temp *= r * r;							\
lnx += 0.142857143 * temp;     /* This term = r**7 / 7 */	\
lnx *= 2.0;							\
								\
/*								\
 * Scale by the exponent.					\
 */								\
t = lnx * y;							\
								\
/*								\
 * Compute the exponential.					\
 */								\
temp = 0.5 * t * t;						\
z = 1 + t + temp;						\
temp *= 0.333333333 * t;					\
z += temp;							\
temp *= 0.25 * t;     /* This term = t**4/4!  */		\
z += temp;							\
temp *= 0.2 * t;						\
z += temp;							\
								\
}

