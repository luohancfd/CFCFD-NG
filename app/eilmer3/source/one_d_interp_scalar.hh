// one_d_interp_scalar.hh

#ifndef ONE_D_INTERP_SCALAR_HH
#define ONE_D_INTERP_SCALAR_HH

#include "../../../lib/util/source/useful.h"

/// \brief One-dimensional reconstruction of a scalar quantity.
///
/// See MBCNS workbook 2000/2 page 36 (26-Jan-2001) for formulation.
/// and MBCNS workbook 2005/Apr page 36 for new index labels
///
inline int one_d_interp_scalar( double qL1, double qL0, double qR0, double qR1, 
				double lenL1, double lenL0, double lenR0, double lenR1, 
				double &qL, double &qR, int apply_limiter,
				int extrema_clipping )
{
    double aL0, aR0, delLminus, del, delRplus, sL, sR;
    double lower_limit, upper_limit;
    // Set up differences and limiter values.
    aL0 = 0.5 * lenL0 / (lenL1 + 2.0*lenL0 + lenR0);
    aR0 = 0.5 * lenR0 / (lenL0 + 2.0*lenR0 + lenR1);
    delLminus = 2.0 * (qL0 - qL1) / (lenL0 + lenL1);
    del = 2.0 * (qR0 - qL0) / (lenR0 + lenL0);
    delRplus = 2.0 * (qR1 - qR0) / (lenR1 + lenR0);
    if ( apply_limiter == 1 ) {
	// val Albada limiter as per Ian Johnston's thesis.
#       define EPSILON 1.0e-12
	sL = (delLminus * del + fabs(delLminus * del)) /
	    (delLminus * delLminus + del * del + EPSILON);
	sR = (del * delRplus + fabs(del * delRplus)) /
	    (del * del + delRplus * delRplus + EPSILON);
    } else {
	sL = 1.0;
	sR = 1.0;
    }
    // The high-order reconstruction, possibly limited.
    qL = qL0 + sL * aL0 * ( del * (2.0*lenL0 + lenL1) + delLminus * lenR0 );
    qR = qR0 - sR * aR0 * ( delRplus * lenL0 + del * (2.0*lenR0 + lenR1) );
    if ( extrema_clipping == 1 ) {
	// An extra limiting filter to make sure that we have not introduced
	// any new extreme values.
	// This was introduced to deal with very sharp transitions in species.
	lower_limit = MINIMUM(qL0, qR0);
	upper_limit = MAXIMUM(qL0, qR0);
	qL = MINIMUM(upper_limit, MAXIMUM(lower_limit, qL)); 
	qR = MINIMUM(upper_limit, MAXIMUM(lower_limit, qR));
    }
    return SUCCESS;
} // end of one_d_interp_scalar()

#endif
