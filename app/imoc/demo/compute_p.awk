# Compute the pressure from Mach number assuming isentropic processing
# of the gas.
# PJ, 25-Feb-2000

function pow( base, ex ) {
   lnb = log(base)
   result = exp(ex * lnb)
   return result
}

BEGIN {
    p1 = 370121.4
    M1 = 2.7226
    g  = 1.4
    ex = g / (g - 1.0)
    base  = 1.0 + 0.5 * (g - 1.0) * M1 * M1
    p0_p1 = pow( base, ex ) 
    print "# x (m), p (Pa)" 
}

$1 != "#" {
    x = $1
    M = $2
    base = 1.0 + 0.5 * (g - 1.0) * M * M
    p0_p = pow( base, ex )
    p = p1 * p0_p1 / p0_p
    print x, p
}
 
