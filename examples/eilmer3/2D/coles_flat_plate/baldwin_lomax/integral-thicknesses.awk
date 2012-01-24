# momentum-deficit.awk
# Compute integral thicknesses for the boundary layer as per AGARD report.
# PJ, October 2007

BEGIN {
    rho_D = 5.81e-2;  # density, kg/m**3 
    u_D = 678.0;      # velocity, m/s
    delta = 0.017;    # boundary-layer edge, distance from wall in metres
    y_offset = 0.24;  # y-position of plate, m
    rhou_D = rho_D * u_D;
    # Subscript D indicates quantity at the outer edge of the boundary layer.
    d_1 = 0.0;        # displacement thickness
    d_2 = 0.0;        # momentum-deficit
    integrating = 0; 
}

$1 != "#" {
    # We are approaching the wall from outside the boundary layer.
    y = y_offset - $2;
    rho = $5;
    u = $6;
    if ( integrating == 0 && y < delta ) integrating = 1;
    if ( integrating ) {
	print "y=", y, "rho=", rho, "u=", u;
	rhou = 0.5 * (rho * u + rho_old * u_old);
	dy = (y_old - y);
	d_1 += (1.0 - rhou/rhou_D) * dy;
	d_2 += rhou/rhou_D * (1.0 - rhou/rhou_D) * dy;
    }
    y_old = y;
    rho_old = rho;
    u_old = u;
}

END {
    print "d_1=", d_1;
    print "d_2=", d_2;
}
