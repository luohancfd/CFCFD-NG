# momentum-deficit.awk
# Compute integral thicknesses for the boundary layer as per AGARD report.
# Wilson Chan, 11-Jul-2011
#  modified from same file in coles-flat-plate example by PJ, October 2007

BEGIN {
    rho_D = 0.17707;  # density, kg/m**3 
    u_D = 712.9;      # velocity, m/s
    delta = 0.00485;  # boundary-layer edge, distance from wall in metres
    y_offset = 0.16;  # y-position of plate, m
    rhou_D = rho_D * u_D;
    # Subscript D indicates quantity at the outer edge of the boundary layer.
    displacement_thickness = 0.0;
    momentum_thickness = 0.0;
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
	displacement_thickness += (1.0 - rhou/rhou_D) * dy;
	momentum_thickness += rhou/rhou_D * (1.0 - rhou/rhou_D) * dy;
    }
    y_old = y;
    rho_old = rho;
    u_old = u;
}

END {
    print "displacement_thickness_metres=", displacement_thickness;
    print "momentum_thickness_metres=", momentum_thickness;
}
