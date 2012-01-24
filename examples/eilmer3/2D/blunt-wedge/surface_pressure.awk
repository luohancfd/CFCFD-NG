# surface_pressure.awk
# Normalise the surface pressure with free-stream dynamic pressure and
# compute the distance around from the stagnation point.

BEGIN {
    q_inf = 1.750e6;   # free-stream dynamic pressure, Pa
    p_inf = 100.0e3;   # free-stream static pressure, Pa
    Rn = 10.0e-3;      # nose radius
    xold = -Rn;        # location of the stagnation point
    yold = 0.0;
    s = 0.0;           # distance around from stagnation point
    count = 0;
    pi = 3.1415927;
    wedge_angle = 10.0/180.0 * pi;
    print "# s/Rn  Cp(CFD)  Cp(Newton)  x(m)  y(m)";
}

$1 != "#" {
    count += 1;
    x = $1;            # cell-centre position
    y = $2;
    p = $9;            # cell-centre pressure
    if ( count == 1 ) p_pitot = p;  # Close enough to the stagnation point.
    dx = x - xold;
    dy = y - yold;
    s += sqrt(dx * dx + dy * dy);
    # Estimate Cp using Modified Newtonian Model.
    theta = 0.5 * pi - (s/Rn);  # local angle of surface
    if (theta < wedge_angle) theta = wedge_angle;
    Cp_MN = (p_pitot - p_inf) / q_inf * sin(theta) * sin(theta);
    print s/Rn, (p - p_inf)/q_inf, Cp_MN, x, y;
    xold = x;
    yold = y;
}
