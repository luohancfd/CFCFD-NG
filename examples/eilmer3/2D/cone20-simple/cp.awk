# cp.awk
# Extract the simulation times and axial force values from the log file.
# The relevant lines in mb_cns.log start with the string "XFORCE"
# and are of the form:
#    XFORCE: t n jb ibndy fx_p fx_v [jb ibndy fx_p fx_v [jb ...]]
# Present the axial force as an average coefficient of pressure to
# compare with that obtained from NACA 1135.

BEGIN {
    p_inf = 95.84e3;  # Pa
    T_inf = 1103.0;   # K
    u_inf = 1000.0;   # m/s
    R = 287;          # J/kg.K
    r_base = 0.29118; # m
    rho_inf = p_inf / (R * T_inf);  # kg/m**3
    q_inf = 0.5 * rho_inf * u_inf * u_inf;  # Pa
    A = 3.14159 * r_base * r_base;  # m**2
    print "# time (ms)  Cp";
    print "# rho_inf= ", rho_inf, " q_inf= ", q_inf, " A= ", A
}

/XFORCE/ {
    # Select just the simulation time and the force on the cone surface.
    t = $3;  # in seconds
    f = $9;  # pressure force in Newtons
    # The coefficient of pressure is based on the difference 
    # between the cone surface pressure and the free-stream pressure.
    Cp = (f / A - p_inf) / q_inf;
    print t*1000.0, Cp;
}
