/*
 * therm_cond.d
 * Interface for all thermal conductivity models.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-19 -- initial cut
 */

import gasmodel;

interface ThermalConductivity {
    ThermalConductivity dup() const;
    final void update_thermal_conductivity(ref GasState Q)
    {
	for ( auto imode = 0; imode < Q.T.length; ++imode) {
	    Q.k[imode] = eval(Q, imode);
	}
    }
    double eval(in GasState Q, int imode);
}
