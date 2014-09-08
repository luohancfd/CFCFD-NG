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
    void update_thermal_conductivity(ref GasState Q);
}
