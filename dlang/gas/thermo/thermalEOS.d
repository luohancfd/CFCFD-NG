/**
 * thermalEOS.d
 * Interface for all thermal equations of state.
 *
 * The thermal equation of state relates the 
 * pressure to the density and temperature.
 * This is also referred to the p-v-T behaviour
 * of the gas, where v is the specific volume,
 * the inverse of density.
 *
 * The thermal equation of state is distinct
 * from the caloric equation of state. The latter
 * relates the internal energy of the gas to the
 * density and temperature. The two equations of
 * state together can be used to specify the 
 * complete thermodynamic state of the gas given
 * two state variables.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-14 -- initial cut at building an interface
 */

import gasmodel;

/++
  ThermalEOS defines the servies provided by a thermal equation
  of state model.

  All thermal equations of state define a relationship between
  the pressure, temperature and density of a gas.
+/
interface ThermalEOS {
    void update_pressure(ref GasState Q) const;
    void update_density(ref GasState Q) const;
    void update_temperature(ref GasState Q) const;
}
