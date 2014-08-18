/**
 * caloricEOS.d
 * Interface for all caloric equations of state.
 *
 * The caloric equation of state relates the
 * internal energy of the gas to the density
 * and temperature.
 *
 * See the counterpart interface file thermalEOS.d
 * for the thermal equation of state (p-v-T behaviour).
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-18 -- initial cut at building the interface
 */

import gasmodel;

/++
  CaloricEOS defines the services provied by a caloric equation
  of state model.

  All caloric equations of state define a relationship between
  the internal energy, temperature and density of a gas.

+/
interface CaloricEOS {
    void update_energy(ref GasState Q) const;
    void update_temperature(ref GasState Q) const;
}

