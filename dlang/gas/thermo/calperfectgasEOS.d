/**
 * calperfectgasEOS.d
 * Implements a calorically perfect gas equation
 * of state. In this model, the specific heats
 * at constant pressure (Cp) and constant volume (Cv)
 * are constant for all temperatures.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-18 -- first cut
 */

import gasmodel;
import caloricEOS;

/++
  Compute the energy of a calorically perfect gas.
  
  Params:
     T = temperature in K
     Cv = specific heat at constant volume in J/(kg.K)
     e0 = reference energy value in J/kg 

  Returns:
     The internal energy in J/kg.
+/
pure double energy(double T, double Cv, double e0) {
    assert(T >= 0.0);
    assert(Cv > 0.0);
    double e = Cv*T + e0;
    return e;
}

/++
  Compute the temperature of a calorically perfect gas.

  Params:
     e = internal energy in J/kg
     Cv = specific heat at constant volume in J/(kg.K)
     e0 = reference energy value in J/kg

  Returns:
     The temperature in K.

+/
pure double temperature(double e, double Cv, double e0) {
    assert(Cv > 0.0);
    double T = (e - e0)/Cv;
    return T;
}

/++
  CaloricallyPerfectGasEOS is a caloric equation of state.

  This model assumes that Cp and Cv are constant with temperature.
+/
class CaloricallyPerfectGasEOS : CaloricEOS {
public:
    this(double Cv, double e0) {
	_Cv = Cv;
	_e0 = e0;
    }

    /++
      Compute the internal energy assuming that temperature
      is up-to-date in GasState Q.
    +/
    override void update_energy(ref GasState Q) const {
	assert(Q.T.length == 1);
	assert(Q.e.length == 1);
	Q.e[0] = energy(Q.T[0], _Cv, _e0);
    }

    /++ Compute the temperature assuming that the internal energy
     is up-to-date in GasState Q.
    +/
    override void update_temperature(ref GasState Q) const {
	assert(Q.T.length == 1);
	assert(Q.e.length == 1);
	Q.T[0] = temperature(Q.e[0], _Cv, _e0);
    }

private:
    double _Cv; /// specific heat at constant volume in J/(kg.K)
    double _e0; /// reference energy value in J/kg
}

unittest {
    import std.math;
    double Cv = 717.0;
    double e0 = 0.0;
    assert(approxEqual(energy(500.0, Cv, e0), 358500.0));
    assert(approxEqual(temperature(358500.0, Cv, e0), 500.0));

    auto cpg = new CaloricallyPerfectGasEOS(Cv, e0);
    auto gd = GasState(1, 1);
    gd.T[0] = 500.0;
    gd.e[0] = 0.0;
    cpg.update_energy(gd);
    assert(approxEqual(gd.e[0], 358500.0));
    gd.e[0] = 358500.0;
    gd.T[0] = 0.0;
    cpg.update_temperature(gd);
    assert(approxEqual(gd.T[0], 500.0));
}
