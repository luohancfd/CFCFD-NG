/**
 * cal_perf_gas_mix_EOS.d
 * Implements a calorically perfect gas mixture equation
 * of state. In this model, the specific heats
 * at constant pressure (Cp) and constant volume (Cv)
 * are functions of the gas composition. The specific heats
 * for each gas component are constant for all temperatures.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-09-07 -- first cut
 */

module gas.thermo.cal_perf_gas_mix_EOS;

import gas.gas_model;
import gas.thermo.caloric_EOS;
import gas.thermo.cal_perf_gas_EOS;

/++
  CaloricallyPerfectGasMixEOS is a caloric equation of state.

  This model assumes that Cp and Cv are constant with temperature
  for all gas components.
+/
class CaloricallyPerfectGasMixEOS : CaloricEOS {
public:
    this(double[] Cv, double[] e0)
    {
	_Cv = Cv.dup;
	_e0 = e0.dup;
    }

    /++
      Compute the internal energy assuming that temperature
      is up-to-date in GasState Q.
    +/
    override void update_energy(ref GasState Q) const {
	assert(Q.T.length == 1);
	assert(Q.e.length == 1);
	double Cv = mass_average(Q, _Cv);
	double e0 = mass_average(Q, _e0);
	Q.e[0] = energy(Q.T[0], Cv, e0);
    }

    /++ Compute the temperature assuming that the internal energy
     is up-to-date in GasState Q.
    +/
    override void update_temperature(ref GasState Q) const {
	assert(Q.T.length == 1);
	assert(Q.e.length == 1);
	double Cv = mass_average(Q, _Cv);
	double e0 = mass_average(Q, _e0);
	Q.T[0] = temperature(Q.e[0], Cv, e0);
    }

private:
    double[] _Cv; /// specific heats at constant volume in J/(kg.K)
    double[] _e0; /// reference energies value in J/kg
}

unittest {
    import std.math;
    double[] Cv = [743.0, 659.0];
    double[] e0 = [0.0, 0.0];
    auto cpg = new CaloricallyPerfectGasMixEOS(Cv, e0);
    auto gd = GasState(2, 1);
    gd.T[0] = 500.0;
    gd.e[0] = 0.0;
    gd.massf[0] = 0.78;
    gd.massf[1] = 0.22;
    cpg.update_energy(gd);
    assert(approxEqual(gd.e[0], 362260.0));
    gd.e[0] = 362260.0;
    gd.T[0] = 0.0;
    cpg.update_temperature(gd);
    assert(approxEqual(gd.T[0], 500.0));
}
