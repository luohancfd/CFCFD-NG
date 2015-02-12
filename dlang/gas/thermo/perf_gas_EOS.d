/**
 * perf_gas_EOS.d
 * Implements the perfect gas equation of state.
 * This module provides simple functions for the
 * the p-v-T behaviour of a perfect gas. It also
 * provides a derived class of ThermalEOS so that
 * the perfect gas model can be used wherever a 
 * ThermalEOS object is valid.
 * 
 * The assumptions of a perfect gas are that:
 *   - the collisions are elastic; and
 *   - the molecules a point masses occupying no volume.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-14 -- first cut
 */

module gas.thermo.perf_gas_EOS;

import gas.gas_model;
import gas.thermo.thermal_EOS;
import util.msg_service;

/++
  Compute the pressure of a perfect gas.
  
  Params:
    rho = The density in kg/m^3
    T = The temperature in K
    R = The specific gas constant in J/(kg.K)

  Returns:
    The pressure in Pa
+/

pure double pressure(double rho, double T, double R)
in {
    assert(rho > 0.0, brokenPreCondition("density", __LINE__, __FILE__));
    assert(T > 0.0, brokenPreCondition("temperature", __LINE__, __FILE__));
    assert(R > 0.0, brokenPreCondition("gas constant", __LINE__, __FILE__));
}
out(result){
    assert(result > 0.0, brokenPostCondition("pressure", __LINE__, __FILE__));
}
body {
    double p = rho*R*T;
    return p;
}

/++
  Compute density for a perfect gas.

  Params:
    p = The pressure in Pa
    T = The temperature in K
    R = The specific gas constant in J/(kg.K)

  Returns:
    The density in kg/m^3
+/
pure double density(double p, double T, double R)
in {
    assert(p > 0.0, brokenPreCondition("pressure", __LINE__, __FILE__));
    assert(T > 0.0, brokenPreCondition("temperature", __LINE__, __FILE__));
    assert(R > 0.0, brokenPreCondition("gas constant", __LINE__, __FILE__));
}
out(result) {
    assert(result > 0.0, brokenPostCondition("density", __LINE__, __FILE__));
}
body {
    double rho = p/(R*T);
    return rho;
}

/++
  Compute the temperature for a perfect gas.
  
  Params:
    rho = The density in kg/m^3
    p = The pressure in Pa
    R = The specific gas constant in J/(kg.K)

  Returns:
    The temperature in K
+/
pure double temperature(double rho, double p, double R)
in {
    assert(rho > 0.0, brokenPreCondition("density", __LINE__, __FILE__));
    assert(p > 0.0, brokenPreCondition("pressure", __LINE__, __FILE__));
    assert(R > 0.0, brokenPreCondition("gas constant", __LINE__, __FILE__));
}
out(result){
    assert(result > 0.0, brokenPostCondition("temperature", __LINE__, __FILE__));
}
body{
    double T = p/(rho*R);
    return T;
}

/++
 PerfectGasEOS is a thermal equation of state.
 
 The perfect gas model assumes point masses and
 elastic collisions.
+/
class PerfectGasEOS : ThermalEOS {
public:
    this(double R) {
	_R = R;
    }

    /++
      Compute the pressure assuming density and temperature
      are up-to-date in GasState Q.
    +/
    override void update_pressure(ref GasState Q) const
    in {
	assert(Q.T.length == 1, brokenPreCondition("Q.T.length", __LINE__, __FILE__));
    }
    body {
	Q.p = pressure(Q.rho, Q.T[0], _R);
    }

    /++
      Compute the density assuming pressure and temperature
      are up-to-date in GasState Q.
    +/
    override void update_density(ref GasState Q) const
    in {
	assert(Q.T.length == 1, brokenPreCondition("Q.T.length", __LINE__, __FILE__));
    }
    body {
	Q.rho = density(Q.p, Q.T[0], _R);
    }

    /++
      Compute the temperature assuming density and pressure
      are up-to-date in GasState Q.
    +/
    override void update_temperature(ref GasState Q) const
    in {
	assert(Q.T.length == 1, brokenPreCondition("Q.T.length", __LINE__, __FILE__));
    }
    body {
	Q.T[0] = temperature(Q.rho, Q.p, _R);
    }

private:
    double _R; /// specific gas constant in J/(kg.K)
}	  

unittest {
    import std.math;
    double R = 287.1;
    assert(approxEqual(pressure(1.2, 300.0, R), 103356.0), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(density(103356.0, 300.0, R), 1.2), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(temperature(1.2, 103356.0, R), 300.0), failedUnitTest(__LINE__, __FILE__));

    auto pg = new PerfectGasEOS(R);
    auto gd = GasState(1, 1);
    gd.T[0] = 300.0;
    gd.rho = 1.2;
    pg.update_pressure(gd);
    assert(approxEqual(gd.p, 103356.0), failedUnitTest(__LINE__, __FILE__));
    gd.p = 103356.0;
    gd.rho = 0.0;
    pg.update_density(gd);
    assert(approxEqual(gd.rho, 1.2), failedUnitTest(__LINE__, __FILE__));
    gd.rho = 1.2;
    gd.T[0] = 0.0;
    pg.update_temperature(gd);
    assert(approxEqual(gd.T[0], 300.0), failedUnitTest(__LINE__, __FILE__));
}



	   


