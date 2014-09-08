/**
 * sutherland_therm_cond.d
 * Implements the Sutherland "law" for
 * thermal conductivity/
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-19 -- initial cut
 */

import std.math;
import gasmodel;
import therm_cond;

/++
  Compute the thermal conductivity using Sutherland's expression.
  
  Params:
     T = temperature of gas in K
     T_ref = reference temperature in K
     k_ref = reference thermal conductivity in W/(m.K)
     S = Sutherland constant in K

   Returns:
     The thermal conductivity in W/(m.K)
+/
pure double sutherland_thermal_conductivity(double T, double T_ref, double k_ref, double S)
in {
    assert(T > 0.0);
    assert(T_ref > 0.0);
    assert(k_ref > 0.0);
}
out(result) {
    assert(result > 0.0);
}
body{
    double k = k_ref*sqrt(T/T_ref)*(T/T_ref)*(T_ref + S)/(T + S);
    return k;
}

/++
  SutherlandThermalConductivity is a thermal conductivity model.
+/
class SutherlandThermalConductivity : ThermalConductivity {
public:
    this(double T_ref, double k_ref, double S) {
	_T_ref = T_ref;
	_k_ref = k_ref;
	_S = S;
    }

    /++
      Compute the thermal conductivity assuming that temperature is
      up-to-date in GasState Q.
    +/
    override void update_thermal_conductivity(ref GasState Q) const {
	assert(Q.T.length == 1);
	assert(Q.k.length == 1);
	Q.k[0] = sutherland_thermal_conductivity(Q.T[0], _T_ref, _k_ref, _S);
    }

private:
    double _T_ref;
    double _k_ref;
    double _S;
}

unittest {
    double T = 300.0;
    double T_ref = 273.0; 
    double k_ref = 0.0241;
    double S = 194.0;
    assert(approxEqual(sutherland_thermal_conductivity(T, T_ref, k_ref, S), 0.0262449));

    auto tcm = new SutherlandThermalConductivity(T_ref, k_ref, S);
    auto gd = GasState(1, 1);
    gd.T[0] = 300.0;
    gd.k[0] = 0.0;
    tcm.update_thermal_conductivity(gd);
    assert(approxEqual(gd.k[0], 0.0262449));
    //    assert(approxEqual(gd.k[0], 0.249));
}
