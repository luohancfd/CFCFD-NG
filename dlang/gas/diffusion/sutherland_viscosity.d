/**
 * sutherland_viscosity.d
 * Implements the Sutherland "law" for
 * viscosity.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-19 -- initial cut
 */

import std.math;
import gasmodel;
import viscosity;

/++
  Compute the viscosity using Sutherland's expression.
  
  Params:
     T = temperature of gas in K
     T_ref = reference temperature in K
     mu_ref = reference viscosity in Pa.s
     S = Sutherland constant in K

   Returns:
     The viscosity in Pa.s.
+/
pure double sutherland_viscosity(double T, double T_ref, double mu_ref, double S)
in {
    assert(T > 0.0);
    assert(T_ref > 0.0);
    assert(mu_ref > 0.0);
}
out(result) {
    assert(result > 0.0);
}
body{
    double mu = mu_ref*sqrt(T/T_ref)*(T/T_ref)*(T_ref + S)/(T + S);
    return mu;
}

/++
  SutherlandViscosity is a viscosity model.
+/
class SutherlandViscosity : Viscosity {
public:
    this(double T_ref, double mu_ref, double S) {
	_T_ref = T_ref;
	_mu_ref = mu_ref;
	_S = S;
    }

    /++
      Compute the viscosity assuming that temperature is
      up-to-date in GasState Q.
    +/
    override void update_viscosity(ref GasState Q) const {
	assert(Q.T.length >= 1);
	Q.mu = sutherland_viscosity(Q.T[0], _T_ref, _mu_ref, _S);
    }

private:
    double _T_ref;
    double _mu_ref;
    double _S;
}

unittest {
    double T = 300.0;
    double T_ref = 273.0; 
    double mu_ref = 1.716e-5;
    double S = 111.0;
    assert(approxEqual(sutherland_viscosity(T, T_ref, mu_ref, S), 1.84691e-05));

    auto vm = new SutherlandViscosity(T_ref, mu_ref, S);
    auto gd = new GasState();
    gd.T ~= 300.0;
    vm.update_viscosity(gd);
    assert(approxEqual(gd.mu, 1.84691e-05));
}
