/**
 * flowstate.d
 * Flow condition class for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module flowstate;

import std.string;
import std.conv;
import geom;
import gasmodel;
import gasmodelutil;

class FlowState {
public:
    GasState gas;  // gas state
    Vector3 vel;   // flow velocity, m/s
    Vector3 B;     // magnetic field strength
    double tke;    // turbulent kinetic energy 0.5(u'^2+v'^2+w'^2)
    double omega;  // turbulence 'frequency' in k-omega model
    double mu_t;   // turbulence viscosity
    double k_t;    // turbulence thermal-conductivity
    int S;         // shock indicator, value 0 or 1

    this(GasModel gm, double p_init, double T_init[], Vector3 vel_init,
	 double[] massf_init=[1.0,], double quality_init=1.0,
	 Vector3 B_init=(0.0,0.0,0.0),
	 double tke_init=0.0, double omega_init=1.0,
	 double mu_t_init=0.0, double k_t_init=0.0,
	 int S_init=0)
    {
	_gm = gm;
	gas = new GasState(_gm, p_init, T_init, massf_init, quality_init);
	vel = vel_init;
	B = B_init;
	tke = tke_init;
	omega = omega_init;
	mu_t = mu_t_init;
	k_t = k_t_init;
	S = S_init;
    }

    this() {} // makes no sense to define the data in the absence of a model

    this(FlowState other)
    {
	_gm = other._gm;
	gas = new GasState(other._gm, other.gas.p, other.gas.T, other.gas.massf,
			   other.gas.quality); 
	vel = other.vel;
	B = other.B;
	tke = other.tke;
	omega = other.omega;
	mu_t = other.mu_t;
	k_t = other.k_t;
	S = other.S;
    }

    override string toString()
    {
	char[] repr;
	repr ~= "FlowState(";
	repr ~= "gas=" ~ to!string(gas);
	repr ~= ", vel=" ~ to!string(vel);
	repr ~= ", B=" ~ to!string(B);
	repr ~= ", tke=" ~ to!string(tke);
	repr ~= ", omega=" ~ to!string(omega);
	repr ~= ", mu_t=" ~ to!string(mu_t);
	repr ~= ", k_t=" ~ to!string(k_t);
	repr ~= ", S=" ~ to!string(S);
	repr ~= ")";
	return to!string(repr);
    }

private:
    GasModel _gm;  // don't want clients to be digging the gas model out of here.
}
