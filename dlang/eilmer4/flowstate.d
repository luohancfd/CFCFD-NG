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

    this(in GasModel gm, in double p_init, in double T_init[], in Vector3 vel_init,
	 in double[] massf_init=[1.0,], in double quality_init=1.0,
	 in Vector3 B_init=(0.0,0.0,0.0),
	 in double tke_init=0.0, in double omega_init=1.0,
	 in double mu_t_init=0.0, in double k_t_init=0.0,
	 in int S_init=0)
    {
	gas = new GasState(gm, p_init, T_init, massf_init, quality_init);
	vel = vel_init;
	B = B_init;
	tke = tke_init;
	omega = omega_init;
	mu_t = mu_t_init;
	k_t = k_t_init;
	S = S_init;
    }

    this() {} // makes no sense to define the data in the absence of a model

    this(in FlowState other, in GasModel gm)
    {
	gas = new GasState(gm, other.gas.p, other.gas.T, other.gas.massf,
			   other.gas.quality); 
	vel = other.vel;
	B = other.B;
	tke = other.tke;
	omega = other.omega;
	mu_t = other.mu_t;
	k_t = other.k_t;
	S = other.S;
    }

    void copy_values_from(in FlowState other)
    {
	gas.copy_values_from(other.gas);
	vel = other.vel;
	B = other.B;
	tke = other.tke;
	omega = other.omega;
	mu_t = other.mu_t;
	k_t = other.k_t;
	S = other.S;
    }

    // Note that we must not send the current object in the others list as well.
    void copy_average_values_from(in FlowState[] others, in GasModel gm)
    {
	uint n = others.length;
	if (n == 0) throw new Error("Need to average from a nonempty array.");
	GasState[] gasList;
	// We need to be honest and not to fiddle with the other gas states.
	foreach(other; others) {
	    if ( this is other ) throw new Error("Must not include destination in source list.");
	    gasList ~= cast(GasState)other.gas;
	}
	gas.copy_average_values_from(gasList, gm);
	// Accumulate from a clean slate and then divide.
	vel.refx = 0.0; vel.refy = 0.0; vel.refz = 0.0;
	B.refx = 0.0; B.refy = 0.0; B.refz = 0.0;
	tke = 0.0;
	omega = 0.0;
	mu_t = 0.0;
	k_t = 0.0;
	S = 0; // remember that shock detector is an integer flag
	foreach(other; others) {
	    vel += other.vel;
	    B += other.B;
	    tke += other.tke;
	    omega += other.omega;
	    mu_t += other.mu_t;
	    k_t += other.k_t;
	    S += other.S;
	}
	vel /= n;
	B /= n;
	tke /= n;
	omega /= n;
	mu_t /= n;
	k_t /= n;
	S = (S > 0) ? 1 : 0;
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
/+ [TODO]
    double * copy_values_to_buffer(double *buf) const;
    double * copy_values_from_buffer(double *buf);
    int BGK_equilibrium(void);
+/
} // end class FlowState
