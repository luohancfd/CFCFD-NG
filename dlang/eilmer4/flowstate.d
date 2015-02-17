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
import std.json;
import json_helper;
import geom;
import gas;

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

    this(in FlowState other, in GasModel gm)
    {
	gas = new GasState(gm, other.gas.p, other.gas.T, other.gas.massf, other.gas.quality); 
	vel = other.vel;
	B = other.B;
	tke = other.tke;
	omega = other.omega;
	mu_t = other.mu_t;
	k_t = other.k_t;
	S = other.S;
    }

    this(in GasModel gm)
    {
	gas = new GasState(gm, 100.0e3, [300.0,], [1.0,], 1.0); 
	vel = Vector3(0.0,0.0,0.0);
	B = Vector3(0.0,0.0,0.0);
	tke = 0.0;
	omega = 1.0;
	mu_t = 0.0;
	k_t = 0.0;
	S = 0;
    }

    this(in JSONValue json_data, in GasModel gm)
    {
	double p = getJSONdouble(json_data, "p", 100.0e3);
	double T[] = getJSONdoublearray(json_data, "T", [300.0,]);
	double[] massf = getJSONdoublearray(json_data, "massf", [1.0,]);
	double quality = 1.0;
	gas = new GasState(gm, p, T, massf, quality);
	double u = getJSONdouble(json_data, "u", 0.0);
	double v = getJSONdouble(json_data, "v", 0.0);
	double w = getJSONdouble(json_data, "w", 0.0);
	vel = Vector3(u,v,w);
	double Bx = getJSONdouble(json_data, "Bx", 0.0);
	double By = getJSONdouble(json_data, "By", 0.0);
	double Bz = getJSONdouble(json_data, "Bz", 0.0);
	B = Vector3(Bx,By,Bz);
	tke = getJSONdouble(json_data, "tke", 0.0);
	omega = getJSONdouble(json_data, "omega", 1.0);
	mu_t = getJSONdouble(json_data, "mu_t", 0.0);
	k_t = getJSONdouble(json_data, "k_t", 0.0);
	S = getJSONint(json_data, "S", 0);
    }

    this() {} // makes no sense to define the data in the absence of a model

    @nogc void copy_values_from(in FlowState other)
    {
	gas.copy_values_from(other.gas);
	vel.refx = other.vel.x; vel.refy = other.vel.y; vel.refz = other.vel.z;
	B.refx = other.B.x; B.refy = other.B.y; B.refz = other.B.z;
	tke = other.tke;
	omega = other.omega;
	mu_t = other.mu_t;
	k_t = other.k_t;
	S = other.S;
    }

    @nogc void copy_average_values_from(in FlowState fs0, in FlowState fs1)
    // Avoids memory allocation, it's all in place.
    {
	gas.copy_average_values_from(fs0.gas, fs1.gas);
	vel.refx = 0.5 * (fs0.vel.x + fs1.vel.x);
	vel.refy = 0.5 * (fs0.vel.y + fs1.vel.y);
	vel.refz = 0.5 * (fs0.vel.z + fs1.vel.z);
	B.refx = 0.5 * (fs0.B.x + fs1.B.x);
	B.refy = 0.5 * (fs0.B.y + fs1.B.y);
	B.refz = 0.5 * (fs0.B.z + fs1.B.z);
	tke = 0.5 * (fs0.tke + fs1.tke);
	omega = 0.5 * (fs0.omega + fs1.omega);
	mu_t = 0.5 * (fs0.mu_t + fs1.mu_t);
	k_t = 0.5 * (fs0.k_t + fs1.k_t);
    } // end copy_average_values_from()

    void copy_average_values_from(in FlowState[] others, in GasModel gm)
    // Note that we must not send the current object in the others list as well.
    // Involves some memory allocation.
    {
	size_t n = others.length;
	if (n == 0) throw new Error("Need to average from a nonempty array.");
	GasState[] gasList;
	// Note that, because we cast away their "const"ness,
	// we need to be honest and not to fiddle with the other gas states.
	foreach(other; others) {
	    if ( this is other ) {
		throw new Error("Must not include destination in source list.");
	    }
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
	S = 0; // Remember that shock detector is an integer flag.
	foreach(other; others) {
	    vel.refx += other.vel.x;
	    vel.refy += other.vel.y;
	    vel.refz += other.vel.z;
	    B.refx += other.B.x;
	    B.refx += other.B.x;
	    B.refx += other.B.x;
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
    } // end copy_average_values_from()

    override string toString() const
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
