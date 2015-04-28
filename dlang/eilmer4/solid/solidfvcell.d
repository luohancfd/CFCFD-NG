/**
 * solidfvcell.d
 *
 * A solid finite-volume cell, to be held by SolidBlock objects.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-22-04
 */

module solidfvcell;

import geom;
import solidfvinterface;
import solidfvvertex;
import solidprops;

class SolidFVCell {
public:
    // Cell properties
    double volume;
    Vector3 pos;
    SolidProps sp;
    // Cell state
    double[] T;
    double[] e;
    double[] dedt;
    // Connections
    SolidFVInterface[] iface;
    SolidFVVertex[] vtx;

    void timeDerivatives(int ftl, int dimensions)
    {
	SolidFVInterface IFn = iface[Face.north];
	SolidFVInterface IFe = iface[Face.east];
	SolidFVInterface IFs = iface[Face.south];
	SolidFVInterface IFw = iface[Face.west];
	SolidFVInterface IFt, IFb;
	if (dimensions == 3) {
	    IFt = iface[Face.top];
	    IFb = iface[Face.bottom];
	}
	// Cell volume (inverted).
	double volInv = 1.0 / volume;
	double integral;
	
	// Sum up fluxes (of form q.n)
	integral = -IFe.flux * IFe.area - IFn.flux * IFn.area
	    + IFw.flux * IFw.area + IFs.flux * IFs.area;

	dedt[ftl] = volInv * integral;
    }
    void stage1Update(double dt)
    {
	double gamma1 = 1.0;
	e[1] = e[0] + dt*gamma1*dedt[0];
	T[1] = updateTemperature(sp, e[1]);
    }
    void stage2Update(double dt)
    {
	// Assuming predictor-corrector
	double gamma1 = 0.5;
	double gamma2 = 0.5;
	e[2] = e[0] + dt*(gamma1*dedt[0] + gamma2*dedt[1]);
	T[2] = updateTemperature(sp, e[2]);
    }

}
