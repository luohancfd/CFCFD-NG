/**
 * conservedquantities.d
 * Class for the vector of conserved quantities, for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module conservedquantities;
import geom;
import gasmodel;

class ConservedQuantities {
public:
    double mass;         // density, kg/m**3
    Vector3 momentum;    // momentum/unit volume
    Vector3 B;           // magnetic field, Tesla
    double total_energy; // total energy
    double[] massf;      // mass fractions of species
    double[] energies;   // modal energies (mode 0 is usually transrotational)
    double tke;          // turbulent kinetic energy
    double omega;        // omega from k-omega turbulence model
    // [TODO] double[] G;          // velocity dist. partial densities, kg/m**3
    // [TODO] double[] H;          // velocity dist. partial densities, (kg*s**2)/(m**5)

    this(in GasModel gm)
    {
	massf.length = gm.n_species;
	energies.length = gm.n_modes;
    }

    this(in ConservedQuantities other)
    {
	mass = other.mass;
	momentum = other.momentum;
	B = other.B;
	total_energy = other.total_energy;
	massf = other.massf.dup;
	energies = other.energies.dup;
	tke = other.tke;
	omega = other.omega;
    }

    void copy_values_from(in ConservedQuantities src)
    {
	mass = src.mass;
	momentum = src.momentum;
	B = src.B;
	total_energy = src.total_energy;
	massf[] = src.massf[];
	energies[] = src.energies[];
	tke = src.tke;
	omega = src.omega;
    }

    void clear_values()
    {
	mass = 0.0;
	momentum.x = 0.0; momentum.y = 0.0; momentum.z = 0.0;
	B.x = 0.0; B.y = 0.0; B.z = 0.0;
	total_energy = 0.0;
	foreach(ref mf; massf) mf = 0.0;
	foreach(ref e; energies) e = 0.0;
	tke = 0.0;
	omega = 0.0;
    }
}
