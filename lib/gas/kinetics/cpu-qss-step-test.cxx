// Authors: Rowan J. Gollan and Kyle Damm
// Date: 11-Mar-2015
//

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <sys/time.h>
#include <cstdlib>

#include "../models/gas-model.hh"
#include "chemical-kinetic-ODE-update.hh"
#include "../../nm/source/ode_step.hh"

using namespace std;

void print_conc(vector<double> &conc)
{
    cout << setprecision(16);
    cout << scientific;
    for ( size_t i = 0; i < conc.size(); ++i ) {
	cout << conc[i] << endl;
    }
}

int main()
{
    // --- Program settings ----
    int nsteps = 1;
    double dt_fixed = 1.0e-10;
    // alpha-QSS settings
    int max_correctors = 4;
    double eps1 = 1.0e-3;
    double c = 2.0;
    double delta = 1.0e-10;
    // --- End: Program settings ---

    Gas_model *gmodel = create_gas_model("gas-model.lua");
    Chemical_kinetic_system *cks = create_chemical_kinetic_system("Evans_Schexnayder.lua", *gmodel);
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();
    QssStep qss_step("qss-step", nsp, max_correctors, eps1, c, delta);

    // --- Initital conditions ---
    vector<double> molef;
    molef.resize(nsp);
    // Moles total = 3.0, stoichiometric ratio -- O2:1, H2:2
    // O2 is 0 in list, H2 is 2 in list
    molef[0] = 1./3.; molef[2] = 2./3.;
    Gas_data Q(gmodel);
    Q.p = 1.0e5;
    Q.T[0] = 1000.0;
    convert_molef2massf(molef, gmodel->M(), Q.massf);
    gmodel->eval_thermo_state_pT(Q);
    cout << "Initial conditions:\n";
    Q.print_values(false);
    // --- End: initial conditions ---

    // --- Main stepping ---
    vector<double> C0, C1;
    C0.resize(nsp); C1.resize(nsp);
    convert_massf2conc(Q.rho, Q.massf, gmodel->M(), C0);
    cks->set_gas_data_ptr(Q);
    
    cout << "Initial concentration:\n";
    print_conc(C0);
    double dt;
    for ( size_t step = 0; step < nsteps; ++step ) {
	dt = dt_fixed;
	if ( !qss_step.advance(*cks, C0, C1, &dt) ) {
	    cout << "Step failed with dt= " << dt_fixed << endl;
	    cout << "Try smaller dt_fixed.\n";
	    exit(1);
	}
	// Copy C1 into C0 for next step
	for ( size_t i = 0; i < C0.size(); ++i ) C0[i] = C1[i];
    }
    cout << "After " << nsteps << " steps :\n";
    print_conc(C1);

    return 0;
}
