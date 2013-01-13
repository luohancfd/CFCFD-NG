/** \file spradian_test.cxx
 *  \ingroup radiation
 *  \brief Test the spradian spectral model
 *  \author Daniel F Potter        
 *  \version 16-Apr-2012
 *
 */
#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <time.h>
#include <sstream>

#include "../../gas/models/gas-model.hh"
#include "../../util/source/useful.h"

#include "spectral_model.hh"
#include "spradian.hh"
#include "LOS_pieces.hh"

using namespace std;

int main()
{
    string gas_input_file = "11sp-air-gm.lua";
    string rad_input_file = "spradian-air-EQ-radiators.lua";

    cout << "Creating an 11 species air gas-model from file: " << gas_input_file << endl;
    Gas_model * gm = create_gas_model(gas_input_file);
    
    cout << "Creating an equilibrium air 'spradian07' radiation model from file: " << rad_input_file << endl;
    Spradian * ssm = new Spradian( rad_input_file );
    
    Gas_data Q(gm);
    
    cout << "Calculating post shock conditions for u=10.254 km/s, p=40 Pa air" << endl;
    // For nowjust explicitly set the gas-state
    // NOTE: see python version of this test program for post-shock calculation
    Q.p = 46174.0; Q.T[0] = 10441.0; Q.T[1] = 10441.0;
    Q.massf[0]= 0.00143390000000;
    Q.massf[1]= 6.50260000000e-05;
    Q.massf[2]= 1.23490000000e-06;
    Q.massf[3]= 4.66490000000e-07;
    Q.massf[4]= 6.12270000000e-05;
    Q.massf[5]= 0.000110000000000;
    Q.massf[6]= 0.224290000000;
    Q.massf[7]= 0.00862060000000;
    Q.massf[8]= 0.721530000000;
    Q.massf[9]= 0.0438860000000;
    Q.massf[10]= 2.01780000000e-06; 

    gm->eval_thermo_state_pT(Q);
    Q.print_values(false);
    
    ssm->setup_parameters(Q);
    ssm->get_params_pointer()->write_to_file("spradian.params.1");
    SpradianParams params("spradian.params.1");
    params.write_to_file("spradian.params.2");

    ostringstream oss;
    oss << "diff spradian.params.1 spradian.params.2";
    system(oss.str().c_str());

    CoeffSpectra X;
    ssm->radiative_spectra_for_gas_state(Q,X);

    delete ssm;
    
    return 0;
}
