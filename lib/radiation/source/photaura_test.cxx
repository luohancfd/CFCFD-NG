/** \file photaura_test.cxx
 *  \ingroup radiation
 *  \brief Test the photaura spectral model
 *  \author Daniel F Potter        
 *  \version 26-Aug-2009
 *
 */
#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <time.h>

#include "../../gas/models/gas-model.hh"
#include "../../util/source/useful.h"

#include "spectral_model.hh"
#include "LOS_pieces.hh"

using namespace std;

int main()
{
    string gas_input_file = "11sp-air-gm.lua";
    string rad_input_file = "air-radiators.lua";

    cout << "Creating an 11 species air gas-model from file: " << gas_input_file << endl;
    Gas_model * gm = create_gas_model(gas_input_file);
    
    cout << "Creating an equilibrium air 'photaura' radiation model from file: " << rad_input_file << endl;
    RadiationSpectralModel * rsm = create_radiation_spectral_model( rad_input_file );
    
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
    
    int n_evals = 100;
    cout << "Performing " << n_evals << " emission coefficient evaluations" << endl;
    clock_t c0 = clock();
    double j_total;
    for( int i=0; i<n_evals; i++ )
    	j_total = rsm->radiative_integrated_emission_for_gas_state( Q, false );
    clock_t c1 = clock();
    cout << "j_total = " << ( j_total * 1.0e-6 )  << " W/cm**3" << endl;
    cout << "elapsed CPU [spectrally_resolved=false]:          " 
    	 <<  (float) (c1 - c0)/CLOCKS_PER_SEC << endl;
    	 
    cout << "Calculating emission and absorption spectra for gas state" << endl;
    CoeffSpectra X;
    rsm->radiative_spectra_for_gas_state( Q, X );
    j_total = X.write_to_file( "coefficient_spectra.txt" );
    cout << "Spectrally integrated emission: j_total = " << ( j_total * 1.0e-6 )  << " W/cm**3" << endl;
    
    TS_data TS(rsm,2);
    double Q_rE_rad;
    TS.set_rad_point(0,&Q,&Q_rE_rad,-0.05,0.1);
    TS.set_rad_point(1,&Q,&Q_rE_rad, 0.05,0.1);
    cout << "testing quick_solve_for_divq" << endl;
    double q_rad = TS.quick_solve_for_divq();
    cout << "q_rad = " << q_rad << " W/cm2" << endl;
    cout << "testing quick_solve_for_divq_with_binning(OPACITY_BINNING, N_bins=10)" << endl;
    q_rad = TS.quick_solve_for_divq_with_binning( OPACITY_BINNING, 10 );
    cout << "q_rad = " << q_rad << " W/cm2" << endl;
    cout << "testing quick_solve_for_divq_with_binning(OPACITY_BINNING, N_bins=100)" << endl;
    q_rad = TS.quick_solve_for_divq_with_binning( OPACITY_BINNING, 100 );
    cout << "q_rad = " << q_rad << " W/cm2" << endl;
    cout << "testing quick_solve_for_divq_with_binning( FREQUENCY_BINNING, N_bins=10)" << endl;
    q_rad = TS.quick_solve_for_divq_with_binning( FREQUENCY_BINNING, 10 );
    cout << "q_rad = " << q_rad << " W/cm2" << endl;
    cout << "testing quick_solve_for_divq_with_binning( FREQUENCY_BINNING, N_bins=95)" << endl;
    q_rad = TS.quick_solve_for_divq_with_binning( FREQUENCY_BINNING, 95 );
    cout << "q_rad = " << q_rad << " W/cm2" << endl;
    cout << "testing quick_solve_for_divq_with_binning( FREQUENCY_BINNING, N_bins=950)" << endl;
    q_rad = TS.quick_solve_for_divq_with_binning( FREQUENCY_BINNING, 950 );
    cout << "q_rad = " << q_rad << " W/cm2" << endl;
    cout << "testing quick_solve_for_divq_with_binning( FREQUENCY_BINNING, N_bins=9500)" << endl;
    q_rad = TS.quick_solve_for_divq_with_binning( FREQUENCY_BINNING, 9500 );
    cout << "q_rad = " << q_rad << " W/cm2" << endl;
    cout << "testing quick_solve_for_divq_with_binning( FREQUENCY_BINNING, N_bins=95000)" << endl;
    q_rad = TS.quick_solve_for_divq_with_binning( FREQUENCY_BINNING, 95000 );
    cout << "q_rad = " << q_rad << " W/cm2" << endl;


    cout << "Clearing gas and radiation models" << endl;
    delete gm;
    delete rsm;
    
    cout << "Finished photaura spectral model testing." << endl;
    
    return 0;
}
