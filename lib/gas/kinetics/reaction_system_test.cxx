/*  \file reaction_system_test.cxx
 *  \brief Testing program for the finite-rate chemistry module.
 *
 *  \author Rowan J Gollan
 *  \version 25-Feb-2006
 *
 */

/** \file tpgm_test.cxx
 * \brief Testing program of the thermally perfect gas mix.
 * \author RJG
 * \version 14-Feb-06
 *
 **/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "physical_constants.hh"
#include "gas.hh"
#include "rr_coeffs.hh"
#include "eq_const_model.hh"
#include "reaction.hh"
#include "../../nm/source/ode_solver.hh"
#include "reaction_scheme.hh"


using namespace std;

void printUsage()
{
    cout << "Usage: reaction_test.x [--verbose|--fp-only|--data]" << endl;
    cout << "Options: \n";
    cout << "    --help      : Print this help.\n";
    cout << "    --verbose   : Print tests in human-readable form (default action)\n";
    cout << "    --fp-only   : Same as verbose but only floating-point values\n";
    cout << "                  are printed.  This is useful for regression comparisons.\n";
    cout << "    --data      : Generate data over large T range for plotting in an output file\n";
    return;
}

void run_standard_test( bool verbose );
void generate_data();

int main( int argc, char *argv[] )
{
    if( argc > 2 ) {
	printUsage();
	return 0;
    }

    if( argc == 1 ) {
	run_standard_test( true );
	return 0;
    }

    string action( argv[1] );

    if( action == "--help" ) {
	printUsage();
	return 0;
    }
    else if( action == "--verbose" ) {
	run_standard_test( true );
	return 0;
    }
    else if( action == "--fp-only" ) {
	run_standard_test( false );
	return 0;
    }
    else if( action == "--data" ) {
	generate_data();
	return 0;
    }
    else {
	cout << "Uknown option: " << action << endl;
	printUsage();
	return 0;
    }

    return 0;

}

void run_standard_test( bool verbose )
{

    cout << setprecision(10) << showpoint;

    if( verbose ) {
	cout << "Begin test of finite-rate chemistry module.\n";
	cout << "-------------------------------------------------------------------\n";
	cout << " Test 1: Adiabatic chemical reactor - constant density and energy  \n";
	cout << "-------------------------------------------------------------------\n";
	cout << "We consider a mixture of N2 and N with the following properties.   \n";
	cout << "p = 100.0 kPa (~ 1 atm)\n";
	cout << "T = 4000.0 K \n";
	cout << "f[N2] = 0.8 \n";
	cout << "f[N] = 0.2 \n";
    }

    // Setup the gas mixture
    set_type_of_gas( string("perf_gas_mix"), string("nitrogen-2sp.pgm") );
    GasModel *g = get_gas_model_pointer();
    Gas_data Q(g);

    Q.T = 4000.0;
    Q.p = 1.0e5;
    Q.f[0] = 0.8;
    Q.f[1] = 0.2;
    Q.dt_chem = -1;
    
    g->EOS_pT( Q, false );
    if( verbose ) {
	cout << "After initial call to EOS_pT... \n";
    }
    print_gas_data( Q );

    // Setup the reaction description.
    // Reaction 1
    GeneralisedArrhenius f_rate1( "_f", 7e+15, -1.6, 113200.0*PC_k_SI );
    GeneralisedArrhenius b_rate1( "_b", 10900, -0.5, 0.0 );
    vector<int> index; index.push_back( 0 ); index.push_back( 1 );
    vector<int> nu; nu.push_back( -1 ); nu.push_back( 2 );
    RateCoeffModel e_rate( "FromEqConst", "_b");
    EqConstThermo e_1;
    e_1.set_index( index );
    e_1.set_nu( nu );
    
    Dissociation reac1("N2-dissociation", "N2 + N2 <=> N + N + N2",
		       &f_rate1, &e_rate, &e_1, 0, 1, 1, 0);

    
    GeneralisedArrhenius f_rate2( "_f", 3e+16, -1.6, 113200.0*PC_k_SI );
    GeneralisedArrhenius b_rate2( "_b", 2.32e+09, -1.5, 0.0 );

    Dissociation reac2("N2-dissociation", "N2 + N <=> N + N + N",
		       &f_rate2, &e_rate, &e_1, 0, 1, 1, 1);

    OdeSolver ode_solver( "reacting system", 2, "qss", 4, 1.15, 0.01, 0.33 );

    ReactionSchemeODE r_system("Dissociation of nitrogen", 2, 0.1, 0.0, &ode_solver );
    r_system.add_reaction( &reac1 );
    r_system.add_reaction( &reac2 );

    double t = 0.0;
    double dt_flow = 1.0e-7;



    while( t < 2.0e-4 ) {
 	r_system.update_gas_state( Q, dt_flow );
 	t += dt_flow;
    }

    if( verbose ) {
	cout << "After call to update_gas_state... \n";
    }

    print_gas_data( Q );
    
//     vector<int> f_index; f_index.push_back(0); f_index.push_back(0);
//     vector<int> f_coeffs; f_coeffs.push_back(1); f_index.push_back(1);
//     vector<int> b_index; b_index.push_back(1); b_index.push_back(1); b_index.push_back(0);
//     vector<int> b_coeffs; b_coeffs.push_back(1); b_coeffs.push_back(1); b_coeffs.push_back(1);

//     GeneralReaction reac3("N2-dissociation", "N2 + N2 <=> N + N + N2",
// 			  &f_rate1, &b_rate1, f_index, f_coeffs, b_index, b_coeffs);
    
//     f_index.clear(); f_index.push_back(0); f_index.push_back(1);
//     f_coeffs.clear(); f_coeffs.push_back(1); f_index.push_back(1);
//     b_index.clear(); b_index.push_back(1); b_index.push_back(1); b_index.push_back(1);

//     GeneralReaction reac4( "N2-dissociation", "N2 + N <=> N + N + N",
// 			   &f_rate2, &b_rate2,
// 			   f_index, f_coeffs, b_index, b_coeffs );

//     reactions.clear();
//     reactions.push_back( &reac1 );
//     reactions.push_back( &reac2 );

//     ReactionSchemeODE r_system2("Dissociation of nitrogen", reactions, 2, 0.1, 0.0, &ode_solver );


//     t = 0.0;
//     dt_flow = 1.0e-7;
//     Q.T = 4000.0;
//     Q.p = 1.0e5;
//     Q.f[0] = 0.8;
//     Q.f[1] = 0.2;
//     Q.dt_chem = -1.0;
    
//     g->EOS_pT( Q, false );

//     if( verbose ) {
// 	cout << "After re-call to EOS_pT... \n";
//     }
    
//     print_gas_data( Q );

//    //  while( t < 2.0e-4 ) {
// // 	r_system2.update_gas_state( Q, dt_flow );
// // 	t += dt_flow;
// //     }
//     if( verbose ) {
// 	cout << "Now describing the reactions by GeneralReaction types.\n";
// 	cout << "After call to update_gas_state... \n";
//     }

//     print_gas_data( Q );

//     if( verbose ) {
// 	cout << "-------------------------------------------------------------------\n";
// 	cout << " Test 2: Hydrogen combustion (adiabatic chemical reactor)          \n";
// 	cout << "-------------------------------------------------------------------\n";
// 	cout << "We consider a stoichiometric mixture of H2 and O2 with the         \n";
// 	cout << "following conditions:                                              \n";
// 	cout << "p = 100.0 kPa (~ 1 atm) \n";
// 	cout << "T = 1200.0 K            \n";
// 	cout << "f[H2] = 0.1189875       \n";
// 	cout << "f[O2] = 0.8810125       \n";
//     }
    
//     // Setup the gas mixture
//     clear_gas_model_pointer();
//     set_type_of_gas( string("perf_gas_mix"), string("h2-combustion-7sp.pgm") );
//     g = get_gas_model_pointer();

//     Q.T = 2000.0;
//     Q.p = 1.0e5;
//     Q.f[0] = 0.0;        Q.f[1] = 0.0;
//     Q.f[2] = 0.1189875;  Q.f[3] = 0.8810125;
//     Q.f[4] = 0.0;        Q.f[5] = 0.0;       Q.f[6] = 0.0;
//     Q.dt_chem = -1.0;
//     g->EOS_pT(Q, false );
    
    
//     // Setup the reaction description.
//     // 1: HO2 + H <=> OH + OH
//     GeneralisedArrhenius frate1( "_f", 6.0e7, 0.0, 0.0 );
//     GeneralisedArrhenius brate1( "_b", 1.7e5, 0.5, 21137.0*PC_k_SI );

//     f_index.clear(); f_index.push_back(1); f_index.push_back(5);
//     f_coeffs.clear(); f_coeffs.push_back(1); f_coeffs.push_back(1);
//     b_index.clear(); b_index.push_back(4); b_index.push_back(4);
//     b_coeffs.clear(); b_coeffs.push_back(1); b_coeffs.push_back(1);

//     GeneralReaction r1( "HO2-decomposition-I", "HO2 + H <=> OH + OH",
// 			&frate1, &brate1, f_index, f_coeffs, b_index, b_coeffs);

//     // 2: HO2 + H <=> H2O + O
//     GeneralisedArrhenius frate2( "_f", 3.0e7, 0.0, 0.0 );
//     GeneralisedArrhenius brate2( "_b", 5.8e5, 0.5, 28686.0*PC_k_SI );

//     f_index.clear(); f_index.push_back(1); f_index.push_back(5);
//     f_coeffs.clear(); f_coeffs.push_back(1); f_coeffs.push_back(1);
//     b_index.clear(); b_index.push_back(0); b_index.push_back(6);
//     b_coeffs.clear(); b_coeffs.push_back(1); b_coeffs.push_back(1);
    
//     GeneralReaction r2( "HO2-decomposition-II", "HO2 + H <=> H2O + O",
// 			&frate2, &brate2, f_index, f_coeffs, b_index, b_coeffs);

//     // 3: H2 + M <=> H + H + M
//     GeneralisedArrhenius frate3( "_f", 5.5e12, -1.0, 51987.0*PC_k_SI);
//     GeneralisedArrhenius brate3( "_b", 1.8e6, -1.0, 0.0 );
    
//     f_index.clear(); f_index.push_back(2);
//     f_coeffs.clear(); f_coeffs.push_back(1);
//     b_index.clear(); b_index.push_back(5); b_index.push_back(5);
//     b_coeffs.clear(); b_coeffs.push_back(1); b_coeffs.push_back(1);

//     ThirdBodyReaction r3( "H+H-recombination", "H2 + M <=> H + H + M",
// 			  &frate3, &brate3, f_index, f_coeffs, b_index, b_coeffs );

//     // 4: O2 + M <=> O + O + M
//     GeneralisedArrhenius frate4( "_f", 7.2e12, -1.0, 59340.0*PC_k_SI);
//     GeneralisedArrhenius brate4( "_b", 4.0e5,  -1.0, 0.0 );

//     f_index.clear(); f_index.push_back(3);
//     f_coeffs.clear(); f_coeffs.push_back(1);
//     b_index.clear(); b_index.push_back(6); b_index.push_back(6);
//     b_coeffs.clear(); b_coeffs.push_back(1); b_coeffs.push_back(1);

//     ThirdBodyReaction r4( "O+O-recombination", "O2 + M <=> O + O + M",
// 			  &frate4, &brate4, f_index, f_coeffs, b_index, b_coeffs );
    
//     // 5: H2O + M <=> OH + H + M
//     GeneralisedArrhenius frate5( "_f", 5.2e15, -1.5, 59386.0*PC_k_SI);
//     GeneralisedArrhenius brate5( "_b", 4.4e8,  -1.5, 0.0 );

//     f_index.clear(); f_index.push_back(0);
//     f_coeffs.clear(); f_coeffs.push_back(1);
//     b_index.clear(); b_index.push_back(4); b_index.push_back(5);
//     b_coeffs.clear(); b_coeffs.push_back(1); b_coeffs.push_back(1);

//     ThirdBodyReaction r5( "OH+H-recombination", "H2O + M <=> OH + H + M",
// 			  &frate5, &brate5, f_index, f_coeffs, b_index, b_coeffs );


//     reactions.clear();
//     reactions.push_back(&r1); reactions.push_back(&r2); reactions.push_back(&r3);
//     reactions.push_back(&r4); reactions.push_back(&r5);

//     OdeSolver ode_solver2( "combusting system", 7, "euler", 4, 1.15, 0.01, 0.33 );

//     ReactionSchemeODE r_system3("Combustion of hydrogen", reactions, 7, 0.1, 0.0, &ode_solver2 );

//     print_gas_data( Q );
//     t = 0.0;
//     dt_flow = 1.0e-10;
//     while( t < 6.0e-5 ) {
// 	r_system3.update_gas_state( Q, dt_flow );
// 	t += dt_flow;
//     }
    
//     print_gas_data( Q );
//     for( int isp = 0; isp < g->nsp; ++isp ) {
// 	cout << "mole_f[" << isp <<"]= " << Q.c[isp]/Q.c_tot << endl;
//     }
    clear_gas_model_pointer();
    return;
}



void generate_data()
{
    cout << "Begin generating data for the reaction_system_test...\n";
    // Setup the gas mixture
    gas_data Q, Q_e;
    set_type_of_gas( string("perf_gas_mix"), string("nitrogen-2sp.pgm") );
    GasModel *g = get_gas_model_pointer();
    set_array_sizes_in_gas_data(Q, 2, 0);
    set_array_sizes_in_gas_data(Q_e, 2, 0);

    Q.T = 4000.0;
    Q.p = 1.0e5;
    Q.f[0] = 0.8;
    Q.f[1] = 0.2;
    Q.dt_chem = -1;
    
    Q_e.T = 4000.0;
    Q_e.p = 1.0e5;
    Q_e.f[0] = 0.8;
    Q_e.f[1] = 0.2;
    Q_e.dt_chem = -1;

    g->EOS_pT( Q, false );
    g->EOS_pT( Q_e, false );

    // Setup the reaction description.
    // Reaction 1
    vector<int> index; index.push_back( 0 ); index.push_back( 1 );
    vector<int> nu; nu.push_back( -1 ); nu.push_back( 2 );
    GeneralisedArrhenius f_rate1( "_f", 7e+15, -1.6, 113200.0*PC_k_SI );
    GeneralisedArrhenius b_rate1( "_b", 10900, -0.5, 0.0 );
    RateCoeffModel e_rate( "FromEqConst", "_b");
    EqConstThermo e_1;
    e_1.set_index( index );
    e_1.set_nu( nu );

    Dissociation reac1("N2-dissociation", "N2 + N2 <=> N + N + N2",
		       &f_rate1, &b_rate1, 0, 0, 1, 1, 0);

    Dissociation reac1e("N2-dissociation", "N2 + N2 <=> N + N + N2",
			&f_rate1, &e_rate, &e_1, 0, 1, 1, 0);


    // Reaction 2
    GeneralisedArrhenius f_rate2( "_f", 3e+16, -1.6, 113200.0*PC_k_SI );
    GeneralisedArrhenius b_rate2( "_b", 2.32e+09, -1.5, 0.0 );

    Dissociation reac2("N2-dissociation", "N2 + N <=> N + N + N",
		       &f_rate2, &b_rate2, 0, 0, 1, 1, 1);
    Dissociation reac2e("N2-dissociation", "N2 + N <=> N + N + N",
			&f_rate2, &e_rate, &e_1, 0, 1, 1, 1);
    vector<Reaction*> reactions;
    reactions.push_back( &reac1 );
    reactions.push_back( &reac2 );
    
    vector<Reaction*> reactions_e;
    reactions_e.push_back( &reac1e );
    reactions_e.push_back( &reac2e );
   
    OdeSolver ode_solver( "reacting system", 2, "qss", 4, 1.15, 0.01, 0.33 );

    ReactionSchemeODE r_system("Dissociation of nitrogen", reactions, 2, 0.1, 0.0, &ode_solver );
    ReactionSchemeODE r_system_e("Dissociation of nitrogen", reactions_e, 2, 0.1, 0.0, &ode_solver );

    double t = 0.0;
    double dt_flow = 1.0e-7;
    ofstream outfile( "nitrogen-system-history.data" );
    outfile << setprecision(10) << showpoint;
    outfile << "# nitrogen-system-history.data\n"
	    << "# Columns:\n"
	    << "# 1: t \n"
	    << "# 2: T \n"
	    << "# 3: f[0] \n"
	    << "# 4: f[1] \n"
	    << "# 5: f[0] - eq calc\n"
	    << "# 6: f[1] - eq calc\n";
	
    while( t < 2.0e-4 ) {
	r_system.update_gas_state( Q, dt_flow );
 	r_system_e.update_gas_state( Q_e, dt_flow );
	
	t += dt_flow;
	outfile << t << "  " << Q.T << "  " << Q.f[0] << "  " << Q.f[1] << "  " << Q_e.f[0] << "  " << Q_e.f[1] << endl;
    }

    outfile.close();
    cout << "File: nitrogen-system-history.data created.\n";
    cout << "Done.\n";
    clear_gas_model_pointer();

    return;
}
