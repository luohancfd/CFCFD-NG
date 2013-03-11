/** \file photaura.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 28-Jul-09 : initial framework
 *
 **/
 
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "photaura.hh"
#include "radiation_constants.hh"
#include "atomic_radiator.hh"
#include "electron_radiator.hh"
#include "../../util/source/lua_service.hh"

using namespace std;

Photaura::
Photaura( lua_State * L )
 : RadiationSpectralModel( L )
{
    lua_getfield(L, -1, "radiators" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Photaura::Photaura():\n";
	ost << "Error in the declaration of radiators: a table is expected.\n";
	input_error(ost);
    }
    
    nrad = lua_objlen(L, -1);
    vector<const char*> rad_names;
    
    for ( int irad = 0; irad < nrad; ++irad ) {
	lua_rawgeti(L, -1, irad+1); // A Lua list is offset one from the C++ vector index
	const char* rad = luaL_checkstring(L, -1);
	rad_names.push_back(rad);
	lua_pop(L, 1);
    }
    
    lua_pop(L,1);	// pop radiators
	
    // Now construct the radiators
    int e_index = -1;
    int e_irad = -1;
    for ( int irad = 0; irad < nrad; ++irad ) {
	radiators.push_back( create_new_radiator( L, rad_names[irad] ) );
	if ( radiators.back()->type == "electron_radiator" ) {
	    e_index = radiators.back()->isp;
	    e_irad = irad;
	}
    }
    
    // Set e_index in all the radiators
    for ( int irad = 0; irad < nrad; ++irad ) {
    	radiators[irad]->set_e_index( e_index );
    }
    
    // Initialise continuum mechanisms now that all radiators have been created
    if ( e_irad!=-1 ) {
    	ElectronRadiator * e = dynamic_cast<ElectronRadiator*>(radiators[e_irad]);
    	e->create_continuum_mechanisms( radiators );
    }
    
    // Set the radiator pointers as required by QSS and FirstOrderLTNE radiators
    for ( int irad = 0; irad < nrad; ++irad ) {
    	if ( radiators[irad]->type=="atomic_radiator" && radiators[irad]->EPM=="QSS" ) {
    	    dynamic_cast<QSSAtomicRadiator*>(radiators[irad])->set_radiator_pointers( radiators );
	}
	else if ( radiators[irad]->type=="atomic_radiator" && radiators[irad]->EPM=="FirstOrderLTNE" ) {
	    dynamic_cast<FirstOrderLTNEAtomicRadiator*>(radiators[irad])->set_radiator_pointers( radiators );
	}
	else if ( radiators[irad]->type=="diatomic_radiator" && radiators[irad]->EPM=="QSS" ) {
	     dynamic_cast<QSSDiatomicRadiator*>(radiators[irad])->set_radiator_pointers( radiators );
	}
    }
}

Photaura::
Photaura( const string input_file )
 : RadiationSpectralModel( input_file )
{
    // 1. Get spectral_model string from lua file
    lua_State *L = initialise_radiation_lua_State();

    if( luaL_dofile(L, input_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "Photaura::Photaura()\n";
	ost << "Error in input file: " << input_file << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, "spectral_data" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Photaura::Photaura():\n";
	ost << "Error in the declaration of spectral_data: a table is expected.\n";
	input_error(ost);
    }
    
    lua_getfield(L, -1, "radiators" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Photaura::Photaura():\n";
	ost << "Error in the declaration of radiators: a table is expected.\n";
	input_error(ost);
    }
    
    nrad = lua_objlen(L, -1);
    vector<const char*> rad_names;
    
    for ( int irad = 0; irad < nrad; ++irad ) {
	lua_rawgeti(L, -1, irad+1); // A Lua list is offset one from the C++ vector index
	const char* rad = luaL_checkstring(L, -1);
	rad_names.push_back(rad);
	lua_pop(L, 1);
    }
    
    lua_pop(L,1);	// pop radiators
	
    // Now construct the radiators
    int e_index = -1;
    int e_irad = -1;
    for ( int irad = 0; irad < nrad; ++irad ) {
	radiators.push_back( create_new_radiator( L, rad_names[irad] ) );
	if ( radiators.back()->type == "electron_radiator" ) {
	    e_index = radiators.back()->isp;
	    e_irad = irad;
	}
    }
    
    // Set e_index in all the radiators
    for ( int irad = 0; irad < nrad; ++irad ) {
    	radiators[irad]->set_e_index( e_index );
    }
    
    // Initialise continuum mechanisms now that all radiators have been created
    if ( e_irad!=-1 ) {
    	ElectronRadiator * e = dynamic_cast<ElectronRadiator*>(radiators[e_irad]);
    	e->create_continuum_mechanisms( radiators );
    }
    
    // Set the radiator pointers as required by QSS and FirstOrderLTNE radiators
    for ( int irad = 0; irad < nrad; ++irad ) {
    	if ( radiators[irad]->type=="atomic_radiator" && radiators[irad]->EPM=="QSS" ) {
    	    dynamic_cast<QSSAtomicRadiator*>(radiators[irad])->set_radiator_pointers( radiators );
	}
	else if ( radiators[irad]->type=="atomic_radiator" && radiators[irad]->EPM=="FirstOrderLTNE" ) {
	    dynamic_cast<FirstOrderLTNEAtomicRadiator*>(radiators[irad])->set_radiator_pointers( radiators );
	}
	else if ( radiators[irad]->type=="diatomic_radiator" && radiators[irad]->EPM=="QSS" ) {
	    dynamic_cast<QSSDiatomicRadiator*>(radiators[irad])->set_radiator_pointers( radiators );
	}
    }
    
    lua_pop(L,1);	// pop spectral_data
    lua_close(L);
}

Photaura::
~Photaura()
{
    // 1. Delete radiators
    for ( size_t irad=0; irad<radiators.size(); ++irad ) 
	delete radiators[irad];
}

string Photaura::str() const
{
    return "Photaura";
}

Radiator *
Photaura::
get_radiator_pointer_from_name( string name )
{
    Radiator * R = 0;
    
    for ( int irad=0; irad<nrad; ++irad ) {
    	if ( radiators[irad]->name == name ) {
    	    R = radiators[irad];
    	    break;
    	}
    }
    
    if ( !R ) {
    	cout << "Photaura::get_radiator_pointer_from_name()" << endl
    	     << "Radiator with name: " << name << " not found." << endl
    	     << "Exiting program." << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    return R;
}

DiatomicRadiator *
Photaura::
get_diatomic_radiator_pointer_from_name( string name )
{
    Radiator * R = get_radiator_pointer_from_name( name );
    
    if ( R->type!="diatomic_radiator" ) {
    	cout << "Photaura::get_diatomic_radiator_pointer_from_name()" << endl
    	     << "Radiator with name: " << name << " has type: " << R->type << endl
    	     << "Exiting program." << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    return dynamic_cast<DiatomicRadiator*>(R);
}

AtomicRadiator *
Photaura::
get_atomic_radiator_pointer_from_name( string name )
{
    Radiator * R = get_radiator_pointer_from_name( name );
    
    if ( R->type!="atomic_radiator" ) {
    	cout << "Photaura::get_atomic_radiator_pointer_from_name()" << endl
    	     << "Radiator with name: " << name << " has type: " << R->type << endl
    	     << "Exiting program." << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    return dynamic_cast<AtomicRadiator*>(R);
}

double
Photaura::
integrate_emission_spectrum_for_radiator( Gas_data &Q, int irad, bool write_to_file )
{
    // 0. Create a CoeffSpectra class to perform the spectral calculations
    CoeffSpectra S;
    
    // 1. Initialise the radiator
    radiators[irad]->calc_partition_functions(Q);
    radiators[irad]->calc_elec_pops(Q);
    radiators[irad]->init_mechanisms(Q);
    
    // 2. Store the radiation spectrum for the specified gas-state
    spectral_distribution_for_gas_state(Q,S.nu);
    
    // 3. Resize j_nu and kappa_nu vectors (and initialise to zero)
    S.j_nu.resize(S.nu.size(),0.0);
    S.kappa_nu.resize(S.nu.size(),0.0);
    
    // 4. Spectral calculation for the requested radiator
    radiators[irad]->calc_spectrum(Q,S);
    
    // 2. Integrate the emission coefficient over frequency using trapezoidal method
    double j_int=0.0;
    for( int inu=1; inu<int(S.nu.size()); ++inu ) {
        j_int += 0.5 * ( S.j_nu[inu] + S.j_nu[inu-1] ) * ( S.nu[inu] - S.nu[inu-1] );
    }
    
    if ( write_to_file ) S.write_to_file(radiators[irad]->name + "-emission-spectra.txt");
    
    return j_int;
}

double
Photaura::
integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved )
{
    double j_total = 0.0;
    
    if ( spectrally_resolved == true ) {
    	j_total = integrate_emission_spectrum(Q);
    }
    else {
    	j_total = sum_emission_coefficients(Q);
    }
    
    return j_total;
}

double
Photaura::
variably_integrated_emission_for_gas_state( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u, bool spectrally_resolved )
{
    double j_total = 0.0;
    
    if ( spectrally_resolved == true ) {
    	cout << "radiative_variably_integrated_emission_for_gas_state()" << endl
    	     << "Not presently available with spectrally resolved option." << endl
    	     << "bailing out!" << endl;
    	exit( BAD_INPUT_ERROR );
    }
    else {
    	j_total = sum_optically_variable_emission_coefficients(Q,wavel_switch,Lambda_l,Lambda_u);
    }
    
    return j_total;
}

double
Photaura::
integrate_emission_spectrum( Gas_data &Q )
{
    // 0. Create a CoeffSpectra class to perform the spectral calculations
    CoeffSpectra S;
    
    // 1. Store the radiation spectrum for the specified gas-state
    spectra_for_gas_state( Q, S );
    if ( DEBUG_RAD > 0 ) cout << "S.nu.size() = " << S.nu.size() << endl;
    
    // 2. Integrate the emission coefficient over frequency using trapezoidal method
    double j_int=0.0;
    for( int inu=1; inu<int(S.nu.size()); ++inu ) {
        j_int += 0.5 * ( S.j_nu[inu] + S.j_nu[inu-1] ) * ( S.nu[inu] - S.nu[inu-1] );
    }
    
    return j_int;
}

double
Photaura::
sum_emission_coefficients( Gas_data &Q )
{
    /* 1. Initialise all radiator mechanisms (partition functions, populations, linewidths) */
    // NOTE: linewidths are not used here
    initialise_all_radiators(Q);
    
    /* 2. Sum emission coefficients */
    double j_total = 0.0;
    for ( int irad=0; irad<nrad; ++irad)  {
    	if ( get_rad_conc(Q,irad) > MIN_CONC ) {
    	    j_total += radiators[irad]->calc_unresolved_emission_coefficient( Q );
    	}
    }
    
    return j_total;
}

double
Photaura::
sum_optically_variable_emission_coefficients( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u )
{
    /* 1. Initialise all radiator mechanisms (partition functions, populations, linewidths) */
    // NOTE: linewidths are not used here
    initialise_all_radiators(Q);
    
    /* 2. Sum emission coefficients */
    double j_total = 0.0;
    for ( int irad=0; irad<nrad; ++irad)  {
    	if ( get_rad_conc(Q,irad) > MIN_CONC ) {
    	    j_total += radiators[irad]->calc_unresolved_OV_emission_coefficient( Q, wavel_switch, Lambda_l, Lambda_u );
    	}
    }
    
    return j_total;
}

void
Photaura::
spectra_for_gas_state( Gas_data &Q, CoeffSpectra &X )
{
    /* 1. Initialise all radiator mechanisms (partition functions, populations, linewidths) */
    initialise_all_radiators(Q);
    
    /* 2. Calculate the frequency distribution to be used for the spectrum. */
    spectral_distribution_for_gas_state(Q,X.nu);
#   if DEBUG_RAD > 0
    cout << "Photaura::spectra_for_gas_state() - X.nu.size() = " << X.nu.size() << endl;
    cout << "lambda_min = " << nu2lambda(X.nu.back()) << ", lambda_max = " << nu2lambda(X.nu.front()) << endl;
#   endif

    /* 2a. Clear, resize and initialise to zero the j_nu and kappa_nu vectors */
    X.j_nu.clear();
    X.j_nu.resize(X.nu.size(),0.0);
    X.kappa_nu.clear();
    X.kappa_nu.resize(X.nu.size(),0.0);
    
    /* 3. Proceed with spectral emission and absorption calculations */
    for ( int irad=0; irad<nrad; ++irad)  {
    	// cout << "calculating spectra for radiator: " << radiators[irad]->name << endl;
    	if ( get_rad_conc(Q,irad) > MIN_CONC ) {
    	    radiators[irad]->calc_spectrum(Q,X);
    	}
#       if DEBUG_RAD > 0
	else {
	    cout << "Radiator: " << radiators[irad]->name 
	         << " has been skipped due to insufficient concentration: " 
	         << get_rad_conc(Q,irad) << " moles" << endl;
	}
#       endif
    }
    
    return;
}

void
Photaura::
spectral_distribution_for_gas_state(Gas_data &Q, vector<double> &nus)
{
    if ( adaptive_spectral_grid ) {
        /* Obtain an optimised spectral distribution */
        nus.reserve( this->get_spectral_points() );
        for( int irad=0; irad<nrad; ++irad ) {
            // impose lower concentration limit
            if ( get_rad_conc( Q, irad ) > MIN_CONC) {
                radiators[irad]->get_spectral_distribution(nus);
            }
        }
        /* Sort the entries into ascending order */
        sort(nus.begin(), nus.end(), less<double>());
        /* Impose spectral limits */
        nus.erase( nus.begin(), lower_bound( nus.begin(), nus.end(), lambda2nu( this->get_lambda_max() ) ) );
        nus.erase( lower_bound( nus.begin(), nus.end(), lambda2nu( this->get_lambda_min() ) ), nus.end() );
        /* Ensure the smallest frequency interval is greater than DELTA_NU_MIN */
        impose_min_interval( nus, DELTA_NU_MIN );
#       if DEBUG_RAD > 0
        cout << "spectral_distribution_for_gas_state()" << endl
             << "nus.size() = " << nus.size() << endl;
#       endif       // DEBUG
    }
    else {
        /* Uniformally distributed spectral points with constant frequency spacing */
        double nu = lambda2nu( this->get_lambda_max() );
        double dnu = ( lambda2nu( this->get_lambda_min() ) - lambda2nu( this->get_lambda_max() ) )
                    / double ( ( this->get_spectral_points() - 1 ) );
        nus.resize( this->get_spectral_points() );
        for( int inu=0; inu<this->get_spectral_points(); ++inu ){
            nus[inu] = nu;
            nu+=dnu;
        }
    }
    return;
}

void
Photaura::
write_line_widths( Gas_data &Q )
{
    // 0. Prepare two outfiles, one for atoms and one for diatoms
    ofstream atom_file;
    atom_file.open("atomic_linewidths.txt",ios::out);
    atom_file << setprecision(12) << scientific << showpoint;
    atom_file << "# atomic_linewidths.txt" << endl
    	      << "# Column 1: Wavelength (nm)" << endl
	      << "# Column 2: Doppler width (nm)" << endl
	      << "# Column 3: Van der Waals width (nm)" << endl
	      << "# Column 4: Resonance width (nm)" << endl
	      << "# Column 5: Stark width (nm)" << endl
	      << "# Column 6: Natural width (nm)" << endl;
    
    ofstream diatom_file;
    diatom_file.open("diatomic_linewidths.txt",ios::out);
    diatom_file << setprecision(12) << scientific << showpoint;
    diatom_file << "# diatomic_linewidths.txt" << endl
    	        << "# Column 1: Wavelength (nm)" << endl
	        << "# Column 2: Doppler width (nm)" << endl
	        << "# Column 3: Collision width (nm)" << endl
	        << "# Column 4: Stark width (nm)" << endl
	        << "# Column 5: Natural width (nm)" << endl;

    // 1. Get line width string from each radiator and patch on to appropriate file
    for ( int irad=0; irad<nrad; irad++ ) {
	if ( radiators[irad]->type=="atomic_radiator" )
	    atom_file << radiators[irad]->get_line_width_string(Q);
	else if ( radiators[irad]->type=="diatomic_radiator" )
	    diatom_file << radiators[irad]->get_line_width_string(Q);
    }
    
    // 2. Close the output file
    atom_file.close();
    diatom_file.close();

    // End of void function
}

void
Photaura::
initialise_all_radiators( Gas_data &Q )
{
    // 1. Calculate partition functions
    for( int irad=0; irad<nrad; ++irad ) {
#       if DEBUG_RAD > 0
	cout << "\nCalculating partition functions for radiator: " << radiators[irad]->name << endl;
#       endif
        // Impose lower concentation limit (only if not an ion in bound-free process)
	if ( radiators[irad]->bf_ion_flag || get_rad_conc(Q,irad) > MIN_CONC ) {
	    // 0. Calculate (and store) partition functions
	    radiators[irad]->calc_partition_functions(Q);
	}
#       if DEBUG_RAD > 0
        cout << "Q_el = " << radiators[irad]->Q_el << ", Q_int = " << radiators[irad]->Q_int << endl;
#       endif
    }
    
    // Calculate electronic populations and initialise radiation mechanisms
    for( int irad=0; irad<nrad; ++irad ) {
#       if DEBUG_RAD > 0
	cout << "\nInitialising radiator: " << radiators[irad]->name << endl;
#       endif
	// Impose lower concentation limit
	if ( get_rad_conc(Q,irad) > MIN_CONC ) {
	    // 1. Calculate (and store) electronic level number densities
	    radiators[irad]->calc_elec_pops(Q);
	    // Now initialise the pseduo-constants for each radiating mechanism.
	    radiators[irad]->init_mechanisms(Q);
	}
	else {
	    if ( DEBUG_RAD > 0 ) 
		cout << "- Skipping radiator due to insufficient species population." << endl;
	}
    }
    
    // Evaluate coupled QSS system if one exists
    // if ( cqns_ ) cqns_->solve_system(Q);
}

double
Photaura::
get_rad_conc( Gas_data &Q, int irad )
{
    int isp = radiators[irad]->isp;
    return Q.massf[isp] * Q.rho / radiators[irad]->m_w;
}

void
Photaura::
prep_rad_pop_files()
{
    for ( size_t irad=0; irad<radiators.size(); ++irad ) {
    	if ( radiators[irad]->name=="e_minus" ) continue;
    	radiators[irad]->prep_spatial_population_file();
    }
    
    return;
}

void
Photaura::
append_current_rad_pops( double x )
{
    for ( size_t irad=0; irad<radiators.size(); ++irad ) {
    	if ( radiators[irad]->name=="e_minus" ) continue;
    	radiators[irad]->append_current_populations(x);
    }
    
    return;
}

void
Photaura::
write_QSS_analysis_files( Gas_data &Q, int index )
{
    for ( size_t irad=0; irad<radiators.size(); ++irad ) {
    	if ( radiators[irad]->EPM!="QSS" ) continue;
    	radiators[irad]->write_level_populations_to_file(Q,index);
    }
    
    return;
}

void impose_min_interval( std::vector<double> &V, double delta_min )
{
    size_t iV = 1, Vsize = V.size();
    double val_prev = V[0];
    while( iV < Vsize ) {
        if ( ( V[iV] - val_prev ) < delta_min ) {
            V.erase(V.begin()+iV); --Vsize;
        }
        else {
            val_prev = V[iV]; iV++;
        }
    }
    
    return;
}

