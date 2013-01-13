/** \file atomic_radiator.cxx
 *  \ingroup radiation2
 *
 *  \author Daniel F. Potter
 *  \version 10-Aug-2009: New radiation2 version
 *
 *  \brief Definitions for atomic radiator classes
 *
 **/

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "../../gas/models/gas_data.hh"

#include "atomic_radiator.hh"
#include "radiation_constants.hh"
#include "spectral_model.hh"

using namespace std;

AtomicElecLev::AtomicElecLev( int ilev, vector<double> lev_data )
 : ElecLev( ilev, lev_data[1]*RC_c*RC_h_SI, int(lev_data[2]) )		// E and g
{
    // Have already checked that lev_data.size() = 3 or 7
    
    n = int(lev_data[0]);
    
    if (lev_data.size()==7) {
	l = int(lev_data[3]);
	L = int(lev_data[4]);
	S = int(lev_data[5]);
	parity = int(lev_data[6]);
    }
    else {
	l = -1;
	L = -1;
	S = -1;
	parity = -1;
    }
}

AtomicElecLev::AtomicElecLev( AtomicElecLev * lev )
: ElecLev( lev->i, lev->get_E(), lev->get_g() )
{
    n = lev->get_n();
    l = lev->get_l();
    L = lev->get_L();
    S = lev->get_S();
    parity = lev->get_parity();
}

string
AtomicElecLev::
string()
{
    ostringstream ost;
    ost << setw(5) << n
        << setw(15) << E / ( RC_c * RC_h_SI )  
        << setw(5) << g
        << setw(5) << l
        << setw(5) << L
        << setw(5) << S
        << setw(10) << parity;
	
    return ost.str();
}

AtomicRadiator::
AtomicRadiator( lua_State * L, string name )
 : Radiator(L,name)
{
    read_elevel_data( L );
    
    read_line_data( L );
    
    read_photoionization_data( L );
}

AtomicRadiator::
~AtomicRadiator()
{
    for ( int ilev=0; ilev<nlevs; ++ilev )
    	delete elevs[ilev];
    
    for ( int iline=0; iline<nlines; ++iline )
    	delete lines[iline];
}

void
AtomicRadiator::
set_e_index( int iel )
{
    Radiator::set_e_index( iel );
}

ElecLev *
AtomicRadiator::
get_elev_pointer( int ie )
{
    return elevs[ie];
}

void
AtomicRadiator::
find_lines_in_spectral_range( double lambda_min, double lambda_max )
{
    for ( int iline=0; iline<nlines; ++iline ) {
    	double lambda_ul = nu2lambda( lines[iline]->nu_ul );
    	if ( lambda_ul >= lambda_min && lambda_ul <= lambda_max ) {
    	    cout << name << " line @ " << lambda_ul << " from ilev " << lines[iline]->elev_u->i << " to ilev " << lines[iline]->elev_l->i << " with A_ul = " << lines[iline]->A_ul << endl;
    	}
    }
    
    return;
}

void
AtomicRadiator::
read_elevel_data( lua_State * L )
{
    lua_getfield(L, -1, "level_data");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "AtomicRadiator::read_elevel_data()\n";
	ost << "Error locating 'level_data' table" << endl;
	input_error(ost);
    }
    
    nlevs = get_int(L, -1, "n_levels");
    elevs.resize( nlevs );

    for ( int ilev=0; ilev<nlevs; ++ilev ) {
	ostringstream lev_oss;
	lev_oss << "ilev_" << ilev;
	lua_getfield(L, -1, lev_oss.str().c_str());
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "AtomicRadiator::read_elevel_data()\n";
	    ost << "Error locating " << lev_oss.str() << " table" << endl;
	    input_error(ost);
	}
	vector<double> lev_data;
	for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i+1);
	    lev_data.push_back( luaL_checknumber(L, -1) );
	    lua_pop(L, 1 );
	}
	lua_pop(L,1);	// pop ilev
	// Check the size of the level data vector
	if ( lev_data.size()!=3 && lev_data.size()!=7 ) {
	    ostringstream oss;
	    oss << "AtomicRadiator::read_elevel_data()" << endl
	        << "Level data expected to have 3 or 7 elements." << endl;
	    input_error( oss );
	}	    
 	// Create the electronic level
	elevs[ilev] = new AtomicElecLev(ilev, lev_data);
    }
    lua_pop(L,1); 	// pop elec_levels
    
    if ( ECHO_RAD_INPUT > 1 ) {
	cout << "nlevs = " << nlevs << endl;
	cout << setw(10) << "[ilev]"
	     << setw(5)  << "[n_el]"
	     << setw(15) << "[E_el (1/cm)]"
	     << setw(5)  << "[g_el]"
	     << setw(5)  << "[l]"
	     << setw(5)  << "[L]"
	     << setw(5)  << "[S]"
	     << setw(10) << "[parity]" << endl;
	for ( int ilev=0; ilev<nlevs; ++ilev )
	    cout << "ilev_" << ilev << " = " << elevs[ilev]->string() << endl;
    }
    
    return;
}

void
AtomicRadiator::
read_line_data( lua_State * L )
{
    lua_getfield(L, -1, "line_data");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "AtomicRadiator::read_line_data()\n";
	ost << "Error locating 'line_data' table" << endl;
	input_error(ost);
    }
    
    nlines = get_int(L, -1, "n_lines");
    lines.resize( nlines );
    
    for ( int iline=0; iline<nlines; ++iline ) {
	ostringstream line_oss;
	line_oss << "iline_" << iline;
	lua_getfield(L, -1, line_oss.str().c_str());
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "AtomicRadiator::read_line_data()\n";
	    ost << "Error locating " << line_oss.str() << " table" << endl;
	    input_error(ost);
	}
	vector<double> line_data;
	for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i+1);
	    line_data.push_back( luaL_checknumber(L, -1) );
	    lua_pop(L, 1 );
	}
	lua_pop(L,1);	// pop iline
	// Check the size of the line data vector
	if ( line_data.size()!=5 && line_data.size()!=8 ) {
	    ostringstream oss;
	    oss << "AtomicRadiator::read_line_data()" << endl
	        << "Line data expected to have 5 or 8 elements." << endl;
	    input_error( oss );
	}	    
 	// Create the electronic level
	lines[iline] = new AtomicLine(line_data, m_w, I);
	
	// Ensure upper and lower state indices are present
	if ( lines[iline]->ie_l < 0 ) lines[iline]->ie_l = find_grouped_E_level( lines[iline]->E_l );
	if ( lines[iline]->ie_u < 0 ) lines[iline]->ie_u = find_grouped_E_level( lines[iline]->E_u );
	
	// Set electronic level pointers for lines
	lines[iline]->elev_l = dynamic_cast<AtomicElecLev*>(elevs[lines[iline]->ie_l]);
	lines[iline]->elev_u = dynamic_cast<AtomicElecLev*>(elevs[lines[iline]->ie_u]);
    }
    lua_pop(L,1); 	// pop line_data
    
    if ( ECHO_RAD_INPUT > 1 ) {
	cout << "line_data = " << nlines << endl;
	cout << setw(10) << "[iline]"
	     << setw(15) << "[E_l]"
	     << setw(15) << "[E_u]"
	     << setw(15) << "[nu_ul]"
	     << setw(5)  << "[g_l]"
	     << setw(5)  << "[g_u]"
	     << setw(15) << "[A_ul]"
	     << setw(10) << "[ie_l]"
	     << setw(10) << "[ie_u]"
	     << setw(10) << "[type]" << endl;
	for ( int iline=0; iline<nlines; ++iline )
	    cout << "iline_" << iline << " = " << lines[iline]->string() << endl;
    }
    
    return;
}

int
AtomicRadiator::
find_grouped_E_level( double E_i )
{
   /* For a specified (arbitrary) energy, find the grouped energy level with 
      the most similiar/nearest given energy level.  This crude, but should work. */
      
   double dE = 1.0;			// E_el values will be MUCH smaller than this
   int this_lev = -1;
      
   for ( int ilev=0; ilev<nlevs; ++ilev) {
       double dE_new = fabs(elevs[ilev]->get_E() - E_i);
       if ( dE_new < dE ) {
	   dE = dE_new;
	   this_lev = ilev;
       }
   }
   
   if ( this_lev < 0 ) {
       cout << "AtomicRadiator::find_grouped_E_level()" << endl
            << "Encountered an error searching for level with E = " << E_i << endl
            << "Bailing out!" << endl;
       exit( FAILURE );
   }
   
   return this_lev;
}

void
AtomicRadiator::
calculate_Q_int(Gas_data &Q)
{
    double T = Q.T[iTe];
    Q_int = 0.0;
    
    // NOTE: Q_int is calculated from simplified electronic level set
    for (int ie=0; ie<nlevs; ie++) {
	Q_int += elevs[ie]->calculate_and_store_Q_el(T);
	elevs[ie]->set_Q_int( elevs[ie]->get_Q_el() );
    }
    
    // For atoms, electronic mode is the only internal energy mode
    Q_el = Q_int;

    return;
}

double
AtomicRadiator::
calculate_total_equil_partition_function( double T )
{
    // 1. Translational contribution
    double Q_tr = pow( 2.0 * M_PI * m_w / RC_Na * RC_k_SI * T / RC_h_SI / RC_h_SI, 1.5 );
    
    // 2. Electronic contribution
    double Q_elec = 0.0;
    for (int ie=0; ie<nlevs; ie++) {
	Q_elec += elevs[ie]->calculate_Q_el(T);
    }
    
    // 3. Return the product of the modal contributions
    return Q_tr * Q_elec;
}

void
AtomicRadiator::
initialise_mechanisms( Gas_data &Q )
{
    // 0. Pre-calculate n_hvy, n_elecs, mw_av
    // NOTE: we are assuming Q.p includes the electron contribution, and that 
    //       Q.p_e is present and correct
    //       [if a 1T model is being used p_e will be zero!]
    double N_elecs = 0.0;
    if ( e_index >= 0 )
    	N_elecs = Q.p_e * RC_Na/(RC_R_u*Q.T[iTe]) * 1.0e-6;
    double N_hvy = (Q.p - Q.p_e ) * RC_Na/(RC_R_u*Q.T[iT]) * 1.0e-6;
    double mw_av = ( Q.rho * RC_Na ) / ( ( N_hvy + N_elecs ) * 1.0e6 ) * 1.0e3;
    
    // 1. Loop through lines, call initialise function.
    //    The following are calculated:
    //      - spectral coefficients
    //      - line-widths 
    for (int iline=0; iline < nlines; iline++) {
	lines[iline]->initialise( Q.T[iT], Q.T[iTe], Q.p, N_hvy, N_elecs, mw_av );
    }
    return;
}


double
AtomicRadiator::
calculate_unresolved_emission_coefficient( Gas_data &Q )
{
    double j_total = 0.0;

    // NOTE: Gas_data structure is not actually needed
    // 1. Loop over all atomic lines (work delegated to lines)
    for ( int iline=0; iline<nlines; ++iline ) {
    	j_total += lines[iline]->j_ul;
    }
    
    return j_total;
}

double
AtomicRadiator::
calculate_unresolved_OV_emission_coefficient( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u )
{
    double j_total = 0.0;

    // NOTE: gas_data structure is not actually needed
    // 1. Loop over all atomic lines (work delegated to lines)
    for ( int iline=0; iline<nlines; ++iline ) {
    	if ( nu2lambda(lines[iline]->nu_ul) < wavel_switch ) {
    	    j_total += Lambda_l * lines[iline]->j_ul;
    	}
    	else {
    	    j_total += Lambda_u * lines[iline]->j_ul;
    	}
    }
    
    return j_total;
}


void
AtomicRadiator::
spectral_distribution( std::vector<double> &nus )
{
    // 1. Loop over all atomic lines (work delegated to lines)
    for ( int iline=0; iline<nlines; ++iline ) {
    	lines[iline]->spectral_distribution(nus);
    }
    
    return;
}

void
AtomicRadiator::
calculate_spectrum( Gas_data &Q, CoeffSpectra &X )
{
    // NOTE: Gas_data structure is not actually needed
    // 1. Loop over all atomic lines (work delegated to lines)
    for ( int iline=0; iline<nlines; ++iline ) {
    	lines[iline]->calculate_spectrum(X);
    }
    
    return;
}

string
AtomicRadiator::
line_width_string( Gas_data &Q )
{
    // 0. Pre-calculate n_hvy, n_elecs, mw_av
    // NOTE: we are assuming Q.p includes the electron contribution
    double N_elecs = 0.0;
    if ( e_index >= 0 )
    	N_elecs = Q.p_e * RC_Na/(RC_R_u*Q.T[iTe]) * 1.0e-6;
    double N_hvy = (Q.p - Q.p_e ) * RC_Na/(RC_R_u*Q.T[iT]) * 1.0e-6;
    double mw_av = ( Q.rho * RC_Na ) / ( ( N_hvy + N_elecs ) * 1.0e6 ) * 1.0e3;
    
    // 1. Loop through lines, call line_width_string function
    string lws = "";
    for (int iline=0; iline<nlines; iline++) {
    	lws += lines[iline]->line_width_string(Q.T[iT], Q.T[iTe], Q.p, N_hvy, N_elecs, mw_av);
    }
    
    return lws;
}

bool
AtomicRadiator::
optically_allowed_transition_test( int ilev_i, int ilev_f )
{
    int type = get_atomic_transition_type( elevs[ilev_i], elevs[ilev_f] );
    
    return ( type==ALLOWED ) ? true : false;
}

/************************** BoltzAtomicRadiator **************************/

BoltzAtomicRadiator::
BoltzAtomicRadiator( lua_State * L, string name )
 : AtomicRadiator(L, name) {}
 
BoltzAtomicRadiator::
~BoltzAtomicRadiator() {}

void
BoltzAtomicRadiator::
calculate_n_e( Gas_data &Q )
{
    double N_total = Q.massf[isp] * Q.rho / m_w * RC_Na;	// convert kg/m**3 -> particles/m**3
    
    for (int ilev=0; ilev<nlevs; ++ilev) {
	elevs[ilev]->set_N( N_total * elevs[ilev]->get_Q_el() / Q_el );
#       if DEBUG_RAD > 0
	cout << "N_el[" << ilev << "] = " << elevs[ilev]->get_N() << endl;
#	endif
    }
    
    // End of void function
}

/*************************** QSSAtomicRadiator ***************************/

QSSAtomicRadiator::
QSSAtomicRadiator( lua_State * L, string name )
 : AtomicRadiator(L, name)
{
    lua_getfield(L, -1, "QSS_model");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "QSSAtomicRadiator::QSSAtomicRadiator()\n";
	ost << "Error locating 'QSS_model' table" << endl;
	input_error(ost);
    }
    
    // 0. Read-in the mininum temperature to apply the CR model
    T_lower = get_number( L, -1, "T_lower" );
    
    // 1. Create the nonequilibrium electronic levels
    lua_getfield(L, -1, "noneq_elevs");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "QSSAtomicRadiator::QSSAtomicRadiator()\n";
	ost << "Error locating 'noneq_elevs' table" << endl;
	input_error(ost);
    }
    
    vector<int> noneq_ilevs;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	noneq_ilevs.push_back( (int) luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L, 1);	// pop noneq_elevs table
    
    int inc_eq_elevs = get_int(L,-1,"inc_eq_elevs");
    
    for ( size_t ne_ilev=0; ne_ilev<noneq_ilevs.size(); ++ne_ilev ) {
    	int ilev = noneq_ilevs[ne_ilev];
    	// a. create vector of pointers to equilibriated levels
    	vector<ElecLev*> eq_elevs;
    	if ( inc_eq_elevs ) {
	    int next_ilev = (int) elevs.size();				// if this is the last noneq level
	    if ( ne_ilev!=noneq_ilevs.size()-1 ) next_ilev = noneq_ilevs[ne_ilev+1];	// if not
	    for ( int eq_ilev=ilev+1; eq_ilev<next_ilev; ++eq_ilev ) 
		eq_elevs.push_back( elevs[eq_ilev] );
	}
    	// b. create the label
    	ostringstream noneq_elev_label;
    	noneq_elev_label << name << "_e" << ilev;
    	// c. create the noneq level
    	noneq_elevs.push_back( new NoneqElecLev( ilev, ne_ilev, noneq_elev_label.str(), elevs[ilev], eq_elevs ) );
    }
    
    // 2. Create the collision-radiative reaction mechanisms
    // 2a. Electron impact excitation
    string eie_model = get_string(L,-1,"electron_impact_excitation");
    if ( eie_model!="none" ) {
    	int n_eie = create_electron_impact_excitation_reactions( L, eie_model );
    	cout << " - Created " << n_eie 
    	     << " electron impact excitation mechanisms for QSSAtomicRadiator: " 
    	     << name << endl;
    }
    // 2b. Electron impact ionization
    string eii_model = get_string(L,-1,"electron_impact_ionization");
    if ( eii_model!="none" ) {
    	int n_eii = create_electron_impact_ionization_reactions( L, eii_model );
    	cout << " - Created " << n_eii 
    	     << " electron impact ionization mechanisms for QSSAtomicRadiator: " 
    	     << name << endl;
    }
    // 2c. Radiative transitions
    string rt_model = get_string(L,-1,"radiative_transitions");
    if ( rt_model!="none" ) {
    	int n_rt = create_radiative_transition_reactions( L, rt_model );
    	cout << " - Created " << n_rt
    	     << " radiative transition mechanisms for QSSAtomicRadiator: " 
    	     << name << endl;
    }
    
    lua_pop(L, 1);	// pop QSS_model table
    
    // 3. Initialise the working matrices and valarrays
    dGdy = new Valmatrix();
    dGdy->resize( noneq_elevs.size(), noneq_elevs.size() );
    C.resize( noneq_elevs.size(), 0.0 );
    y_out.resize( noneq_elevs.size(), 0.0 );
}
 
QSSAtomicRadiator::
~QSSAtomicRadiator()
{
    for ( size_t i=0; i<noneq_elevs.size(); ++i )
    	delete noneq_elevs[i];
    
    for ( size_t i=0; i<reactions.size(); ++i )
    	delete reactions[i];
    
    delete dGdy;
}

void
QSSAtomicRadiator::
set_radiator_pointers( std::vector<Radiator*> radiators )
{
    ion = get_radiator_pointer_from_name( radiators, name + "_plus" );
    elec = get_radiator_pointer_from_name( radiators, "e_minus" );
    
    for ( size_t ir=0; ir<reactions.size(); ++ir ) {
    	if ( reactions[ir]->get_type()=="ElectronImpactExcitation" ) {
    	    dynamic_cast<ElectronImpactExcitation*>(reactions[ir])->set_electron_pointer( elec );
    	}
    	else if ( reactions[ir]->get_type()=="ElectronImpactIonization" ) {
    	    dynamic_cast<ElectronImpactIonization*>(reactions[ir])->set_ion_pointer( ion );
    	    dynamic_cast<ElectronImpactIonization*>(reactions[ir])->set_electron_pointer( elec );
    	}
    }
    
    return;
}

void
QSSAtomicRadiator::
level_population_file( Gas_data &Q, int index )
{
    elec->calc_partition_functions(Q);
    ion->calc_partition_functions(Q);
    calculate_Q_int(Q);
    calculate_n_e(Q);
    
    ofstream ofile;
    ostringstream oss;
    oss << name << "_level_populations-" << index << ".txt";
    string fname = oss.str();
    ofile.open(fname.c_str(),ios::out);
    ofile << setprecision(12) << scientific << showpoint;
    ofile << name + "_level_populations.txt" << endl
    	  << "# Column 1: Level" << endl
    	  << "# Column 2: Level energy (eV)" << endl
	  << "# Column 3: Number density (particles/m**3)" << endl
	  << "# Column 4: Number density / degeneracy (particles/m**3)" << endl
	  << "# Column 5: Boltzmann population / degeneracy" << endl
	  << "# Column 6: Saha population / degeneracy" << endl;
double delta_N_saha_N_boltz = 0.0;
double delta_N_QSS_N_boltz = 0.0;
    for ( int ilev=0; ilev<nlevs; ++ilev ) {
    	ofile << setw(20) << ilev
    	      << setw(20) << get_elev_pointer(ilev)->get_E() / RC_e_SI 
    	      << setw(20) << get_elev_pointer(ilev)->get_N()
    	      << setw(20) << get_elev_pointer(ilev)->get_N() / get_elev_pointer(ilev)->get_g()
    	      << setw(20) << eval_Boltzmann_population_for_level(Q,ilev) / get_elev_pointer(ilev)->get_g()
    	      << setw(20) << eval_Saha_population_for_level(Q,ilev) / get_elev_pointer(ilev)->get_g()
    	      << endl;
delta_N_saha_N_boltz += eval_Saha_population_for_level(Q,ilev) / eval_Boltzmann_population_for_level(Q,ilev) - 1.0;
delta_N_QSS_N_boltz += get_elev_pointer(ilev)->get_N() / eval_Boltzmann_population_for_level(Q,ilev) - 1.0;
    }
    ofile.close();
    
cout << name << ": delta N_saha / N_boltz = " << delta_N_saha_N_boltz / double(nlevs) << ", delta N_QSS / N_boltz = " << delta_N_QSS_N_boltz / double(nlevs) << endl;
    return;
}

double
QSSAtomicRadiator::
eval_Boltzmann_population_for_level( Gas_data &Q, int ilev )
{
    double N_total = Q.massf[isp] * Q.rho / m_w * RC_Na;	// convert kg/m**3 -> particles/m**3
    
    return N_total * elevs[ilev]->get_Q_el() / Q_el;
}

double
QSSAtomicRadiator::
eval_Saha_population_for_level( Gas_data &Q, int ilev )
{
    // 1. Number densities
    double N_elecs = Q.massf[e_index] * Q.rho / elec->m_w * RC_Na;		// convert kg/m**3 -> particles/m**3
    double N_ions = Q.massf[ion->isp] * Q.rho / ion->m_w * RC_Na;
    
    // 2. Partition functions
    double Q_level = this->eval_translational_partition_function(Q) * elevs[ilev]->get_Q_el() * exp( I / RC_k_SI / Q.T[ion->iTe] );
    double Q_elec = elec->eval_translational_partition_function(Q) * elec->Q_int;
    double Q_ion = ion->eval_translational_partition_function(Q) * ion->Q_int;
    
    // 3. Evaluate the Saha equation
    return N_ions * N_elecs * Q_level / Q_ion / Q_elec;
}

int
QSSAtomicRadiator::
create_electron_impact_excitation_reactions( lua_State * L, string model )
{
    // Loop over all possible excitation transitions
    int n_reactions = 0;
    string _model = model;
    
    // Option to ommit this mechanism
    if ( model=="none" ) return n_reactions;
    
    for ( size_t ilev=0; ilev<noneq_elevs.size(); ++ilev ) {
    	for ( size_t flev=ilev+1; flev<noneq_elevs.size(); ++flev ) {
    	    // Check that the two levels have different energies
    	    if ( noneq_elevs[ilev]->elev->E==noneq_elevs[flev]->elev->E )
    	        continue;
    	    // Can only use Frost data for certain transitions
    	    if ( model=="Frost and Drawin" ) {
    	    	if ( ilev <= 2 && flev <= 20 ) _model = "Frost";
    	    	else _model = "Drawin";
    	    }
    	    else if ( model=="Frost and Gryzinski" ) {
    	    	if ( ilev <= 2 && flev <= 20 ) _model = "Frost";
    	    	else _model = "Gryzinski";
    	    }
    	    // Can only use Bultel and Zatsarinny and Tayal data for certain transitions
    	    if ( model=="ZatsarinnyTayal and Drawin" ) {
    	    	if ( ilev==0 && ( flev<=17 || flev==19 || flev==20 || flev==21 ) ) _model = "ZatsarinnyTayal";
    	    	else if ( ilev==1 && ( flev== 2 || flev== 4 || flev== 6 || flev==8  || 
    	    	    		       flev==10 || flev==12 || flev==13 || flev==15 ||
    	    	    		       flev==17 || flev==19 || flev==21 ) )
    	    	    _model = "ZatsarinnyTayal";
    	    	else if ( ilev==2 && ( flev== 6 || flev==10 || flev==12 || flev==13 || 
    	    	    		       flev==17 || flev==19 || flev==21 ) )
    	    	    _model = "ZatsarinnyTayal";
    	    	else
    	    	    _model = "Drawin";
    	    }
    	    else if ( model=="ZatsarinnyTayal and Gryzinski" ) {
    	    	if ( ilev==0 && ( flev<=17 || flev==19 || flev==20 || flev==21 ) ) _model = "ZatsarinnyTayal";
    	    	else if ( ilev==1 && ( flev== 2 || flev== 4 || flev== 6 || flev==8  || 
    	    	    		       flev==10 || flev==12 || flev==13 || flev==15 ||
    	    	    		       flev==17 || flev==19 || flev==21 ) )
    	    	    _model = "ZatsarinnyTayal";
    	    	else if ( ilev==2 && ( flev== 6 || flev==10 || flev==12 || flev==13 || 
    	    	    		       flev==17 || flev==19 || flev==21 ) )
    	    	    _model = "ZatsarinnyTayal";
    	    	else
    	    	    _model = "Gryzinski";
    	    }
    	    // Can only use Suno and Kato data for certain transitions
    	    if ( model=="SunoKato and Drawin" ) {
    	    	if ( ilev==0 && ( flev <= 14 || flev==16 || flev==17 || flev==19 || flev==21 ) ) _model = "SunoKato";
    	    	else if ( ilev==1 && ( flev==2 || flev==4 || flev==5 || flev==6 ) ) _model = "SunoKato";
    	    	else _model = "Drawin";
    	    }
    	    else if ( model=="SunoKato and Gryzinski" ) {
    	    	if ( ilev==0 && ( flev <= 14 || flev==16 || flev==17 || flev==19 || flev==21 ) ) _model = "SunoKato";
    	    	else if ( ilev==1 && ( flev==2 || flev==4 || flev==5 || flev==6 ) ) _model = "SunoKato";
    	    	else _model = "Gryzinski";
    	    }
    	    reactions.push_back( new ElectronImpactExcitation( L, _model, this, noneq_elevs[ilev], noneq_elevs[flev] ) );
    	    n_reactions++;
    	}
    }
    
    return n_reactions;
}

int
QSSAtomicRadiator::
create_electron_impact_ionization_reactions( lua_State * L, string model )
{
    // Loop over all possible ionization transition lower states
    int n_reactions = 0;
    string _model = model;
    
    // Option to ommit this mechanism
    if ( model=="none" ) return n_reactions;
    
    for ( size_t ilev=0; ilev<noneq_elevs.size(); ++ilev ) {
    	// N and O preference the curve fits of Kunc and Soon
    	if ( model=="KuncSoon and CJDrawin" ) {
    	    if ( ilev>=0 && ilev<=2 ) _model = "KuncSoon";
    	    else _model = "CJDrawin";
    	}
    	else if ( model=="KuncSoon and Drawin" ) {
    	    if ( ilev>=0 && ilev<=2 ) _model = "KuncSoon";
    	    else _model = "Drawin";
    	}
    	// Panesi preferences the GA fits of Bultel for N and O
    	if ( model=="Bultel and CJDrawin" ) {
    	    if ( ilev>=0 && ilev<=2 ) _model = "Bultel";
    	    else _model = "CJDrawin";
    	}
    	else if ( model=="Bultel and Drawin" ) {
    	    if ( ilev>=0 && ilev<=2 ) _model = "Bultel";
    	    else _model = "Drawin";
    	}
    	// C preferences the data of Suno and Kato
    	if ( model=="SunoKato and CJDrawin" ) {
    	    if ( ilev==0 ) _model = "SunoKato";
    	    else  _model = "CJDrawin";
    	}
    	else if ( model=="SunoKato and Drawin" ) {
    	    if ( ilev==0 ) _model = "SunoKato";
    	    else  _model = "Drawin";
    	}
    	
    	reactions.push_back( new ElectronImpactIonization( L, _model, this, noneq_elevs[ilev] ) );
    	n_reactions++;
    }
    
    return n_reactions;
}

int
QSSAtomicRadiator::
create_radiative_transition_reactions( lua_State * L, string model )
{
    // Loop over all possible radiative transitions between the noneq levels
    int n_reactions = 0;
    
    // Option to ommit this mechanism
    if ( model=="none" ) return n_reactions;
    
    for ( size_t ne_ilev_u=1; ne_ilev_u<noneq_elevs.size(); ++ne_ilev_u ) {
    	for ( size_t ne_ilev_l=0; ne_ilev_l<ne_ilev_u; ++ne_ilev_l ) {
	    double gA_sum = 0.0; int g_sum = 0;
	    int ilev_u = noneq_elevs[ne_ilev_u]->ilev;
	    int ilev_l = noneq_elevs[ne_ilev_l]->ilev;
	    // Loop over all lines, find those that correspond to this transition
	    for (size_t iline=0; iline<lines.size(); ++iline) {
		if ( lines[iline]->ie_u==ilev_u && lines[iline]->ie_l==ilev_l ) {
		    gA_sum += double(lines[iline]->g_u)*lines[iline]->A_ul;
		    g_sum += lines[iline]->g_u;
		}
	    }
	    // If degeneracy is non-zero and transition probability is
	    // sufficiently high, create a new radiative transition
	    double A_ul = gA_sum / double(g_sum);
	    if ( g_sum > 0 && A_ul > MIN_TRANSITION_PROBABILITY ) {
		reactions.push_back( new RadiativeTransition( model, noneq_elevs[ne_ilev_u], noneq_elevs[ne_ilev_l], A_ul ) );
		n_reactions++;
	    }
    	}
    }
    
    return n_reactions;
}

void
QSSAtomicRadiator::
calculate_n_e( Gas_data &Q )
{
    /* Linear approach, solve system via Gaussian Elimination */
    
    // Use Boltzmann populations if temperature is low
    if ( Q.T.back() < T_lower ) {
    	for ( int ilev=0; ilev<nlevs; ++ilev ) {
	    elevs[ilev]->set_N( eval_Boltzmann_population_for_level(Q,ilev) );
    	}
    	return;
    }
    
    // 0. Reset Jacobian matrix and source and solution vectors
    for ( size_t i=0; i<noneq_elevs.size(); ++i ) {
	C[i] = 0.0;
	y_out[i] = 0.0;
	for ( size_t j=0; j<noneq_elevs.size(); ++j ) {
	    dGdy->set( i, j, 0.0 );
	}
    }
    
    // 1.  Population summations, first matrix row
    for ( size_t ne_ilev=0; ne_ilev<noneq_elevs.size(); ++ne_ilev ) {
    	double bf_acc = 0.0;	// accumulated boltzmann fractions
    	for ( size_t eq_ilev=0; eq_ilev<noneq_elevs[ne_ilev]->eq_elevs.size(); ++eq_ilev ) {
    	    bf_acc += noneq_elevs[ne_ilev]->eq_elevs[eq_ilev]->get_Q_el() / noneq_elevs[ne_ilev]->elev->get_Q_el();
    	}
    	double tmp = dGdy->get(0,ne_ilev) + 1.0 + bf_acc;
    	dGdy->set(0,ne_ilev,tmp);
    }
    
    // 2. Contributions from reactions
    for ( size_t ir=0; ir<reactions.size(); ++ir )
	reactions[ir]->add_jacobian_contributions( Q, *dGdy );
    
    // 3. Construct source vector
    C[0] = Q.massf[isp] * Q.rho / m_w / 1.0e6;		// Convert moles/m**3 to moles/cm**3
    // NOTE: some reactions such as EII have terms that will not be a function of the unkown
    //       populations, therefore they need to go in the source vector for this method
    for ( size_t ir=0; ir<reactions.size(); ++ir )
    	reactions[ir]->add_source_vector_contributions( Q, C );
    
    // 4. Solve the system
    if( dGdy->gaussian_elimination( y_out, C ) ) {
        cout << "QSSAtomicRadiator::calculate_n_e()" << endl
             << "Gaussian elimination failed for QSSAtomicRadiator: " << name << endl
             << "The gas-state was: " << endl
             << "Q.T[0] = " << Q.T[0] << endl
             << "Q.p = " << Q.p << endl
             << "Q.p_e = " << Q.p_e << endl;
        exit( FAILURE );
    }
    
    // 5.  Map results back onto radiator
    for ( size_t ne_ilev=0; ne_ilev<noneq_elevs.size(); ++ne_ilev ) {
    	// 5a. Firstly noneq levels
    	double N_ne_ilev = y_out[ne_ilev] * 1.0e6 * RC_Na;	// Convert moles/cm**3 -> particles/m**3
    	if ( N_ne_ilev < 0.0 ) {
    	    cout << "QSSDAtomicRadiator::calculate_n_e()" << endl
    	         << name << ".N[" << noneq_elevs[ne_ilev]->elev->i << "] = " << N_ne_ilev << endl
    	         << "Bailing out!" << endl;
    	    exit( FAILURE );
    	}
    	noneq_elevs[ne_ilev]->elev->set_N( N_ne_ilev );
	// 5b. Now equilibriated levels
	for ( size_t eq_ilev=0; eq_ilev<noneq_elevs[ne_ilev]->eq_elevs.size(); ++eq_ilev ) {
	    double bf = noneq_elevs[ne_ilev]->eq_elevs[eq_ilev]->get_Q_el() / noneq_elevs[ne_ilev]->elev->get_Q_el();
	    noneq_elevs[ne_ilev]->eq_elevs[eq_ilev]->set_N( N_ne_ilev * bf );
	}
    }
    
#   if ENFORCE_ATOMIC_BOLTZ_QSS_LIMIT
    // 6. Enforce the Boltzmann limit
    for ( int ilev=0; ilev<nlevs; ++ilev ) {
    	double N_boltz = eval_Boltzmann_population_for_level(Q,ilev);
    	// Boltzmann as the upper limit
    	if ( !finite(elevs[ilev]->get_N()) || elevs[ilev]->get_N() > N_boltz )
    	    elevs[ilev]->set_N(N_boltz);
    	// Check for finiteness
    	/*
    	if ( !finite( elevs[ilev]->get_N() ) ) {
    	    cout << "QSSAtomicRadiator::calculate_n_e()" << endl
    	         << "Radiator: " << name << " QSS calculation failed." << endl;
    	    cout << "Source vector: \n"; print_valarray(C);
    	    cout << "y_out: \n"; print_valarray(y_out);
    	    cout << "Bailing out!" << endl;
    	    exit( FAILURE );
    	}
    	*/
    }
#   endif
    
#   if DEBUG_RAD>0
    // 7. Check that the total species density has been conserved
    double N_total = 0.0;
    for ( int ilev=0; ilev<nlevs; ++ilev ) {
	N_total += elevs[ilev]->get_N();
	cout << "N_el[" << ilev << "] = " << elevs[ilev]->get_N() << endl;
    }
    cout << "N_total(QSS) = " << N_total << ", N_total(CFD) = " << Q.massf[isp] * Q.rho / m_w * RC_Na << endl;
    
    // 8. Check that the G vector is zero as required
    valarray<double> G;
    G.resize( noneq_elevs.size(), 0.0 );
    for ( size_t ir=0; ir<reactions.size(); ++ir )
    	reactions[ir]->add_eval_contributions( Q, G );
    for ( size_t i=0; i<G.size(); ++i )
    	cout << "G[" << i << "] = " << G[i] << endl;
#   endif
}

/********************* FirstOrderLTNEAtomicRadiator *********************/

FirstOrderLTNEAtomicRadiator::
FirstOrderLTNEAtomicRadiator( lua_State * L, string name )
 : AtomicRadiator(L, name)
{
    // Set the number of Boltzmann levels
    nlevs_boltz = 3;
    
    // Ensure there are more electronic levels than this
    if ( nlevs < 3 ) {
    	cout << "FirstOrderLTNEAtomicRadiator::FirstOrderLTNEAtomicRadiator()" << endl
    	     << "nlevs (" << nlevs << ") is less than nlevs_boltz (" << nlevs_boltz << ")" << endl
    	     << "Exiting program." << endl;
    	exit( BAD_INPUT_ERROR );
    }
}
 
FirstOrderLTNEAtomicRadiator::
~FirstOrderLTNEAtomicRadiator() {}

void
FirstOrderLTNEAtomicRadiator::
set_radiator_pointers( std::vector<Radiator*> radiators )
{
    ion = get_radiator_pointer_from_name( radiators, name + "_plus" );
    elec = get_radiator_pointer_from_name( radiators, "e_minus" );
    
    return;
}

void
FirstOrderLTNEAtomicRadiator::
calculate_n_e( Gas_data &Q )
{
    // 1. Calculate the Boltzmann populations    
    double N_total = Q.massf[isp] * Q.rho / m_w * RC_Na;	// convert kg/m**3 -> particles/m**3
    
    for (int ilev=0; ilev<nlevs_boltz; ++ilev) {
	elevs[ilev]->set_N( N_total * elevs[ilev]->get_Q_el() / Q_el );
    }
    
    // 2. Calculate the Saha populations
    // 2a. Number densities
    double N_elecs = Q.massf[e_index] * Q.rho / RC_m_SI;
    double N_ions = Q.massf[ion->isp] * Q.rho / ion->m_w * RC_Na;
    
    // 2b. Parition functions
    double Q_ion = ion->eval_translational_partition_function(Q) * ion->Q_int * exp( - I / RC_k_SI / Q.T[ion->iTe] );
    double Q_elec = elec->eval_translational_partition_function(Q) * elec->Q_int;
    double Q_trans = this->eval_translational_partition_function(Q);
    
    // 2c. Loop over remaining levels and apply the Saha equation
    for ( int ilev=nlevs_boltz; ilev<nlevs; ++ilev ) {
    	elevs[ilev]->set_N( N_ions * N_elecs * elevs[ilev]->get_Q_el() * Q_trans / Q_ion / Q_elec );
    }
    
    // End of void function
}

/************************** NoneqAtomicRadiator **************************/

NoneqAtomicRadiator::
NoneqAtomicRadiator( lua_State * L, string name )
 : AtomicRadiator(L, name)
{
    lua_getfield(L, -1, "level_data");
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "NoneqAtomicRadiator::NoneqAtomicRadiator()\n";
        ost << "Error locating 'level_data' table" << endl;
        input_error(ost);
    }

    lua_getfield(L, -1, "isp_list");
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "NoneqAtomicRadiator::NoneqAtomicRadiator()\n";
        ost << "Error locating isp_list table" << endl;
        input_error(ost);
    }
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
        lua_rawgeti(L, -1, i+1);
        isp_list.push_back( luaL_checknumber(L, -1) );
        lua_pop(L, 1 );
    }

    lua_pop(L,1);       // pop isp_list

    lua_pop(L,1);       // pop level_data
}

NoneqAtomicRadiator::
~NoneqAtomicRadiator() {}

void
NoneqAtomicRadiator::
calculate_n_e( Gas_data &Q )
{
    const double constA = Q.rho / m_w * RC_Na;        // convert kg/m**3 -> particles/m**3

    int last_neq_lev = -1;

    for (int ilev=0; ilev<nlevs; ++ilev) {
        int lev_isp = isp_list[ilev];
        if ( lev_isp >= 0 ) {
            // directly set the level population
            elevs[ilev]->set_N( Q.massf[lev_isp] * constA );
            last_neq_lev = ilev;
        }
        else if ( last_neq_lev >= 0 ) {
            // boltzmann equilibrium with last noneq level
            const double Q_ratio = elevs[ilev]->get_Q_int() / elevs[last_neq_lev]->get_Q_int();
            elevs[ilev]->set_N( Q_ratio * elevs[last_neq_lev]->get_N() );
        }
        else {
            cout << "NoneqAtomicRadiator::calculate_n_e()" << endl
                 << "Level %d does not have a nonequilibrium level below it." << endl;
            exit( FAILURE );
        }

#       if DEBUG_RAD > 0
        cout << "N_el[" << ilev << "] = " << elevs[ilev]->get_N() << endl;
#       endif
    }
    // End of void function
}

int get_atomic_transition_type( AtomicElecLev * lev_i, AtomicElecLev * lev_f )
{
    // Check if data is provided
    if ( lev_i->get_l()==-1      || lev_f->get_l()==-1 || 
         lev_i->get_L()==-1      || lev_f->get_L()==-1 || 
         lev_i->get_S()==-1      || lev_f->get_S()==-1 || 
         lev_i->get_parity()==-1 || lev_f->get_parity()==-1 ) return NO_DATA;
    
    // Determine transition type (optically allowed or optically forbidden)
    // orbital quantum number l: +/- 1
    if ( abs(lev_i->get_l()-lev_f->get_l())!=1 ) return FORBIDDEN;
    // OAM quantum number L: +/-1,0 but not 0-0
    if ( ( abs(lev_i->get_L()-lev_f->get_L())!=1 && abs(lev_i->get_L()-lev_f->get_L())!=0 ) ||
	 ( lev_i->get_L()==0 && lev_f->get_L()==0 ) ) return FORBIDDEN;
    // resultant spin quantum number S: 0
    if ( lev_i->get_S()!=lev_f->get_S() ) return FORBIDDEN;
    // parity: must change
    if ( lev_i->get_parity()==lev_f->get_parity() ) return FORBIDDEN;
    
    // If we got to here the transition is optically allowed
    return ALLOWED;
}
