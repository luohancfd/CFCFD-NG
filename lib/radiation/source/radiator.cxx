/** \file radiator.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 21-Aug-07: Initial implementation
 *           06-Jul-09: Improved port from lib/radiation
 *
 *  \brief Definitions for class describing a generic radiating atom/molecule/electron
 *
 **/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>

#include "../../util/source/lua_service.hh"

#include "radiator.hh"
#include "atomic_radiator.hh"
#include "diatomic_radiator.hh"
#include "polyatomic_radiator.hh"
#include "electron_radiator.hh"
#include "planck_radiator.hh"
#include "radiation_constants.hh"
#include "photaura.hh"

using namespace std;

ElecLev::ElecLev( int i_val, double E_val, int g_val )
 : i(i_val), E(E_val), g(g_val) {}
 
ElecLev::~ElecLev()
{
    delete PICS_model;
}

string
ElecLev::
string()
{
    ostringstream ost;
    ost << setw(15) << E / ( RC_c * RC_h_SI )    
        << setw(5)  << g;
	
    return ost.str();
}

double
ElecLev::
calculate_equilibrium_Q_total( double T )
{
    return calculate_Q_el(T);
}

double
ElecLev::
calculate_Q_el( double T )
{
    return double(g) * exp( - E / RC_k_SI / T );
}

double
ElecLev::
calculate_and_store_Q_el( double T )
{
    Q_el = double(g) * exp( - E / RC_k_SI / T );
    
    return Q_el;
}

double
ElecLev::
calculate_boltz_N(double T, int g_i, double E_i)
{
    // NOTE: Here we are assuming sub-groupings are always in Boltzmann equilibrium

    double N_i = N * ( double(g_i) / double(g) ) * exp( - ( E_i - E ) / ( RC_k_SI * T ) );
    
    return N_i;
}


/************************** NoneqAtomicElecLev ***************************/

NoneqElecLev::NoneqElecLev( int ilev, int ne_ilev, string label, ElecLev * elev, vector<ElecLev*> eq_elevs )
: ilev( ilev ), ne_ilev( ne_ilev ), label( label ), elev( elev ), eq_elevs( eq_elevs )
{}

NoneqElecLev::~NoneqElecLev() {}

/************************** Radiator ***************************/

Radiator::Radiator( lua_State * L, string name )
: name( name )
{
    // Basic radiator data
    type = get_string( L, -1, "type" );
    
    EPM = get_string( L, -1, "E_pop_method" );
    
    isp = get_int( L, -1, "isp" );
        
    m_w = get_positive_number( L, -1, "mol_weight" );
    
    h_f = get_number( L, -1, "h_f" );
    h_f *= m_w / RC_Na;			// Convert J/kg -> J/particle
    
    I = get_number( L, -1, "eta_I" );
    I *=  RC_h_SI * RC_c;			// Convert cm^-1 to J
    
    Z = get_int( L, -1, "Z" );
    
    iT = get_int( L, -1, "iT" );
    
    iTe = get_int( L, -1, "iTe" );
    
    // Initialise the bf_ion_flag to false
    bf_ion_flag = false;
    
    if ( ECHO_RAD_INPUT > 1 ) {
	cout << "isp = " << isp << endl
	     << "mol_weight = " << m_w << " kg/mol" << endl
	     << "h_f = " << h_f << " J/kg" << endl
	     << "eta_ionise = " << I << " cm-1" << endl
	     << "Z = " << Z << endl
	     << "iTe = " << iTe << endl;
    }
}

Radiator::
~Radiator() {}

void
Radiator::set_e_index( int iel )
{
    e_index = iel;
}

double
Radiator::eval_translational_partition_function(Gas_data &Q)
{
    // Ref: Vincenti and Kruger 
    return pow( 2.0 * M_PI * m_w / RC_Na * RC_k_SI * Q.T[iT] / RC_h_SI / RC_h_SI, 1.5 );
}

double
Radiator::eval_translational_partition_function_from_T( double T )
{
    // Ref: Vincenti and Kruger 
    return pow( 2.0 * M_PI * m_w / RC_Na * RC_k_SI * T / RC_h_SI / RC_h_SI, 1.5 );
}

void
Radiator::read_photoionization_data( lua_State * L )
{
    // Create the photoionization model
    lua_getfield(L, -1, "photoionXsection_model");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Radiator::read_photoionization_data()\n";
	ost << "Error locating 'photoionXsection_model' table" << endl;
	input_error(ost);
    }
    
    for ( int ilev=0; ilev<nlevs; ++ ilev ) {
    	// hydrogenic effective principal quantum number
    	double n_eff = sqrt( RC_H_ionise_J / ( I - get_elev_pointer(ilev)->E ) );
    	if ( !isfinite(n_eff) ) {
    	    get_elev_pointer(ilev)->PICS_model = new NoPICSModel();
    	    continue;
    	}
    	// use the create function to work out what model to use
    	// NOTE: using Z+1 to represent ion charge
    	get_elev_pointer(ilev)->PICS_model = create_new_PICS_model( L, ilev, n_eff, Z+1, I );
    }
    
    lua_pop(L,1);	// pop photoionXsection_model
    
    return;
}

void
Radiator::level_population_file( Gas_data &Q, int index )
{
    calculate_Q_int(Q);
    calculate_n_e(Q);
    
    ofstream ofile;
    ostringstream oss;
    oss << name << "_level_populations-" << index << ".txt";
    string fname = oss.str();
    ofile.open(fname.c_str(),ios::out);
    ofile << "# " + name + "_level_populations.txt" << endl
    	  << "# Column 1: Level (nm)" << endl
	  << "# Column 2: Number density (particles/m**3)" << endl
	  << "# Column 3: Number density / degeneracy (particles/m**3)" << endl;
    for ( int ilev=0; ilev<nlevs; ++ilev ) {
    	ofile << setw(20) << ilev << setw(20) << get_elev_pointer(ilev)->get_N() << setw(20) << get_elev_pointer(ilev)->get_N() / get_elev_pointer(ilev)->get_g() << endl;
    }
    ofile.close();
    
    return;
}

void
Radiator::prep_x_pop_file()
{
    string fname = name + "-" + EPM + "_spatial_populations.txt";
    ofstream ofile;
    ofile.open(fname.c_str(),ios::out);
    ofile << setprecision(12) << scientific << showpoint;
    ofile << "# " + fname << endl
          << "# Column 1: Spatial location, x (m)" << endl;
    for ( int ilev=0; ilev<nlevs; ++ilev )
    	ofile << "# Column " << ilev + 2 << ": Level population, N (particles/m**3)" << endl;
    
    ofile.close();
}

void
Radiator::append_current_pops( double x )
{
    string fname = name + "-" + EPM + "_spatial_populations.txt";
    ofstream ofile;
    ofile.open(fname.c_str(),ios::app);
    ofile << setprecision(12) << scientific << showpoint;
    ofile << setw(20) << x;
    for ( int ilev=0; ilev<nlevs; ++ilev )
    	ofile << setw(20) << get_elev_pointer(ilev)->get_N();
    ofile << "\n";
    
    ofile.close();
}

bool
Radiator::optically_allowed_transition_test( int ilev_i, int ilev_f )
{
    cout << "Radiator::optically_allowed_transition_test()" << endl
         << "This function is only currently implemented for atomic radiators." << endl
         << "Exiting program." << endl;
    exit( BAD_INPUT_ERROR );
}

double
Radiator::sum_level_populations()
{
    double N_total = 0;

    for ( int ilev=0; ilev<nlevs; ++ilev )
        N_total += get_elev_pointer(ilev)->get_N();

    return N_total;
}

Radiator * create_new_radiator( lua_State * L, const std::string name )
{
    Radiator * new_radiator = 0;
    
    // Get E_pop_method to determine sub-type
    
    lua_getglobal(L, name.c_str());
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "create_new_radiator()\n";
	ost << "Error locating information table for radiator: " << name << endl;
	input_error(ost);
    }

    string type = get_string(L, -1, "type");
    
    string E_pop_method = "";
    if ( type!="electron_radiator" && type!="planck_radiator" ) {
    	E_pop_method = get_string(L, -1, "E_pop_method");
    }
    
    if ( ECHO_RAD_INPUT > 0 )
	cout << "-> Creating " << name << " as a " << E_pop_method << " " << type << endl;

    if ( type == "polyatomic_radiator" ) {
        int linear = get_int(L, -1, "linear");
        if ( linear ) {
            if ( E_pop_method=="boltzmann" )
                new_radiator = new BoltzLinearPolyatomicRadiator(L,name);
            else {
                ostringstream ost;
                ost << "create_new_radiator()" << endl
                    << "Error creating the radiator: " << name << endl
                    << "The population method: " << E_pop_method
                    << "is not yet implemented for polyatomics" << endl;
                input_error(ost);
            }
        }
        else {
            ostringstream ost;
            ost << "create_new_radiator()" << endl
                << "Error creating the radiator: " << name << endl
                << "nonlinear polyatomics are not implemented yet" << endl;
            input_error(ost);
        }
    }
    else if( type == "diatomic_radiator" ) {
	if ( E_pop_method=="boltzmann" )
	    new_radiator = new BoltzDiatomicRadiator(L,name);
	else if ( E_pop_method=="QSS" )
	    new_radiator = new QSSDiatomicRadiator(L,name);
        else if ( E_pop_method=="noneq" )
            new_radiator = new NoneqDiatomicRadiator(L,name);
    }
    else if( type == "atomic_radiator" ) {
	if ( E_pop_method=="boltzmann" )
	    new_radiator = new BoltzAtomicRadiator(L,name);
	else if ( E_pop_method=="QSS" )
	    new_radiator = new QSSAtomicRadiator(L,name);
	else if ( E_pop_method=="FirstOrderLTNE" )
	    new_radiator = new FirstOrderLTNEAtomicRadiator(L,name);
        else if ( E_pop_method=="noneq" )
            new_radiator = new NoneqAtomicRadiator(L,name);
    }
    else if ( type == "electron_radiator" ) {
	new_radiator = new ElectronRadiator(L,name);
    }
    else if( type == "planck_radiator" ) {
	new_radiator = new PlanckRadiator(L,name);
    }
    else {
	ostringstream ost;
	ost << "The specified radiator type: " << type << endl
	    << "is not available or no yet implemented.\n" << endl
	    << "Bailing Out!\n";
	input_error(ost);
    }
    
    lua_pop(L,1); 	// pop radiator
    
    return new_radiator;
}

Radiator * get_radiator_pointer_from_name( vector<Radiator*> radiators, std::string name )
{
    for ( size_t irad=0; irad<radiators.size(); ++irad ) {
    	if ( radiators[irad]->name == name ) return radiators[irad];
    }
    
    cout << "get_radiator_pointer_from_name()" << endl
         << "Radiator: " << name << " not found, bailing out!" << endl;
    exit( BAD_INPUT_ERROR );
}
