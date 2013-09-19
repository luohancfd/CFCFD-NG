/** \file cr_reactions.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 03-04-10: Improved port from old lib/radiation
 *  \brief Declarations for nonequilibirum population model transitions
 *
 **/
 
#include <sstream>
 
#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
 
#include "cr_reactions.hh"
#include "radiation_constants.hh"
#include "atomic_radiator.hh"
#include "diatomic_radiator.hh"
#include "spectra_pieces.hh"

using namespace std;

/******************************* CR_Reaction ********************************/

CR_Reaction::CR_Reaction( string type )
: type( type )
{}

CR_Reaction::~CR_Reaction()
{
    delete forward_rate_coeff;
    delete backward_rate_coeff;
}

string CR_Reaction::get_type()
{
    return type;
}

string CR_Reaction::get_equation()
{
    return equation;
}

string CR_Reaction::get_forward_rate_coeff_type()
{
    return forward_rate_coeff->get_type_str();
}

int
CR_Reaction::
eval_reaction_rates( double T_f, double T_b, Gas_data &Q, double &k_f, double &k_b )
{
    if ( backward_rate_coeff->get_equilibrium_flag() ) {
	k_f = forward_rate_coeff->get_rate( T_f, Q );
	double Kc = this->eval_equilibrium_constant( T_b );
	double k_f_star;
	if ( T_b != T_f ) k_f_star = forward_rate_coeff->get_rate( T_b, Q );
	else              k_f_star = k_f;
	k_b =  k_f_star / Kc;
    }
    else if ( forward_rate_coeff->get_equilibrium_flag() ) {
	k_b = backward_rate_coeff->get_rate( T_b, Q );
	double Kc = this->eval_equilibrium_constant( T_f );
	double k_b_star;
	if ( T_f != T_b ) k_b_star = backward_rate_coeff->get_rate( T_b, Q );
	else              k_b_star = k_b;
	k_f = k_b_star * Kc;
    }
    else {
	k_f = forward_rate_coeff->get_rate( T_f, Q );
	k_b = backward_rate_coeff->get_rate( T_b, Q );
    }
    
    if ( !isfinite(k_f) || !isfinite(k_b) ) {
    	cout << "CR_Reaction::eval_reaction_rates()" << endl
    	     << "k_f = " << k_f << ", k_b = " << k_b << ", T[0] = " << Q.T[0] << endl
    	     << this->get_type() << endl
    	     << this->get_equation() << endl
    	     << this->get_equation() << endl;
    	exit( FAILURE );
    }
    
    return SUCCESS;
}

/************************ HeavyParticleImpactExcitation **************************/

HeavyParticleImpactExcitation::
HeavyParticleImpactExcitation( lua_State * L, Radiator * rad )
: CR_Reaction( "HeavyParticleImpactExcitation" ), rad( rad )
{
    QSSDiatomicRadiator * QSS_rad = dynamic_cast<QSSDiatomicRadiator*>(rad);
    
    iTv = QSS_rad->iTv;
    iT = QSS_rad->iT;
    
    equation = get_string( L, -1, "equation" );
    cout << "- Creating a HeavyParticleImpactExcitation reaction with equation: " << equation << endl;
    
    vector<string> equation_tks;
    tokenize_equation_string( equation, equation_tks );
    
    M_name = equation_tks[2];
    
    ne_elev_l = 0; ne_elev_u = 0;
    
    for ( size_t i=0; i<QSS_rad->noneq_elevs.size(); ++i ) {
    	if ( QSS_rad->name + "_" + QSS_rad->noneq_elevs[i]->label == equation_tks[0] )
    	    ne_elev_l = QSS_rad->noneq_elevs[i];
    	if ( QSS_rad->name + "_" + QSS_rad->noneq_elevs[i]->label == equation_tks[4] )
    	    ne_elev_u = QSS_rad->noneq_elevs[i];
    }
    
    if ( ne_elev_l==0 ) {
    	ostringstream oss;
    	oss << "HeavyParticleImpactExcitation::HeavyParticleImpactExcitation()" << endl
    	    << "Lower level: " << equation_tks[0] << " not found." << endl;
    	input_error( oss );
    }
    if ( ne_elev_u==0 ) {
    	ostringstream oss;
    	oss << "HeavyParticleImpactExcitation::HeavyParticleImpactExcitation()" << endl
    	    << "Upper level: " << equation_tks[4] << " not found." << endl;
    	input_error( oss );
    }
    
    forward_rate_coeff = create_explicit_rate_coeff( L, "forward" );
    backward_rate_coeff = create_explicit_rate_coeff( L, "backward" );
}

double
HeavyParticleImpactExcitation::
eval_equilibrium_constant( double T )
{
    // Kc = PI(Q)_products / PI(Q)_reactants
    
    double Kc = ne_elev_u->elev->calculate_equilibrium_Q_total(T) / ne_elev_l->elev->calculate_equilibrium_Q_total(T);

    return Kc;
}

int
HeavyParticleImpactExcitation::
add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using the assumed temperatures
    double T_f = sqrt( Q.T[iT] * Q.T[iTv] );
    double T_b = Q.T[iT];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of third-bodies
    double c_M = 0.0;
    if ( M ) c_M = Q.massf[M->isp] * Q.rho / M->m_w / 1.0e6;
    else c_M = ( Q.p - Q.p_e ) / ( RC_k_SI * Q.T[iT] )  / 1.0e6;
    
    // 1. Compute derivates and fill in the Jacobian matrix
    double domega_dc, tmp;
    
    // 1a. Lower state row, but not if it is the ground state
    if ( ne_elev_l->ne_ilev != 0 ) {
    	// derivative of omega wrt lower state population
    	domega_dc = - k_f * c_M;
    	tmp = dGdy.get( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev );
    	dGdy.set( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev, (tmp+domega_dc) );
    	// derivative of omega wrt upper state population
    	domega_dc = k_b * c_M;
    	tmp = dGdy.get( ne_elev_l->ne_ilev, ne_elev_u->ne_ilev );
    	dGdy.set( ne_elev_l->ne_ilev, ne_elev_u->ne_ilev, (tmp+domega_dc) );
    }
    
    // 1b. Upper state row (should never be the ground state)
    // derivative of omega wrt lower state population
    domega_dc = k_f * c_M;
    tmp = dGdy.get( ne_elev_u->ne_ilev, ne_elev_l->ne_ilev );
    dGdy.set( ne_elev_u->ne_ilev, ne_elev_l->ne_ilev, (tmp+domega_dc) );
    // derivative of omega wrt upper state population
    domega_dc = - k_b * c_M;
    tmp = dGdy.get( ne_elev_u->ne_ilev, ne_elev_u->ne_ilev );
    dGdy.set( ne_elev_u->ne_ilev, ne_elev_u->ne_ilev, (tmp+domega_dc) );
    
    return SUCCESS;
}

int
HeavyParticleImpactExcitation::
add_eval_contributions( Gas_data &Q, valarray<double> &G )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using the assumed temperatures
    double T_f = sqrt( Q.T[iT] * Q.T[iTv] );
    double T_b = Q.T[iT];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of third-bodies
    double c_M = 0.0;
    if ( M ) c_M = Q.massf[M->isp] * Q.rho / M->m_w / 1.0e6;
    else c_M = ( Q.p - Q.p_e ) / ( RC_k_SI * Q.T[iT] )  / 1.0e6;
    
    // 0c. moles per cm**3 of upper and lower states
    double c_u = ne_elev_u->elev->N / 1.0e6 / RC_Na;
    double c_l = ne_elev_l->elev->N / 1.0e6 / RC_Na;
    
    // 1. Fill in the G vector
    
    // 1a. Upper state
    G[ne_elev_u->ne_ilev] += k_f * c_l * c_M - k_b * c_u * c_M;
    
    // 1b. Lower state - if not hte ground state
    if ( ne_elev_l->ne_ilev != 0 )
    	G[ne_elev_l->ne_ilev] += k_b * c_u * c_M - k_f * c_l * c_M;
    
    return SUCCESS;
}

int
HeavyParticleImpactExcitation::
add_source_vector_contributions( Gas_data &Q, valarray<double> &C )
{
    // Nothing to be done
    return SUCCESS;
}

string 
HeavyParticleImpactExcitation::
get_latex_string()
{
    string lname, lev_prefix;
    get_latex_species_piecies( rad->name, lname, lev_prefix );
    
    string lM_name;
    if ( M_name!="M" ) {
    	string tmp;
    	get_latex_species_piecies( M_name, lM_name, tmp );
    }
    else {
    	lM_name = M_name;
    }
    
    string indent = "           ";
    
    ostringstream oss;
    oss << indent << "\\PIEreac{" << lname << "}{" << lM_name << "}{" << lev_prefix << "lev" << ne_elev_l->label << "}{" << lev_prefix << "lev" << ne_elev_u->label << "}";
    oss << "\t\t&\t\t" << forward_rate_coeff->get_latex_string() << "\t\\\\" << endl;
    
    return oss.str();
}

/************************ ElectronImpactExcitation **************************/

ElectronImpactExcitation::
ElectronImpactExcitation( lua_State * L, Radiator * rad )
: CR_Reaction( "ElectronImpactExcitation" )
{
    QSSDiatomicRadiator * QSS_rad = dynamic_cast<QSSDiatomicRadiator*>(rad);
    
    iTe = QSS_rad->iTe;
    iT = QSS_rad->iT;
    
    equation = get_string( L, -1, "equation" );
    cout << "- Creating a ElectronImpactExcitation reaction with equation: " << equation << endl;
    
    vector<string> equation_tks;
    tokenize_equation_string( equation, equation_tks );
    
    ne_elev_l = 0; ne_elev_u = 0;
    
    for ( size_t i=0; i<QSS_rad->noneq_elevs.size(); ++i ) {
    	if ( QSS_rad->name + "_" + QSS_rad->noneq_elevs[i]->label == equation_tks[0] )
    	    ne_elev_l = QSS_rad->noneq_elevs[i];
    	if ( QSS_rad->name + "_" + QSS_rad->noneq_elevs[i]->label == equation_tks[4] )
    	    ne_elev_u = QSS_rad->noneq_elevs[i];
    }
    
    if ( ne_elev_l==0 ) {
    	ostringstream oss;
    	oss << "ElectronImpactExcitation::ElectronImpactExcitation()" << endl
    	    << "Lower level: " << equation_tks[0] << " not found." << endl;
    	input_error( oss );
    }
    if ( ne_elev_u==0 ) {
    	ostringstream oss;
    	oss << "ElectronImpactExcitation::ElectronImpactExcitation()" << endl
    	    << "Upper level: " << equation_tks[4] << " not found." << endl;
    	input_error( oss );
    }
    
    forward_rate_coeff = create_explicit_rate_coeff( L, "forward" );
    backward_rate_coeff = create_explicit_rate_coeff( L, "backward" );

    rad_name = rad->name;
}

ElectronImpactExcitation::
ElectronImpactExcitation( lua_State * L, std::string model, Radiator * rad, NoneqElecLev * ne_elev_l, NoneqElecLev * ne_elev_u )
: CR_Reaction( "ElectronImpactExcitation" ), iTe( rad->iTe ), iT( rad->iT ), ne_elev_l( ne_elev_l ), ne_elev_u( ne_elev_u )
{
    // Create the equation string
    equation = ne_elev_l->label + " + e_minus <=> " + ne_elev_u->label + " + e_minus";
    
    if ( model=="Drawin" ) {
    	int transition_type = get_atomic_transition_type( dynamic_cast<AtomicElecLev*>(ne_elev_l->elev),
    	    						  dynamic_cast<AtomicElecLev*>(ne_elev_u->elev) );
    	if ( transition_type == ALLOWED ) 
    	    forward_rate_coeff = new DrawinOpticallyAllowedElectronImpactExcitation(ne_elev_l->elev->get_E(),ne_elev_u->elev->get_E());
    	else
    	    forward_rate_coeff = new DrawinOpticallyForbiddenElectronImpactExcitation(ne_elev_l->elev->get_E(),ne_elev_u->elev->get_E());
    	backward_rate_coeff = new FromEquilibriumConstant();
    }
    else if ( model=="Frost" ) {
    	// Check that this is N
    	if ( rad->name != "N" ) {
    	    ostringstream oss;
    	    oss << "ElectronImpactExcitation::ElectronImpactExcitation()" << endl
    	        << "The Frost electron-impact excitation model is meant for use with N not " << rad->name << "." << endl;
    	    input_error( oss );
    	}
    	lua_getfield(L,-1,"frost_data");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactIonization::ElectronImpactIonization()\n";
	    ost << "Error locating 'frost_data' table" << endl;
	    input_error(ost);
	}
    	forward_rate_coeff = new FrostNitrogenElectronImpactExcitation( L, ne_elev_l->elev, ne_elev_u->elev );
    	backward_rate_coeff = new FromEquilibriumConstant();
    	lua_pop(L,1); // pop frost_data
    }
    else if ( model=="Gryzinski" ) {
    	forward_rate_coeff = new GryzinskiElectronImpactExcitation( rad, ne_elev_l->elev, ne_elev_u->elev );
    	backward_rate_coeff = new FromEquilibriumConstant();
    }
    else if ( model=="ZatsarinnyTayal" ) {
    	// Check that this is O
    	if ( rad->name != "O" ) {
    	    ostringstream oss;
    	    oss << "ElectronImpactExcitation::ElectronImpactExcitation()" << endl
    	        << "The Zatsarinny-Tayal electron-impact excitation model is meant for use with O not " << rad->name << "." << endl;
    	    input_error( oss );
    	}
    	lua_getfield(L,-1,"zatsarinny_tayal_data");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactIonization::ElectronImpactIonization()\n";
	    ost << "Error locating 'zatsarinny_tayal_data' table" << endl;
	    input_error(ost);
	}
    	forward_rate_coeff = new ZatsarinnyTayalOxygenElectronImpactExcitation( L, ne_elev_l->elev, ne_elev_u->elev );
    	backward_rate_coeff = new FromEquilibriumConstant();
    	lua_pop(L,1); // pop zatsarinny_tayal_data
    }
    else if ( model=="Bultel" ) {
    	lua_getfield(L,-1,"bultel_data");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactIonization::ElectronImpactIonization()\n";
	    ost << "Error locating 'bultel_data' table" << endl;
	    input_error(ost);
	}
    	ostringstream oss;
    	oss << "reac_EIE_ilev" << ne_elev_l->elev->i << "_to_ilev" << ne_elev_u->elev->i;
    	lua_getfield(L,-1,oss.str().c_str());
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactExcitation::ElectronImpactExcitation()\n";
	    ost << "Error locating " << oss.str() << " table" << endl;
	    input_error(ost);
	}
    	forward_rate_coeff = create_explicit_rate_coeff( L, "forward" );
    	backward_rate_coeff = create_explicit_rate_coeff( L, "backward" );
    	lua_pop(L,1); // pop reaction
    	lua_pop(L,1); // pop bultel_data
    }
    else if ( model=="SunoKato" ) {
    	lua_getfield(L,-1,"suno_kato_data");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactIonization::ElectronImpactIonization()\n";
	    ost << "Error locating 'suno_kato_data' table" << endl;
	    input_error(ost);
	}
    	ostringstream oss;
    	oss << "reac_EIE_ilev" << ne_elev_l->elev->i << "_to_ilev" << ne_elev_u->elev->i;
    	lua_getfield(L,-1,oss.str().c_str());
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactExcitation::ElectronImpactExcitation()\n";
	    ost << "Error locating " << oss.str() << " table" << endl;
	    input_error(ost);
	}
    	forward_rate_coeff = new SunoKatoCarbonElectronImpactExcitation( L, ne_elev_l->elev, ne_elev_u->elev );
    	backward_rate_coeff = new FromEquilibriumConstant();
    	lua_pop(L,1); // pop reaction
    	lua_pop(L,1); // pop suno_kato_data
    }
    else {
    	ostringstream ost;
    	ost << "ElectronImpactExcitation::ElectronImpactExcitation()" << endl
    	    << "The " << model << " rate coefficient model is not implemented." << endl;
    	input_error(ost);
    }
    
    rad_name = rad->name;
}

double
ElectronImpactExcitation::
eval_equilibrium_constant( double T )
{
    // Kc = PI(Q)_products / PI(Q)_reactants
    
    double Kc = ne_elev_u->elev->calculate_equilibrium_Q_total(T) / ne_elev_l->elev->calculate_equilibrium_Q_total(T);

    return Kc;
}

int
ElectronImpactExcitation::
add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using the electronic temperature
    double T_f = Q.T[iTe];
    double T_b = Q.T[iTe];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of electrons
    double c_e = Q.massf[elec->isp] * Q.rho / elec->m_w / 1.0e6;
    
    // 1. Compute derivates and fill in the Jacobian matrix
    double domega_dc, tmp;
    
    // 1a. Lower state row, but not if it is the ground state
    if ( ne_elev_l->ne_ilev != 0 ) {
    	// derivative of omega wrt lower state population
    	domega_dc = - k_f * c_e;
    	tmp = dGdy.get( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev );
    	dGdy.set( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev, (tmp+domega_dc) );
    	// derivative of omega wrt upper state population
    	domega_dc = k_b * c_e;
    	tmp = dGdy.get( ne_elev_l->ne_ilev, ne_elev_u->ne_ilev );
    	dGdy.set( ne_elev_l->ne_ilev, ne_elev_u->ne_ilev, (tmp+domega_dc) );
    }
    
    // 1b. Upper state row (should never be the ground state)
    // derivative of omega wrt lower state population
    domega_dc = k_f * c_e;
    tmp = dGdy.get( ne_elev_u->ne_ilev, ne_elev_l->ne_ilev );
    dGdy.set( ne_elev_u->ne_ilev, ne_elev_l->ne_ilev, (tmp+domega_dc) );
    // derivative of omega wrt upper state population
    domega_dc = - k_b * c_e;
    tmp = dGdy.get( ne_elev_u->ne_ilev, ne_elev_u->ne_ilev );
    dGdy.set( ne_elev_u->ne_ilev, ne_elev_u->ne_ilev, (tmp+domega_dc) );
    
    return SUCCESS;
}

int
ElectronImpactExcitation::
add_eval_contributions( Gas_data &Q, valarray<double> &G )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using the electronic temperature
    double T_f = Q.T[iTe];
    double T_b = Q.T[iTe];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of electrons
    double c_e = Q.massf[elec->isp] * Q.rho / elec->m_w / 1.0e6;
    
    // 0c. moles per cm**3 of upper and lower states
    double c_u = ne_elev_u->elev->N / 1.0e6 / RC_Na;
    double c_l = ne_elev_l->elev->N / 1.0e6 / RC_Na;
    
    // 1. Fill in the G vector
    
    // 1a. Upper state
    G[ne_elev_u->ne_ilev] += k_f * c_l * c_e - k_b * c_u * c_e;
    
    // 1b. Lower state - if not hte ground state
    if ( ne_elev_l->ne_ilev != 0 )
    	G[ne_elev_l->ne_ilev] += k_b * c_u * c_e - k_f * c_l * c_e;
    
    return SUCCESS;
}

int
ElectronImpactExcitation::
add_source_vector_contributions( Gas_data &Q, valarray<double> &C )
{
    // Nothing to be done
    return SUCCESS;
}

string 
ElectronImpactExcitation::
get_latex_string()
{
    string lname, lev_prefix;
    get_latex_species_piecies( rad_name, lname, lev_prefix );
    
    string lM_name = "e$^-$";
    
    string indent = "           ";
    
    ostringstream oss;
    oss << indent << "\\PIEreac{" << lname << "}{" << lM_name << "}{" << lev_prefix << "lev" << ne_elev_l->label << "}{" << lev_prefix << "lev" << ne_elev_u->label << "}";
    oss << "\t\t&\t\t" << forward_rate_coeff->get_latex_string() << "\t\\\\" << endl;
    
    return oss.str();
}

/************************ HeavyParticleImpactDissociation **************************/

HeavyParticleImpactDissociation::
HeavyParticleImpactDissociation( lua_State * L, Radiator * rad )
: CR_Reaction( "HeavyParticleImpactDissociation" ), rad( rad )
{
    QSSDiatomicRadiator * QSS_rad = dynamic_cast<QSSDiatomicRadiator*>(rad);
    
    D = QSS_rad->get_D();
    iTv = QSS_rad->iTv;
    iT = QSS_rad->iT;
    
    equation = get_string( L, -1, "equation" );
    cout << "- Creating a HeavyParticleImpactDissociation reaction with equation: " << equation << endl;
    
    vector<string> equation_tks;
    tokenize_equation_string( equation, equation_tks );
    
    ne_elev_l = 0;
    
    for ( size_t i=0; i<QSS_rad->noneq_elevs.size(); ++i ) {
    	if ( QSS_rad->name + "_" + QSS_rad->noneq_elevs[i]->label == equation_tks[0] )
    	    ne_elev_l = QSS_rad->noneq_elevs[i];
    }
    
    if ( ne_elev_l==0 ) {
    	ostringstream oss;
    	oss << "HeavyParticleImpactDissociation::HeavyParticleImpactDissociation()" << endl
    	    << "Lower level: " << equation_tks[0] << " not found." << endl;
    	input_error( oss );
    }
    
    M_name = equation_tks[2];
    atom_A_name = equation_tks[4];
    atom_B_name = equation_tks[6];
    
    forward_rate_coeff = create_explicit_rate_coeff( L, "forward" );
    backward_rate_coeff = create_explicit_rate_coeff( L, "backward" );

}

double
HeavyParticleImpactDissociation::
eval_equilibrium_constant( double T )
{
    double Q_atom_A = atom_A->calc_total_equil_partition_function(T);
    double Q_atom_B = atom_B->calc_total_equil_partition_function(T);
    double Q_elev = ne_elev_l->elev->calculate_equilibrium_Q_total(T) * rad->eval_translational_partition_function_from_T(T) * exp( D / RC_k_SI / T );
    
    // Kc = PI(Q)_products / PI(Q)_reactants
    
    double Kc = Q_atom_A * Q_atom_B / Q_elev;
    
    Kc /= RC_Na * 1.0e6;		// Convert particles / m**3 -> moles / cm**3
    
    return Kc;
}

int
HeavyParticleImpactDissociation::
add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy )
{
    // 0. Prepare data
    // 0a. Forwardreaction rates using the geometric average temperature
    double T_f = sqrt( Q.T[iT] * Q.T[iTv] );
    double T_b = Q.T[iT];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of third-bodies
    double c_M = 0.0;
    if ( M ) c_M = Q.massf[M->isp] * Q.rho / M->m_w / 1.0e6;
    else c_M = ( Q.p - Q.p_e ) / ( RC_k_SI * Q.T[iT] )  / 1.0e6;
    
    // 1. Compute derivates and fill in the Jacobian matrix
    double domega_dc, tmp;
    
    // 1a. Lower state row, but not if it is the ground state
    if ( ne_elev_l->ne_ilev != 0 ) {
    	// derivative of omega wrt lower state population
    	domega_dc = - k_f * c_M;
    	tmp = dGdy.get( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev );
    	dGdy.set( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev, (tmp+domega_dc) );
    }
    
    return SUCCESS;
}

int
HeavyParticleImpactDissociation::
add_eval_contributions( Gas_data &Q, valarray<double> &G )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using the geometric average temperature
    double T_f = sqrt( Q.T[iT] * Q.T[iTv] );
    double T_b = Q.T[iT];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of third-bodies
    double c_M = 0.0;
    if ( M ) c_M = Q.massf[M->isp] * Q.rho / M->m_w / 1.0e6;
    else c_M = ( Q.p - Q.p_e ) / ( RC_k_SI * Q.T[iT] )  / 1.0e6;
    
    // 0c. moles per cm**3 of the atoms
    double c_A = Q.massf[atom_A->isp] * Q.rho / atom_A->m_w / 1.0e6;
    double c_B = Q.massf[atom_B->isp] * Q.rho / atom_B->m_w / 1.0e6;
    
    // 0d. moles per cm**3 of the lower state
    double c_l = ne_elev_l->elev->N / RC_Na / 1.0e6;
    
    // 1. Fill in the G vector
    
    // 1a. Lower state - if not hte ground state
    if ( ne_elev_l->ne_ilev != 0 )
    	G[ne_elev_l->ne_ilev] += k_b * c_A * c_B * c_M - k_f * c_l * c_M;
    
    return SUCCESS;
}

int
HeavyParticleImpactDissociation::
add_source_vector_contributions( Gas_data &Q, valarray<double> &C )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using the geometric average temperature
    double T_f = sqrt( Q.T[iT] * Q.T[iTv] );
    double T_b = Q.T[iT];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of third-bodies
    double c_M = 0.0;
    if ( M ) c_M = Q.massf[M->isp] * Q.rho / M->m_w / 1.0e6;
    else c_M = ( Q.p - Q.p_e ) / ( RC_k_SI * Q.T[iT] )  / 1.0e6;
    
    // 0c. moles per cm**3 of the atoms
    double c_A = Q.massf[atom_A->isp] * Q.rho / atom_A->m_w / 1.0e6;
    double c_B = Q.massf[atom_B->isp] * Q.rho / atom_B->m_w / 1.0e6;
    
    // 1. Fill in the C vector
    
    // 1a. Lower state - if not the ground state
    if ( ne_elev_l->ne_ilev != 0 )
    	C[ne_elev_l->ne_ilev] -= k_b * c_A * c_B * c_M;
    
    return SUCCESS;
}

string 
HeavyParticleImpactDissociation::
get_latex_string()
{
    return "";
}

/************************ ElectronImpactDissociation **************************/

ElectronImpactDissociation::
ElectronImpactDissociation( lua_State * L, Radiator * rad )
: CR_Reaction( "ElectronImpactDissociation" ), rad( rad )
{
    QSSDiatomicRadiator * QSS_rad = dynamic_cast<QSSDiatomicRadiator*>(rad);
    
    D = QSS_rad->get_D();
    iTv = QSS_rad->iTv;
    iTe = QSS_rad->iTe;
    iT = QSS_rad->iT;
    
    equation = get_string( L, -1, "equation" );
    cout << "- Creating a ElectronImpactDissociation reaction with equation: " << equation << endl;
    
    vector<string> equation_tks;
    tokenize_equation_string( equation, equation_tks );
    
    ne_elev_l = 0;
    
    for ( size_t i=0; i<QSS_rad->noneq_elevs.size(); ++i ) {
    	if ( QSS_rad->name + "_" + QSS_rad->noneq_elevs[i]->label == equation_tks[0] )
    	    ne_elev_l = QSS_rad->noneq_elevs[i];
    }
    
    if ( ne_elev_l==0 ) {
    	ostringstream oss;
    	oss << "ElectronImpactExcitation::ElectronImpactExcitation()" << endl
    	    << "Lower level: " << equation_tks[0] << " not found." << endl;
    	input_error( oss );
    }
    
    atom_A_name = equation_tks[4];
    atom_B_name = equation_tks[6];
    
    forward_rate_coeff = create_explicit_rate_coeff( L, "forward" );
    backward_rate_coeff = create_explicit_rate_coeff( L, "backward" );
    
    rad_name = rad->name;
}

double
ElectronImpactDissociation::
eval_equilibrium_constant( double T )
{
    double Q_atom_A = atom_A->calc_total_equil_partition_function(T);
    double Q_atom_B = atom_B->calc_total_equil_partition_function(T);
    double Q_elev = ne_elev_l->elev->calculate_equilibrium_Q_total(T) * rad->eval_translational_partition_function_from_T(T) * exp( D / RC_k_SI / T );
    
    // Kc = PI(Q)_products / PI(Q)_reactants
    
    double Kc = Q_atom_A * Q_atom_B / Q_elev;
    
    Kc /= RC_Na * 1.0e6;		// Convert particles / m**3 -> moles / cm**3
    
    return Kc;
}

int
ElectronImpactDissociation::
add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using the geometric average temperature
    double T_f = sqrt( Q.T[iTv] * Q.T[iTe] );
    double T_b = Q.T[iTe];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of electrons
    double c_e = Q.massf[elec->isp] * Q.rho / elec->m_w / 1.0e6;
    
    // 1. Compute derivates and fill in the Jacobian matrix
    double domega_dc, tmp;
    
    // 1a. Lower state row, but not if it is the ground state
    if ( ne_elev_l->ne_ilev != 0 ) {
    	// derivative of omega wrt lower state population
    	domega_dc = - k_f * c_e;
    	tmp = dGdy.get( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev );
    	dGdy.set( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev, (tmp+domega_dc) );
    }
    
    return SUCCESS;
}

int
ElectronImpactDissociation::
add_eval_contributions( Gas_data &Q, valarray<double> &G )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using the geometric average temperature
    double T_f = sqrt( Q.T[iTv] * Q.T[iTe] );
    double T_b = Q.T[iTe];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of electrons
    double c_e = Q.massf[elec->isp] * Q.rho / elec->m_w / 1.0e6;
    
    // 0c. moles per cm**3 of the atoms
    double c_A = Q.massf[atom_A->isp] * Q.rho / atom_A->m_w / 1.0e6;
    double c_B = Q.massf[atom_B->isp] * Q.rho / atom_B->m_w / 1.0e6;
    
    // 0d. moles per cm**3 of the lower state
    double c_l = ne_elev_l->elev->N / RC_Na / 1.0e6;
    
    // 1. Fill in the G vector
    
    // 1a. Lower state - if not hte ground state
    if ( ne_elev_l->ne_ilev != 0 )
    	G[ne_elev_l->ne_ilev] += k_b * c_A * c_B * c_e - k_f * c_l * c_e;
    
    return SUCCESS;
}

int
ElectronImpactDissociation::
add_source_vector_contributions( Gas_data &Q, valarray<double> &C )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using the geometric average temperature
    double T_f = sqrt( Q.T[iTv] * Q.T[iTe] );
    double T_b = Q.T[iTe];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of electrons
    double c_e = Q.massf[elec->isp] * Q.rho / elec->m_w / 1.0e6;
    
    // 0c. moles per cm**3 of the atoms
    double c_A = Q.massf[atom_A->isp] * Q.rho / atom_A->m_w / 1.0e6;
    double c_B = Q.massf[atom_B->isp] * Q.rho / atom_B->m_w / 1.0e6;
    
    // 1. Fill in the C vector
    
    // 1a. Lower state - if not the ground state
    if ( ne_elev_l->ne_ilev != 0 )
    	C[ne_elev_l->ne_ilev] -= k_b * c_A * c_B * c_e;
    
    return SUCCESS;
}

string 
ElectronImpactDissociation::
get_latex_string()
{
    string lname, lev_prefix;
    get_latex_species_piecies( rad_name, lname, lev_prefix );
    
    string lM_name = "e$^-$";
    
    string indent = "           ";
    
    ostringstream oss;
    oss << indent << "\\PIDreac{" << lname << "}{" << lM_name << "}{" << lev_prefix << "lev" << ne_elev_l->label << "}{" << atom_A->name << "}{" << atom_B->name << "}";
    oss << "\t\t&\t\t" << forward_rate_coeff->get_latex_string() << "\t\\\\" << endl;
    
    return oss.str();
}

/************************ ElectronImpactIonization **************************/

ElectronImpactIonization::
ElectronImpactIonization( lua_State *L, std::string model, Radiator * rad, NoneqElecLev * ne_elev_l )
: CR_Reaction( "ElectronImpactIonization" ),   rad( rad ), iTe( rad->iTe ), iT( rad->iT ), ne_elev_l( ne_elev_l )
{
    // Create the equation string
    equation = ne_elev_l->label + " + e_minus <=> " + rad->name + "_plus + e_minus + e_minus";
    
    if ( model=="Drawin" ) {
    	forward_rate_coeff = new DrawinElectronImpactIonization(ne_elev_l->elev->get_E(),rad->I);
    	backward_rate_coeff = new FromEquilibriumConstant();
    }
    else if ( model=="CJDrawin" ) {
    	forward_rate_coeff = new CJDrawinElectronImpactIonization(ne_elev_l->elev->get_E(),rad->I);
    	backward_rate_coeff = new FromEquilibriumConstant();
    }
    else if ( model=="Bultel" ) {
    	lua_getfield(L,-1,"bultel_data");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactIonization::ElectronImpactIonization()\n";
	    ost << "Error locating 'bultel_data' table" << endl;
	    input_error(ost);
	}
    	ostringstream oss;
    	oss << "reac_EII_ilev" << ne_elev_l->elev->i;
    	lua_getfield(L,-1,oss.str().c_str());
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactIonization::ElectronImpactIonization()\n";
	    ost << "Error locating " << oss.str() << " table" << endl;
	    input_error(ost);
	}
    	forward_rate_coeff = create_explicit_rate_coeff( L, "forward" );
    	backward_rate_coeff = create_explicit_rate_coeff( L, "backward" );
    	lua_pop(L,1); // pop reaction
    	lua_pop(L,1); // pop bultel_data
    }
    else if ( model=="SunoKato" ) {
    	lua_getfield(L,-1,"suno_kato_data");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactIonization::ElectronImpactIonization()\n";
	    ost << "Error locating 'suno_kato_data' table" << endl;
	    input_error(ost);
	}
    	ostringstream oss;
    	oss << "reac_EII_ilev" << ne_elev_l->elev->i;
    	lua_getfield(L,-1,oss.str().c_str());
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactIonization::ElectronImpactIonization()\n";
	    ost << "Error locating " << oss.str() << " table" << endl;
	    input_error(ost);
	}
    	forward_rate_coeff = new SunoKatoCarbonElectronImpactIonization( L, ne_elev_l->elev, rad->I );
    	backward_rate_coeff = new FromEquilibriumConstant();
    	lua_pop(L,1); // pop reaction
    	lua_pop(L,1); // pop suno_kato_data
    }
    else if ( model=="KuncSoon" ) {
    	lua_getfield(L,-1,"kunc_soon_data");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactIonization::ElectronImpactIonization()\n";
	    ost << "Error locating 'kunc_soon_data' table" << endl;
	    input_error(ost);
	}
    	ostringstream oss;
    	oss << "reac_EII_ilev" << ne_elev_l->elev->i;
    	lua_getfield(L,-1,oss.str().c_str());
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "ElectronImpactIonization::ElectronImpactIonization()\n";
	    ost << "Error locating " << oss.str() << " table" << endl;
	    input_error(ost);
	}
    	forward_rate_coeff = new KuncSoonElectronImpactIonization( L, ne_elev_l->elev, rad->I );
    	backward_rate_coeff = new FromEquilibriumConstant();
    	lua_pop(L,1); // pop reaction
    	lua_pop(L,1); // pop kunc_soon_data
    }
    else {
    	ostringstream ost;
    	ost << "ElectronImpactIonization::ElectronImpactIonization()" << endl
    	    << "The " << model << " model is not implemented." << endl;
    	input_error(ost);
    }
}

double
ElectronImpactIonization::
eval_equilibrium_constant( double T )
{
    // NOTE: not including h_f in partition functions
    double Q_ion = ion->calc_total_equil_partition_function(T);
    double Q_elec = elec->calc_total_equil_partition_function(T);
    double Q_elev = ne_elev_l->elev->calculate_equilibrium_Q_total(T) * rad->eval_translational_partition_function_from_T(T) * exp( rad->I / RC_k_SI / T );
    
    // Kc = PI(Q)_products / PI(Q)_reactants
    
    double Kc = Q_ion * Q_elec / Q_elev;
    
    Kc /= RC_Na * 1.0e6;		// Convert particles / m**3 -> moles / cm**3
    
    return Kc;
}

int
ElectronImpactIonization::
add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using the electronic temperature
    double T_f = Q.T[iTe];
    double T_b = Q.T[iTe];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of electrons
    double c_e = Q.massf[elec->isp] * Q.rho / elec->m_w / 1.0e6;
    
    // 1. Compute derivates and fill in the Jacobian matrix
    double domega_dc, tmp;
    
    // 1a. Lower state row, but not if it is the ground state
    if ( ne_elev_l->ne_ilev != 0 ) {
    	// derivative of omega wrt lower state population
    	domega_dc = - k_f * c_e;
    	tmp = dGdy.get( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev );
    	dGdy.set( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev, (tmp+domega_dc) );
    }
    
    return SUCCESS;
}

int
ElectronImpactIonization::
add_eval_contributions( Gas_data &Q, valarray<double> &G )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using the electronic temperature
    double T_f = Q.T[iTe];
    double T_b = Q.T[iTe];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of electrons
    double c_e = Q.massf[elec->isp] * Q.rho / elec->m_w / 1.0e6;
    
    // 0c. moles per cm**3 of the ion
    double c_i = Q.massf[ion->isp] * Q.rho / ion->m_w / 1.0e6;
    
    // 0d. moles per cm**3 of the lower state
    double c_l = ne_elev_l->elev->N / RC_Na / 1.0e6;
    
    // 1. Fill in the G vector
    
    // 1a. Lower state - if not hte ground state
    if ( ne_elev_l->ne_ilev != 0 )
    	G[ne_elev_l->ne_ilev] += k_b * c_i * c_e * c_e - k_f * c_l * c_e;
    
    return SUCCESS;
}

int
ElectronImpactIonization::
add_source_vector_contributions( Gas_data &Q, valarray<double> &C )
{
    // 0. Prepare data
    // 0a. Forward reverse reaction rates using the electronic temperature
    double T_f = Q.T[iTe];
    double T_b = Q.T[iTe];
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 0b. moles per cm**3 of electrons
    double c_e = Q.massf[elec->isp] * Q.rho / elec->m_w / 1.0e6;
    
    // 0c. moles per cm**3 of the ion
    double c_i = Q.massf[ion->isp] * Q.rho / ion->m_w / 1.0e6;
    
    // 1. Fill in the G vector
    
    // 1a. Lower state - if not the ground state
    if ( ne_elev_l->ne_ilev != 0 )
    	C[ne_elev_l->ne_ilev] -= k_b * c_i * c_e * c_e;
    
    return SUCCESS;
}

string 
ElectronImpactIonization::
get_latex_string()
{
    return "";
}

/************************** RadiativeTransition ****************************/

RadiativeTransition::
RadiativeTransition( lua_State * L, Radiator * rad )
: CR_Reaction( "RadiativeTransition" )
{
    QSSDiatomicRadiator * QSS_rad = dynamic_cast<QSSDiatomicRadiator*>(rad);
    
    equation = get_string( L, -1, "equation" );
    vector<string> equation_tks;
    tokenize_equation_string( equation, equation_tks );
    
    ne_elev_u = 0; ne_elev_l = 0;
    
    for ( size_t i=0; i<QSS_rad->noneq_elevs.size(); ++i ) {
    	if ( QSS_rad->name + "_" + QSS_rad->noneq_elevs[i]->label == equation_tks[0] )
    	    ne_elev_u = QSS_rad->noneq_elevs[i];
    	if ( QSS_rad->name + "_" + QSS_rad->noneq_elevs[i]->label == equation_tks[2] )
    	    ne_elev_l = QSS_rad->noneq_elevs[i];
    }
    
    if ( ne_elev_u==0 ) {
    	ostringstream oss;
    	oss << "RadiativeTransition::RadiativeTransition()" << endl
    	    << "Upper level: " << equation_tks[0] << " not found." << endl;
    	input_error( oss );
    }
    if ( ne_elev_l==0 ) {
    	ostringstream oss;
    	oss << "RadiativeTransition::RadiativeTransition()" << endl
    	    << "Lower level: " << equation_tks[2] << " not found." << endl;
    	input_error( oss );
    }
    
    forward_rate_coeff = create_explicit_rate_coeff( L, "forward" );
    backward_rate_coeff = create_explicit_rate_coeff( L, "backward" );

}

RadiativeTransition::
RadiativeTransition( std::string model, NoneqElecLev * ne_elev_u, NoneqElecLev * ne_elev_l, double A_ul )
: CR_Reaction( "RadiativeTransition" ), ne_elev_u( ne_elev_u ), ne_elev_l( ne_elev_l )
{
    // Create the equation string
    equation = ne_elev_u->label + " <=> " + ne_elev_u->label + " + hv";
    
    // Decay constant
    double tau = 1.0 / A_ul;
    
    // Forward rate coefficient model
    if ( model=="OpticallyThin" ) {
    	forward_rate_coeff = new OpticallyThinExponentialDecay( tau );
    }
    else if ( model=="OpticallyVariable" ) {
    	double wavel_nm = nu2lambda( ( ne_elev_u->elev->get_E() - ne_elev_l->elev->get_E() ) / RC_h_SI );
    	// use default parameters for the moments
    	forward_rate_coeff = new OpticallyVariableExponentialDecay( tau, wavel_nm );
    }
    else if ( model=="CurveFit" ) {
    	forward_rate_coeff = new CurveFitExponentialDecay( tau );
    }
    else {
    	ostringstream oss;
    	oss << "RadiativeTransition::RadiativeTransition()" << endl
    	    << "Model: " << model << " not recognised." << endl;
    	input_error( oss );
    }
    
    // Dummy backward rate coefficient model
    backward_rate_coeff = new ZeroRate();
}

int
RadiativeTransition::
eval_reaction_rates( double T_f, double T_b, Gas_data &Q, double &k_f, double &k_b )
{
    k_f = forward_rate_coeff->get_rate( T_f, Q );
    k_b = 0.0;
    
    return SUCCESS;
}

double
RadiativeTransition::
eval_equilibrium_constant( double T )
{
    UNUSED_VARIABLE(T);
    
    cout << "RadiativeTransition::eval_equilibrium_constant()" << endl
         << "This function is not meant for use, bailing out!" << endl;
    exit( FAILURE );
}

int
RadiativeTransition::
add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy )
{
    // 0. Prepare data
    // 0a. Forward reaction rates using dummy temperatures
    double T_f = 0.0;
    double T_b = 0.0;
    double k_f, k_b;
    this->eval_reaction_rates( T_f, T_b, Q, k_f, k_b );
    
    // 1. Compute derivates and fill in the Jacobian matrix
    //    NOTE: - the derivatives wrt to the lower state currently should equate to zero (k_b=0)
    //          - ...but should implemented a correct method to make k_b non-zero in the presence of absorption
    double domega_dc, tmp;
    
    // 1b. Upper state row (should never be the ground state)
    // derivative of omega wrt lower state population
    domega_dc = k_b;
    tmp = dGdy.get( ne_elev_u->ne_ilev, ne_elev_l->ne_ilev );
    dGdy.set( ne_elev_u->ne_ilev, ne_elev_l->ne_ilev, (tmp+domega_dc) );
    // derivative of omega wrt upper state population
    domega_dc = - k_f;
    tmp = dGdy.get( ne_elev_u->ne_ilev, ne_elev_u->ne_ilev );
    dGdy.set( ne_elev_u->ne_ilev, ne_elev_u->ne_ilev, (tmp+domega_dc) );
    
    // 1a. Lower state row, but not if it is the ground state
    if ( ne_elev_l->ne_ilev != 0 ) {
    	// derivative of omega wrt lower state population
    	domega_dc = - k_b;
    	tmp = dGdy.get( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev );
    	dGdy.set( ne_elev_l->ne_ilev, ne_elev_l->ne_ilev, (tmp+domega_dc) );
    	// derivative of omega wrt upper state population
    	domega_dc = k_f;
    	tmp = dGdy.get( ne_elev_l->ne_ilev, ne_elev_u->ne_ilev );
    	dGdy.set( ne_elev_l->ne_ilev, ne_elev_u->ne_ilev, (tmp+domega_dc) );
    }
    
    return SUCCESS;
}

int
RadiativeTransition::
add_eval_contributions( Gas_data &Q, valarray<double> &G )
{
    return SUCCESS;
}

int
RadiativeTransition::
add_source_vector_contributions( Gas_data &Q, valarray<double> &C )
{
    // Nothing to be done
    return SUCCESS;
}

string 
RadiativeTransition::
get_latex_string()
{
    string lname, lev_prefix;
    get_latex_species_piecies( "CO", lname, lev_prefix );
    
    string indent = "           ";
    
    ostringstream oss;
    oss << indent << "\\PIDreac{" << lname << "}{" << lev_prefix << "lev" << ne_elev_l->label << "}{" << lev_prefix << "lev" << ne_elev_u->label << "}";
    oss << "\t\t&\t\t" << forward_rate_coeff->get_latex_string() << "\t\\\\" << endl;
    
    return oss.str();
}

/************************** Helper functions *****************************/


void tokenize_equation_string( string equation, vector<string> &tks )
{
    stringstream ss(equation);
    string s;
    while (getline(ss, s, ' ')) tks.push_back( s );
    
    return;
}

void get_latex_species_piecies( string name, string &lname, string &lev_prefix )
{
    if ( name=="Ar" ) {
    	lname = "Ar";
    	lev_prefix = "\\Ar";
    }
    else if ( name=="Ar_plus" ) {
    	lname = "Ar$^+$";
    	lev_prefix = "\\Arp";
    }
    else if ( name=="C" ) {
    	lname = "C";
    	lev_prefix = "\\C";
    }
    else if ( name=="C_plus" ) {
    	lname = "C$^+$";
    	lev_prefix = "\\Cp";
    }
    else if ( name=="N" ) {
    	lname = "N";
    	lev_prefix = "\\N";
    }
    else if ( name=="N_plus" ) {
    	lname = "N$^+$";
    	lev_prefix = "\\Np";
    }
    else if ( name=="O" ) {
    	lname = "O";
    	lev_prefix = "\\O";
    }
    else if ( name=="O_plus" ) {
    	lname = "O$^+$";
    	lev_prefix = "\\Op";
    }
    else if ( name=="C2" ) {
    	lname = "C$_2$";
    	lev_prefix = "\\CC";
    }
    else if ( name=="CN" ) {
    	lname = "CN";
    	lev_prefix = "\\CN";
    }
    else if ( name=="CO" ) {
    	lname = "CO";
    	lev_prefix = "\\CO";
    }
    else if ( name=="CO_plus" ) {
    	lname = "CO$^+$";
    	lev_prefix = "\\COp";
    }
    else if ( name=="N2" ) {
    	lname = "N$_2$";
    	lev_prefix = "\\NN";
    }
    else if ( name=="N2_plus" ) {
    	lname = "N$_2^+$";
    	lev_prefix = "\\NNp";
    }
    else if ( name=="NO" ) {
    	lname = "NO";
    	lev_prefix = "\\NO";
    }
    else if ( name=="NO_plus" ) {
    	lname = "NO$^+$";
    	lev_prefix = "\\NOp";
    }
    else if ( name=="O2" ) {
    	lname = "O$_2$";
    	lev_prefix = "\\OO";
    }
    else if ( name=="O2_plus" ) {
    	lname = "O$_2^+$";
    	lev_prefix = "\\OOp";
    }
    else if ( name=="e-" || name=="e" ) {
    	lname = "e$^-$";
    	lev_prefix = "\\em";
    }
    else {
    	cout << "get_latex_species_piecies()" << endl
    	     << "species: " << name << " not in the library" << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    return;
}
