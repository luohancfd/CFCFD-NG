/** \file coupled_noneq_system.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 24-Mar-09: First attempt at coupled QSS model
 *  \brief Definitions for coupled nonequilibirum population model
 *
 **/

#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "qss_noneq_system.hh"
#include "radiation_constants.hh"
#include "../../util/source/useful.h"
#include "amalgam.hh"

using namespace std;

QSSNoneqSystem::QSSNoneqSystem( const string rad_name, const string input_file ) 
{
    if ( ECHO_INPUT > 0 ) {
	cout << "-> Creating a QSSNoneqSystem for " << rad_name << endl;
    }
    
    ConfigParser cfg( input_file );
    string section = rad_name + "-QSS-data";
    
    ne_rad_ = create_new_noneq_radiator( rad_name, section, &cfg );
    size_t n_states = ne_rad_->elec_states_.size();
    
    // Fill in ne_index_'s - just the level number
    for ( size_t ilev=0; ilev<n_states; ++ilev ) {
	ne_rad_->elec_states_[ilev]->ne_index_ = ilev;
    }
    
    // Initialise the (stored) Jacobian matrix, sets elements to zero
    dGdy_.resize( n_states, n_states );
    // and the source vector...
    C_.resize( n_states, 0.0 );
    // and the output vector...
    y_out_.resize( n_states, 0.0 );
    
    if ( ECHO_INPUT > 1 ) {
	cout << "- Boltzmann equilibriations of non-nonequilibrium levels requested." << endl;
    }
    
    // Reactions are initialised separately so that all radiators are pre-initialised
}

QSSNoneqSystem::~QSSNoneqSystem()
{
    // delete NoneqRadiator pointed at by ne_rad_ pointer
    delete ne_rad_;
    
    for ( size_t ireac=0; ireac<reactions_.size(); ++ireac )
	delete reactions_[ireac];
}

int
QSSNoneqSystem::complete_initialisation( ConfigParser *cfg )
{
    ne_rad_->irad_ = get_irad_from_name( ne_rad_->name_ );
    ne_rad_->rad_pointer_ = get_radiator_pointer( ne_rad_->irad_ );
    if ( ECHO_INPUT > 0 )
	cout << "Completing QSS initialisation for radiator: " << ne_rad_->name_ << endl;
    
    // auto-create all electronic levels as noneq levels if requested
    if ( ne_rad_->all_levels_noneq_ ) {
	if ( ECHO_INPUT > 0 ) {
	    cout << "- All " << ne_rad_->rad_pointer_->ie_n << " levels for radiator: " 
	         << ne_rad_->rad_pointer_->rad_name << " are noneq states" << endl;
        }
	for ( int ie=0; ie<ne_rad_->rad_pointer_->ie_n; ++ie ) {
	    ostringstream label;
	    label << "e" << ie;
	    ne_rad_->elec_states_.push_back( new NoneqElecState( ie, label.str() ) );
	}
    }
    
    // Initialise boltz_fractions_ if required
    // NOTE: -1's will remain in base positions for whole simulations
    if ( ne_rad_->boltz_eqs_ ) ne_rad_->boltz_fractions_.resize(ne_rad_->rad_pointer_->ie_n, -1.0);
    
    string section = ne_rad_->rad_pointer_->rad_name + "-QSS-data/reactions"; 

    int n_explicit_reactions;
    cfg->parse_int( section, "n_explicit_reactions", n_explicit_reactions, 0 );
    
    for ( int ireac=0; ireac<n_explicit_reactions; ++ireac ) {
	ostringstream reac_section;
	reac_section << section << "/reac-" << ireac;
	reactions_.push_back( create_new_explicit_CR_reaction( reac_section.str(), cfg ) );
    }
    
    bool auto_create_reactions;
    cfg->parse_boolean( section, "auto_create_reactions", auto_create_reactions, false );
    
    if ( auto_create_reactions ) {
	// NOTE: only considering these three processes for the moment
	vector<string> reac_types;
	reac_types.push_back( "electron_impact_excitation" );
	reac_types.push_back( "electron_impact_ionization" );
	reac_types.push_back( "radiative_transition" );
	vector<string> models(reac_types.size());
	for ( size_t im=0; im<reac_types.size(); ++im ) {
	    string key = reac_types[im] + "_model";
	    cfg->parse_string( section, key, models[im], "default" );
	}
	autocreate_implicit_CR_reactions_for_radiator( cfg, ne_rad_, reac_types, models, reactions_ );
    }
    
    // Echo reactions to screen if requested
    if ( ECHO_INPUT > 0 ) {
	for ( size_t ireac=0; ireac<reactions_.size(); ++ireac ) {
	    cout << "- CR_Reaction " << ireac << ": " << reactions_[ireac]->get_type_str()
	         << " - " << reactions_[ireac]->get_equation_str() << endl;
	}
    }
    
    if ( reactions_.size()==0 ) {
	cout << "QSSNoneqSystem::complete_initialisation()" << endl
	     << "No reactions have been initialised for radiator:" << ne_rad_->rad_pointer_->rad_name << endl
	     << "Bailing out!" << endl;
	exit( BAD_INPUT_ERROR );
    }
    else {
	cout << "QSSNoneqSystem::complete_initialisation()" << endl
	     << reactions_.size() << " reactions have been initialised for radiator:" << ne_rad_->rad_pointer_->rad_name << endl;
    }
    
    return SUCCESS;
}

int
QSSNoneqSystem::solve_system(Gas_data &Q)
{
    /* Linear approach, solve system via Gaussian Elimination */
    
    size_t n_states = ne_rad_->elec_states_.size();
    
    // 0. Setup Jacobian matrix and source and solution vectors
    // Now part of the class - just clear all elements
    // CHECK-ME: is it better just to create new instances each time?
    for ( size_t i=0; i<n_states; ++i ) {
	C_[i] = 0.0;
	y_out_[i] = 0.0;
	for ( size_t j=0; j<n_states; ++j ) {
	    dGdy_.set( i, j, 0.0 );
	}
    }
    
    // 1.  Population summations, first matrix row
    if ( ne_rad_->boltz_eqs_ ) {
	// Include boltzmann fractions from non-noneq levels
	// NOTE: also doing direct contributions to prevent repeating loops
	//       and saving the fractions for use at the end 
	for ( size_t is=0; is<n_states; ++is ) {
	    int ie_base = ne_rad_->elec_states_[is]->index_;
	    int ie_limit;
	    if ( is==(n_states-1) ) ie_limit = ne_rad_->rad_pointer_->ie_n;
	    else ie_limit = ne_rad_->elec_states_[is+1]->index_;
	    // Assume Te is the electronic temperature
	    double Te = Q.T_e;
	    double Q_base = ne_rad_->rad_pointer_->calculate_Q_for_state( Te, ie_base );
	    double bf_acc = 0.0;
	    for ( int ie=ie_base+1; ie<ie_limit; ++ie ) {
		double bf = ne_rad_->rad_pointer_->calculate_Q_for_state( Te, ie ) / Q_base;
		ne_rad_->boltz_fractions_[ie] = bf;
		bf_acc += bf;
	    }
	    double tmp = dGdy_.get(0,is) + 1.0 + bf_acc;
	    dGdy_.set(0,is,tmp);
	}
    }
    else {
	// Just Direct contributions from noneq levels
	for ( size_t is=0; is<n_states; ++is ) dGdy_.set(0,is,1.0);
    }
    
    // 2. Contributions from reactions
    for ( size_t ir=0; ir<reactions_.size(); ++ir )
	reactions_[ir]->add_jacobian_contributions( Q, dGdy_ );
    
    // 3. Construct source vector
    C_[0] = Q.c[ne_rad_->rad_pointer_->isp] / 1.0e6;		// Convert moles/m**3 to moles/cm**3
    
    // 4. Solve the system
    if( gaussian_elimination( dGdy_, y_out_, C_ ) ) {
        cout << "QSSNoneqSystem::solve_system() - gaussian_elimination failed!" << endl;
        return FAILURE;
    }
    
    // 5.  Map results back onto radiator
    // 5a. Firstly noneq states
    for ( size_t is=0; is<n_states; ++is ) {
	int ie = ne_rad_->elec_states_[is]->index_;
	ne_rad_->rad_pointer_->N_el[ie] = y_out_[is] * 1.0e6 * RC_Na;	// Convert moles/cm**3 -> particles/m**3
    }
    
    // 5b. non-noneq states
    if ( ne_rad_->boltz_eqs_ ) {
	int ie_base = 0;
	for ( int ie=0; ie<ne_rad_->rad_pointer_->ie_n; ++ie ) {
	    if ( ne_rad_->boltz_fractions_[ie]<0.0 ) {
		// This is a base
		ie_base = ie;
	    }
	    else {
		// This level requires equilibriation
		ne_rad_->rad_pointer_->N_el[ie] = ne_rad_->rad_pointer_->N_el[ie_base] * ne_rad_->boltz_fractions_[ie];
	    }
	}
    }
    
#   if DEBUG_RAD > 0
    // Check that the total species density has been conserved
    double n_sum = 0.0;
    for ( int ie=0; ie<ne_rad_->rad_pointer_->ie_n; ++ie ) {
	n_sum += ne_rad_->rad_pointer_->N_el[ie];
    }
    cout << "QSSNoneqSystem::solve_system()" << endl
         << "Radiator: " << ne_rad_->rad_pointer_->rad_name << endl
         << "n(QSS) = " << n_sum << ", n(CFD) = " << Q.c[ne_rad_->rad_pointer_->isp] * RC_Na << endl;
#   endif
    
    return SUCCESS;
}

