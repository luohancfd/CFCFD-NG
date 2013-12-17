// Author: Daniel F. Potter
// Date: 06-Oct-2009
// Place: Brisbane, Queendland, AUST
#include <cstdlib>
#include <cmath>
#include <sstream>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"

#include "chemical-equilibrium-system.hh"
#include "gas-model.hh"

using namespace std;

Partial_equilibrium_reaction::
Partial_equilibrium_reaction( vector<int> betas, vector<Chemical_species*> species )
{
    betas_.assign( betas.begin(), betas.end() );
    species_.assign( species.begin(), species.end() );
}

string
Partial_equilibrium_reaction::
str()
{
    string reactants = "", products = "";
    for ( size_t isp=0; isp<species_.size(); ++isp ) {
    	string participant = "";
	for ( int i=0; i<abs(betas_[isp]); ++i ) 
	    participant += species_[isp]->get_name() + " + ";
    	if ( betas_[isp] < 0 ) reactants += participant;
        else if ( betas_[isp] > 0 ) products += participant;
    }
    
    reactants.erase( reactants.size() - 3, 3 );
    products.erase( products.size() - 3, 3 );
    
    return reactants + " <=> " + products;
}

int
Partial_equilibrium_reaction::
G_sum()
{
    int G_sum = 0;
    for ( size_t irsp=0; irsp<betas_.size(); ++irsp )
    	G_sum -= betas_[irsp];
    
    return G_sum;
}

double
Partial_equilibrium_reaction::
s_eval_K_p( double T )
{
    // NOTE: this the equilibrium constant for partial pressures in atmospheres
    
    double tmp = 0.0;
    for ( size_t irsp=0; irsp<species_.size(); ++irsp ) {
    	double G = species_[irsp]->eval_gibbs_free_energy(T);
    	double R = species_[irsp]->get_R();
    	tmp += double(betas_[irsp]) * G / ( R * T );
    }
    
    return exp( - tmp );
}

double
Partial_equilibrium_reaction::
s_eval_K_c( double T )
{
    // NOTE: this the equilibrium constant for molar concentration
    
    double K_p = s_eval_K_p(T);
    double K_c = K_p*pow(PC_P_atm/(PC_R_u*T),-G_sum());
    
    return K_c;
}

double
Partial_equilibrium_reaction::
s_eval_equilibrium_constant( double T )
{
    // NOTE: this the equilibrium constant for the log of mass density
    
    double K_m = 0.0;
    for ( size_t irsp=0; irsp<species_.size(); ++irsp ) {
    	double G = species_[irsp]->eval_gibbs_free_energy(T);
    	double R = species_[irsp]->get_R();
    	K_m += double(betas_[irsp]) * ( G / ( R * T ) + log( R * T / PC_P_atm ) );
    }
    
    return K_m;
}

Chemical_equilibrium_system::
Chemical_equilibrium_system( double min_massf, vector<Chemical_species*> &species )
 : nsp_( species.size() ), min_massf_( min_massf )
{
    // 0. Determine the set of base elements and create partial equilibrium reactions
    ions_ = false;
    for ( size_t isp=0; isp<nsp_; ++isp ) {
    	if ( species[isp]->get_Z()!=0 ) ions_ = true;
    	vector<int> betas;
    	vector<string> participants;
    	species[isp]->partial_equilibrium_participants( betas, participants );
    	// If this species a base-element then two identical participants will be present
    	if ( participants.size()==2 && participants[0]==participants[1] ) {
    	    cout << "- " << species[isp]->get_name() << " is a base element" << endl;
    	    be_indices_.push_back( isp );
    	}
    	else if ( species[isp]->get_name()!="e_minus" ) {
    	    // Create a partial equilibrium reaction (e- has charge balance)
    	    vector<Chemical_species*> participant_pointers;
	    for ( size_t ip=0; ip<participants.size(); ++ip ) {
		for ( size_t jsp=0; jsp<species.size(); ++jsp ) {
		    if ( participants[ip].compare( species[jsp]->get_name() )==0 ) {
			participant_pointers.push_back( species[jsp] );
			break;
		    }
		    else if ( jsp==(species.size()-1) ) {
		    	ostringstream ost;
		    	ost << "Chemical_equilibrium_system::Chemical_equilibrium_system()\n";
		    	ost << "Base element: " << participants[ip] << " is not a species" << endl;
		    	input_error(ost);
		    }
		}
	    }
	    pe_reactions_.push_back( new Partial_equilibrium_reaction( betas, participant_pointers ) );
	    cout << "- Successfully created new partial equilibrium reaction: " << pe_reactions_.back()->str() << endl;
	}
    }
    
    // 1. Calculate base-element weighting matrix
    be_weightings_.resize( be_indices_.size() );
    for ( size_t ibe = 0; ibe < be_indices_.size(); ++ibe ) {
    	string be_name = species[ be_indices_[ibe] ]->get_name();
    	be_weightings_[ibe].resize( nsp_ );
    	for ( size_t isp=0; isp<nsp_; ++isp ) {
    	    be_weightings_[ibe][isp] = (double) species[isp]->get_element_count( be_name ) \
    	    	* species[be_indices_[ibe]]->get_M() / species[isp]->get_M();
    	}
    }
    
    // 2. Calculate charge weighting vector
    // double M_elec = 0.000548579903e-3;
    charge_weightings_.resize( nsp_ );
    for ( size_t isp=0; isp<nsp_; ++isp ) {
    	charge_weightings_[isp] = (double) species[isp]->get_Z() \
    		* 1.0 / species[isp]->get_M();
    }
    
    // 3. Initialise vectors
    Q_.resize( nsp_ );
    yguess_.resize( nsp_ );
    yout_.resize( nsp_ );
    
    // 4. Initialise zero-solver
    double tol = 1.0e-12;
    int max_iter = 100000;
    bool use_jacobian = true;
    zero_solver_.set_constants(nsp_, tol, max_iter, use_jacobian);
    // f_jac_ is the Jacobian scale factor
    f_jac_ = 1.0;
   
}

Chemical_equilibrium_system::~Chemical_equilibrium_system()
{
    for ( size_t ir=0; ir<pe_reactions_.size(); ++ir )
    	delete pe_reactions_[ir];
}

int
Chemical_equilibrium_system::solve_system( double T, double rho, vector<double> &massf )
{  
    // FIXME: just for testing purposes
    
    // 0. Store local copy of T and log(rho)
    log_rho_ = log(rho);

    // 1. Calculate and store mass density equilibrium constants
    for ( size_t ir=0; ir<pe_reactions_.size(); ++ir )
    	pe_reactions_[ir]->store_equilibrium_constant(T);
    
    // 2. Compute (constant) source terms
    compute_source_terms( massf );

    // 3. Fill out yguess vector
    for ( size_t isp=0; isp<nsp_; ++isp )
    	yguess_[isp] = ( massf[isp] > min_massf_ ) ? log(massf[isp]) : log(min_massf_);
    
    // 4. Solve the system of equations with Rowan's Newton iterator
    if ( zero_solver_.solve( *this, yguess_, yout_ ) ) {
	cout << "Chemical_equilibrium_system::solve_system()" << endl
	     << "Zero solver has failed, bailing out!" << endl;
	exit( FAILURE );
    }
    
    // 5.  Map results back onto massf vector
    for ( size_t isp=0; isp<nsp_; ++isp )
    	massf[isp] = ( yout_[isp] > log(min_massf_) ) ? exp(yout_[isp]) : 0.0;
    
    return SUCCESS;
}

#define WITH_MASSF_SUM 0

int
Chemical_equilibrium_system::compute_source_terms( vector<double> &massf )
{
    // Initialise Q_ vector with zeros
    for ( size_t isp=0; isp<nsp_; ++isp ) {
	 Q_[isp] = 0.0;
    }
    
    int iQ=0;	// current element index for the Q vector
    
    // 1. charge conservation line (if ions_ flag is set to true)
    if ( ions_ ) {
    	// Count the charge rather than set to 0, as charge neutrality may not
    	// be enforced by the CFD solver
    	for ( size_t isp=0; isp<nsp_; ++isp ) {
    	    Q_[iQ] += massf[isp] * charge_weightings_[isp];
    	    iQ++;
    	}
    }
    
#   if WITH_MASSF_SUM
    // 2a. mass-fraction summation line
    Q_[iQ] = 1.0;
    iQ++;
    
    // 2b. elemental conservation lines (one for each base-element, but not the first one)
    for ( size_t ibe=1; ibe<be_indices_.size(); ++ibe ) {
    	for ( size_t isp=0; isp<nsp_; ++isp ) {
    	    Q_[iQ] += massf[isp] * be_weightings_[ibe][isp];
    	}
    	iQ++;
    }
#   else
    // 2. elemental conservation lines (one for each base-element)
    for ( size_t ibe=0; ibe<be_indices_.size(); ++ibe ) {
    	for ( size_t isp=0; isp<nsp_; ++isp ) {
    	    Q_[iQ] += massf[isp] * be_weightings_[ibe][isp];
    	}
    	iQ++;
    }
#   endif
    
    // 3. Equilibrium constants for non-base-element species
    for ( size_t ir=0; ir<pe_reactions_.size(); ++ir ) {
    	Q_[iQ] += pe_reactions_[ir]->get_K_m();
    	// Now subtract out weighted ln(rho) term to use mass-fractions
    	Q_[iQ] -= double(pe_reactions_[ir]->G_sum())*log_rho_;
    	iQ++;
    }
    
#   if 0
    // print out the Q vector to the screen
    cout << "Chemical_equilibrium_system::compute_source_terms()" << endl;
    for ( size_t isp=0; isp<nsp_; ++isp ) {
	cout << "Q[" << isp << "] = " << Q_[isp] << endl;
    }
#   endif

    return SUCCESS;
}

int
Chemical_equilibrium_system::f( const vector<double> &y, vector<double> &G )
{
    /* Create the equation system for the ZeroSystem for a given y vector */
    
    // 0.  Apply source terms (as negatives as we are creating a zero system)
    for ( size_t isp=0; isp<nsp_; ++isp ) {
	 G[isp] = -Q_[isp];
    }
    
    int iG=0;	// current element index for the G vector
    
    // 1. charge conservation line (if ions_ flag is set to true)
    if ( ions_ ) {
    	for ( size_t isp=0; isp<nsp_; ++isp ) {
    	    G[iG] += exp(y[isp]) * charge_weightings_[isp];
    	}
    	iG++;
    }
    
#   if WITH_MASSF_SUM
    // 2a. mass-fraction summation line
    for ( size_t isp=0; isp<nsp_; ++isp )
    	G[iG] += exp(y[isp]);
    iG++;
    
    // 2b. elemental conservation lines (one for each base-element, but not the first one)
    for ( size_t ibe=1; ibe<be_indices_.size(); ++ibe ) {
    	for ( size_t isp=0; isp<nsp_; ++isp ) {
    	    G[iG] += exp(y[isp]) * be_weightings_[ibe][isp];
    	}
    	iG++;
    }
#   else
    // 2. elemental conservation lines (one for each base-element)
    for ( size_t ibe=0; ibe<be_indices_.size(); ++ibe ) {
    	for ( size_t isp=0; isp<nsp_; ++isp ) {
    	    G[iG] += exp(y[isp]) * be_weightings_[ibe][isp];
    	}
    	iG++;
    }
#   endif

    // 3. Equilibrium constants for non-base-element species
    for ( size_t ir=0; ir<pe_reactions_.size(); ++ir ) {
    	for ( int irsp=0; irsp<pe_reactions_[ir]->get_nrsp(); ++irsp ) {
    	    int isp = pe_reactions_[ir]->get_isp( irsp );
    	    G[iG] += pe_reactions_[ir]->get_G( irsp ) * y[isp];
    	}
    	iG++;
    }
    
#   if 0
    // print out the G vector to the screen
    cout << "Chemical_equilibrium_system::f()" << endl;
    for ( size_t is=0; is<nsp_; ++is ) {
	cout << "G[" << is << "] = " << G[is] << endl;
    }
#   endif

    return SUCCESS;
}

int
Chemical_equilibrium_system::Jac( const vector<double> &y, Valmatrix &dGdy )
{
    /* Create the Jacobian matrix for the ZeroSystem for a given y vector */
    
    // 0.  Clear the jacobian matrix
    for ( size_t i=0; i<nsp_; ++i ) {
	for ( size_t j=0; j<nsp_; ++j ) {
	    dGdy.set(i,j,0.0);
	}
    }

    int iG = 0;		// current matrix line

    // 1. charge conservation line (if ions_ flag is set to true)
    if ( ions_ ) {
    	for ( size_t isp=0; isp<nsp_; ++isp )
    	    dGdy.set(iG,isp,exp(y[isp])*charge_weightings_[isp]*f_jac_);
    	iG++;
    }
    
#   if WITH_MASSF_SUM
    // 2a. mass-fraction summation line
    for ( size_t isp=0; isp<nsp_; ++isp )
    	dGdy.set(iG,isp,1.0*f_jac_);
    iG++;
    
    // 2b. elemental conservation lines (one for each base-element, but not the first one)
    for ( size_t ibe=1; ibe<be_indices_.size(); ++ibe ) {
    	for ( size_t isp=0; isp<nsp_; ++isp )
    	    dGdy.set(iG,isp,exp(y[isp])*be_weightings_[ibe][isp]*f_jac_);
    	iG++;
    }
#   else
    // 2. elemental conservation lines (one for each base-element)
    for ( size_t ibe=0; ibe<be_indices_.size(); ++ibe ) {
    	for ( size_t isp=0; isp<nsp_; ++isp )
    	    dGdy.set(iG,isp,exp(y[isp])*be_weightings_[ibe][isp]*f_jac_);
    	iG++;
    }
#   endif
    
    // 3. Equilibrium constant line (one for every non-base-element species)
    for ( size_t ir=0; ir<pe_reactions_.size(); ++ir ) {
    	for ( int irsp=0; irsp<pe_reactions_[ir]->get_nrsp(); ++irsp ) {
    	    int isp = pe_reactions_[ir]->get_isp( irsp );
    	    double G_isp = pe_reactions_[ir]->get_G( irsp );
    	    // d sum( G_i * ln( f_i ) ) / d f_i = G_i / f_i
    	    dGdy.set(iG,isp,G_isp*f_jac_);
    	}
    	iG++;
    }
    
#   if 0
    // print out the Jacobian matrix to the screen
    cout << "Chemical_equilibrium_system::Jac()" << endl
         << "Jacobian matrix dGdy: " << endl
         << dGdy.str() << endl;
#   endif
    
    return SUCCESS;
}

int
Chemical_equilibrium_system::
test_system( double T, double p, std::vector<double> molef )
{
    cout << "T = " << T << "\t";
    for ( size_t ir=0; ir<pe_reactions_.size(); ++ir ) {
    	// compute K_p from thermo relations
    	double K_p_thermo = pe_reactions_[ir]->eval_K_p(T);
    	// compute K_p from given mole-fractions
    	double K_p_given = 1.0;
    	for ( int irsp=0; irsp<pe_reactions_[ir]->get_nrsp(); ++irsp ) {
    	    int isp = pe_reactions_[ir]->get_isp( irsp );
    	    double beta = pe_reactions_[ir]->get_beta(irsp);
    	    K_p_given *= pow( molef[isp]*p/PC_P_atm, beta);
    	}
    	// cout << pe_reactions_[ir]->str() << ": " << fabs(K_p_thermo-K_p_given)/K_p_given;
    	cout << fabs(K_p_thermo-K_p_given)/K_p_given << " \t ";
    }
    
    cout << endl;
    
    return SUCCESS;
}

No_chemical_equilibrium_system::
No_chemical_equilibrium_system( double min_massf, vector<Chemical_species*> &species )
{
    cout << "No_chemical_equilibrium_system::No_chemical_equilibrium_system()" << endl
         << "Initialising a dummy chemical equilibrium system." << endl;
}

No_chemical_equilibrium_system::~No_chemical_equilibrium_system() {}

int
No_chemical_equilibrium_system::solve_system( double T, double rho, vector<double> &massf )
{
    cout << "No_chemical_equilibrium_system::solve_system()" << endl
         << "This class has no functionality - exiting program." << endl;
         
    exit( FAILURE );
}

int
No_chemical_equilibrium_system::compute_source_terms( vector<double> &massf )
{
    cout << "No_chemical_equilibrium_system::compute_source_terms()" << endl
         << "This class has no functionality - exiting program." << endl;
         
    exit( FAILURE );
}

int
No_chemical_equilibrium_system::f( const vector<double> &y, vector<double> &G )
{
    cout << "No_chemical_equilibrium_system::f()" << endl
         << "This class has no functionality - exiting program." << endl;
         
    exit( FAILURE );
}

int
No_chemical_equilibrium_system::Jac( const vector<double> &y, Valmatrix &dGdy )
{
    cout << "No_chemical_equilibrium_system::Jac()" << endl
         << "This class has no functionality - exiting program." << endl;
         
    exit( FAILURE );
}

int
No_chemical_equilibrium_system::
test_system( double T, double p, std::vector<double> molef )
{
    cout << "No_chemical_equilibrium_system::test_system()" << endl
         << "This class has no functionality - exiting program." << endl;
         
    exit( FAILURE );
}

Partial_equilibrium_reaction *
No_chemical_equilibrium_system::
get_partial_equilibrium_reaction_pointer( size_t index )
{
    cout << "No_chemical_equilibrium_system::get_partial_equilibrium_reaction_pointer()" << endl
         << "This class has no functionality - exiting program." << endl;
         
    exit( FAILURE );
}

