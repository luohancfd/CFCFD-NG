/** \file conservation_systems.cxx
 *  \brief Conservation systems for post-shock relaxation
 *
 *  \author Daniel F Potter
 *  \date 30-Jun-2010
 *
 **/

#include <valarray>

#include "conservation_systems.hh"

using namespace std;

/* FrozenConservationSystem (Rankine-Hugonoit equations) */

FrozenConservationSystem::FrozenConservationSystem() {}

FrozenConservationSystem::FrozenConservationSystem( Gas_model * gm, Gas_data * Q,
                                                    double u )
    :  gmodel_( gm ), Q_(Q)
{
    A_ = Q->rho*u;
    B_ = Q->rho*u*u + Q->p;
    C_ = Q->rho*u*( Q->e[0] + 0.5*u*u  ) + Q->p*u;
   
}

FrozenConservationSystem::FrozenConservationSystem( const FrozenConservationSystem &c )
    : gmodel_( c.gmodel_ ), Q_ ( c.Q_ ), A_( c.A_ ), B_( c.B_ ), C_( c.C_ ) {}

FrozenConservationSystem::~FrozenConservationSystem() {}

void
FrozenConservationSystem::
initialise(Gas_model * gm, Gas_data * Q, double u )
{
    gmodel_ = gm;
    
    Q_ = Q;
    
    A_ = Q->rho*u;
    B_ = Q->rho*u*u + Q->p;
    C_ = Q->rho*u*( Q->e[0] + 0.5*u*u ) + Q->p*u;
}

int FrozenConservationSystem::f( const valarray<double> &y, valarray<double> &G )
{

    double rho, T, u;
    rho = y[0];
    T = y[1];
    u = y[2];
    
    Q_->T[0] = T;
    Q_->rho = rho;
    gmodel_->eval_thermo_state_rhoT( *Q_ );
    
    G[0] = rho*u - A_;
    G[1] = A_*u + Q_->p - B_;
    G[2] = A_*( Q_->e[0] + 0.5*u*u ) + u*Q_->p - C_;
    
    return 0;
}

int FrozenConservationSystem::Jac( const valarray<double> &y, Valmatrix &dGdy )
{
    double rho, T, u;
    rho = y[0];
    T = y[1];
    u = y[2];
    
    // NOTE: f() should have been just called, so Q_ should be correct
    int status;
    double dpdrho = gmodel_->dpdrho_const_T( *Q_, status );
    double dpdT = 1.0 / gmodel_->dTdp_const_rho( *Q_, status );
    double Cv = gmodel_->modal_Cv(*Q_,0);

    dGdy.set(0,0,u);           dGdy.set(0,1,0.0);    dGdy.set(0,2,rho);
    dGdy.set(1,0,dpdrho);      dGdy.set(1,1,dpdT);   dGdy.set(1,2,A_);
    dGdy.set(2,0,u*dpdrho);    dGdy.set(2,1,A_*Cv);   dGdy.set(2,2,A_*u+Q_->p);

    return 0;
}

/* NoneqConservationSystem (multiple species multiple temperatures) */

NoneqConservationSystem::NoneqConservationSystem() {}

NoneqConservationSystem::NoneqConservationSystem( Gas_model * gm, Gas_data * Q,
                                                  double u )
    :  gmodel_( gm ), Q_( Q ), nsp_( gm->get_number_of_species() ), 
       ntm_( gm->get_number_of_modes() )
{
    ndim_ = nsp_ + 1 + ntm_;
    e_index_ = -1;
    for ( int isp=0; isp<nsp_; ++isp ) {
    	if ( gmodel_->species_name(isp)=="e_minus" )
    	    e_index_ = isp;
    }
    A_.resize( ndim_ );
    
    encode_conserved( A_, *Q_, u );
   
}

NoneqConservationSystem::NoneqConservationSystem( const NoneqConservationSystem &c )
    : gmodel_( c.gmodel_ ), Q_( c.Q_ ), nsp_( c.nsp_ ), ntm_( c.ntm_ ),
      ndim_( c.ndim_ ), e_index_( c.e_index_ ), A_( c.A_ ) {}

NoneqConservationSystem::~NoneqConservationSystem() {}

void
NoneqConservationSystem::
initialise(Gas_model * gm, Gas_data * Q, double u)
{
    gmodel_ = gm;
    nsp_ = gmodel_->get_number_of_species();
    ntm_ = gmodel_->get_number_of_modes();    
    ndim_ = nsp_ + 1 + ntm_;
    e_index_ = -1;
    for ( int isp=0; isp<nsp_; ++isp ) {
    	if ( gmodel_->species_name(isp)=="e_minus" )
    	    e_index_ = isp;
    }
    A_.resize( ndim_ );
    
    Q_ = Q;
    
    encode_conserved( A_, *Q_, u );
}

int NoneqConservationSystem::f( const valarray<double> &y, valarray<double> &G )
{
    // 1. Map y vector (solution guess) onto the gas-data structure
    // 1a. Total dendisty and species mass-fractions
    Q_->rho = 0.0;
    for ( int isp=0; isp<nsp_; ++isp )
    	Q_->rho += y[isp];
    for ( int isp=0; isp<nsp_; ++isp )
    	Q_->massf[isp] = y[isp] / Q_->rho;
    // 1b. Modal temperatures
    for ( int itm=0; itm<ntm_; ++itm )
    	Q_->T[itm] = y[nsp_+itm];
    // 1c. Velocity
    double u = y[nsp_+ntm_];
    
    // 2. Fill out the remaining quantities with an EOS evaluation
    gmodel_->eval_thermo_state_rhoT( *Q_ );
    
    // 3. Compute the zero vector elements
    int iG = 0;
    // 3a. Species density flux
    for ( int isp=0; isp<nsp_; ++isp ) {
    	G[iG] = Q_->rho*Q_->massf[isp]*u - A_[iG];
    	++iG;
    }
    // 3b. Total momentum flux
    G[iG] = Q_->rho * u * u + Q_->p - A_[iG];
    ++iG;
    // 3c. Total energy flux
    double E = Q_->e[0] + 0.5 * u * u;
    G[iG] = u * ( Q_->rho * E + Q_->p ) - A_[iG];
    ++iG;
    if ( ntm_ > 1 ) {
	// 3d. Modal energy flux (but not first or last mode)
	for ( int itm=1; itm<(ntm_-1); ++itm ) {
	    G[iG] = u * ( Q_->rho * Q_->e[itm] ) - A_[iG];
	    ++iG;
	}
	// 3e. Modal energy flux for last mode (which is assumed to govern electrons)
	G[iG] = u * ( Q_->rho * Q_->e[ntm_-1] + Q_->p_e ) - A_[iG];
    }
    
    // print_valarray(G);
    
    return 0;
}

int NoneqConservationSystem::Jac( const valarray<double> &y, Valmatrix &dGdy )
{
    // NOTE: f() should have been just called, so *Q_ should be correct
    double u = y[nsp_+ntm_];
    
    int iG = 0;
    int status;

    // 1. Species mass flux rows
    for ( int isp=0; isp<nsp_; ++isp ) {
    	// 1a. deriv. wrt rho_j
    	for ( int jsp=0; jsp<nsp_; ++jsp ) {
    	    if ( isp==jsp ) dGdy.set(iG,jsp,u);
    	    else dGdy.set(iG,jsp,0.0);
    	}
    	// 1b. deriv. wrt T_i
    	for ( int itm=0; itm<ntm_; ++itm )
    	    dGdy.set(iG,nsp_+itm,0.0);
    	// 1c. deriv. wrt u
    	dGdy.set(iG,nsp_+ntm_,Q_->massf[isp]*Q_->rho);
    	++iG;
    }
    
    // 2. Total momentum flux row
    // 2a. deriv. wrt rho_i
    for ( int isp=0; isp<nsp_; ++isp )
    	dGdy.set(iG,isp,u*u+gmodel_->dpdrho_i_const_T( *Q_, isp, status ));
    // 2b. deriv. wrt T_i
    for ( int itm=0; itm<ntm_; ++itm )
    	dGdy.set(iG,nsp_+itm,gmodel_->dpdT_i_const_rho( *Q_, itm, status ));
    // 2c. deriv wrt u
    dGdy.set(iG,nsp_+ntm_,2.0*Q_->rho*u);
    ++iG;
    
    // 3. Total energy flux row
    double E = Q_->e[0] + 0.5 * u * u;
    // 3a. deriv. wrt rho_i
    for ( int isp=0; isp<nsp_; ++isp )
    	dGdy.set(iG,isp,u*E+u*gmodel_->dpdrho_i_const_T( *Q_, isp, status ));
    // 3b. deriv. wrt T_i
    for ( int itm=0; itm<ntm_; ++itm )
    	dGdy.set(iG,nsp_+itm,Q_->rho*u*gmodel_->modal_Cv( *Q_, itm ) +
    	    u*gmodel_->dpdT_i_const_rho( *Q_, itm, status ));
    // 3c. deriv wrt u
    dGdy.set(iG,nsp_+ntm_,Q_->rho*Q_->e[0]+Q_->rho*1.5*u*u + Q_->p);
    ++iG;
    
    if ( ntm_ > 1 ) {
	// 4. Modal energy flux rows (but not first or last)
	for ( int itm=1; itm<(ntm_-1); ++itm ) {
	    // 4a. deriv. wrt rho_i
	    for ( int isp=0; isp<nsp_; ++isp )
		dGdy.set(iG,isp,u*Q_->e[itm]);
	    // 4b. deriv. wrt T_j
	    for ( int jtm=0; jtm<ntm_; ++jtm ) {
		if ( itm==jtm ) dGdy.set(iG,nsp_+jtm,
		                      Q_->rho*u*gmodel_->modal_Cv( *Q_, itm ));
		else dGdy.set(iG,nsp_+jtm,0.0);
	    }
	    // 4c. deriv wrt u
	    dGdy.set(iG,nsp_+ntm_,Q_->rho*Q_->e[itm]);
	    ++iG;
	}
	
	// 5. Last modal energy flux row (assumed to be govern electrons)
	int itm = ntm_ - 1;
	// 5a. deriv. wrt rho_i
	for ( int isp=0; isp<nsp_; ++isp ) {
	    if ( isp==e_index_ )
	    	dGdy.set(iG,isp,u*Q_->e[itm]+
	    	    u*gmodel_->dpdrho_i_const_T( *Q_, isp, status ));
	    else dGdy.set(iG,isp,u*Q_->e[itm]);
	}
	// 5b. deriv. wrt T_j
	for ( int jtm=0; jtm<ntm_; ++jtm ) {
	    if ( itm==jtm )
	    	dGdy.set(iG,nsp_+jtm,Q_->rho*u*gmodel_->modal_Cv( *Q_, itm ) + 
	    	    u*gmodel_->dpdT_i_const_rho( *Q_, itm, status ));
	    else dGdy.set(iG,nsp_+jtm,0.0);
	}
	// 5c. deriv wrt u
	dGdy.set(iG,nsp_+ntm_,Q_->rho*Q_->e[itm]+Q_->p_e);
    }
    
    // cout << "dGdy = " << dGdy.str() << endl;

    return 0;
}

int NoneqConservationSystem::encode_conserved( valarray<double> &y, 
                                               const Gas_data &Q,
                                               const double u )
{
    int iy=0;
    // 1a. Species density flux
    for ( int isp=0; isp<nsp_; ++isp ) {
    	y[iy] = Q.massf[isp] * Q.rho * u;
    	++iy;
    }
    // 1b. Total momentum flux
    y[iy] = Q.rho * u * u + Q.p;
    ++iy;
    // 1c. Total energy flux
    double E = Q.e[0] + 0.5 * u * u;
    y[iy] = u * ( Q.rho * E + Q.p );
    ++iy;
    if ( ntm_ > 1 ) {
	// 1d. Modal energy flux (but not first or last mode)
	for ( int itm=1; itm<(ntm_-1); ++itm ) {
	    y[iy] = u * ( Q.rho * Q.e[itm] );
	    ++iy;
	}
	// 1e. Modal energy flux for last mode (which is assumed to govern electrons)
	y[iy] = u * ( Q.rho * Q.e[ntm_-1] + Q.p_e );
    }
    
    return 0;
}

int NoneqConservationSystem::set_constants( const valarray<double> &A )
{
    for ( size_t i=0; i<A.size(); ++i )
    	A_[i] = A[i];
    
    return 0;
}

