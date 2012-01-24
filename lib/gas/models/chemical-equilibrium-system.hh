// Author: Daniel F. Potter
// Date: 06-Oct-2009
// Place: Brisbane, Queendland, AUST

#ifndef CHEMICAL_EQUILIBRIUM_SYSTEM_HH
#define CHEMICAL_EQUILIBRIUM_SYSTEM_HH
        
#include <vector>
#include <string>
#include <valarray>

#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"
#include "../../../lib/nm/source/zero_system.hh"
#include "../../../lib/nm/source/zero_finders.hh"
#include "chemical-species.hh"

class Partial_equilibrium_reaction {
public:
    Partial_equilibrium_reaction( std::vector<int> betas, 
    	 			  std::vector<Chemical_species*> species );
    ~Partial_equilibrium_reaction() {};
    
    std::string str();
    
    double eval_K_p( double T )
    { return s_eval_K_p(T); }
    
    double eval_K_c( double T )
    { return s_eval_K_c(T); }
    
    // NOTE: these are actually not the equilibrium constants, but a pseudo-eq-constant...
    double eval_equilibrium_constant( double T )
    { return s_eval_equilibrium_constant(T); }
    
    void store_equilibrium_constant( double T )
    { K_m_ = s_eval_equilibrium_constant(T); }
    
    double get_K_m()
    { return K_m_; }
    
    int get_nsp()
    { return (int) species_.size(); }
    
    int get_G( int irsp )
    { return - betas_[irsp]; }
    
    int G_sum();
    
    int get_beta( int irsp )
    { return betas_[irsp]; }
    
    int get_isp( int irsp )
    { return species_[irsp]->get_isp(); }
    
    int get_nrsp()
    { return (int) species_.size(); }
    
private:
    std::vector<int> betas_;
    std::vector<Chemical_species*> species_;
    double K_m_;
    
    double s_eval_K_p( double T );
    
    double s_eval_K_c( double T );
    
    double s_eval_equilibrium_constant( double T );
};

class Chemical_equilibrium_system : public ZeroSystem {
public:
    Chemical_equilibrium_system() {}
    Chemical_equilibrium_system( double min_massf, std::vector<Chemical_species*> &species );
    ~Chemical_equilibrium_system();
    
public:
    virtual int solve_system( double T, double rho, std::vector<double> &massf );
    virtual int test_system( double T, double p, std::vector<double> molef );
    
    virtual int f( const std::valarray<double> &y, std::valarray<double> &G );
    virtual int Jac( const std::valarray<double> &y, Valmatrix &dGdy );
    
    virtual Partial_equilibrium_reaction * get_partial_equilibrium_reaction_pointer( size_t index )
    { return pe_reactions_[index]; }
    
private:
    std::vector<Partial_equilibrium_reaction*> pe_reactions_;
    std::vector<int> be_indices_;
    std::vector<double> charge_weightings_;
    matrix be_weightings_;
    size_t nsp_;
    bool ions_;
    double log_rho_;
    double min_massf_;
    
    NewtonRaphsonZF zero_solver_;
    double f_jac_;
    
    std::valarray<double> Q_;
    std::valarray<double> yguess_;
    std::valarray<double> yout_;
    
    virtual int compute_source_terms( std::vector<double> &massf );
};

class No_chemical_equilibrium_system : public Chemical_equilibrium_system {
public:
    No_chemical_equilibrium_system( double min_massf, std::vector<Chemical_species*> &species );
    ~No_chemical_equilibrium_system();
    
public:
    virtual int solve_system( double T, double rho, std::vector<double> &massf );
    virtual int test_system( double T, double p, std::vector<double> molef );
    
    virtual int f( const std::valarray<double> &y, std::valarray<double> &G );
    virtual int Jac( const std::valarray<double> &y, Valmatrix &dGdy );
    
    virtual Partial_equilibrium_reaction * get_partial_equilibrium_reaction_pointer( size_t index );

private:
    int compute_source_terms( std::vector<double> &massf );

};

#endif
