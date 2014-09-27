/*  \file reaction_scheme.hh
 *  \brief Declarations for the ReactionScheme class
 *
 *  \author Rowan J Gollan
 *  \version 23-Feb-2006
 **/

#ifndef R_SCHEME_HH
#define R_SCHEME_HH

#include <string>
#include <vector>

#include "../../nm/source/no_fuss_linear_algebra.hh"
#include "../../nm/source/ode_system.hh"
#include "../../nm/source/ode_solver.hh"
#include "../models/gas_data.hh"
#include "../models/thermally_perfect_gas_mix.hh"
#include "reaction.hh"
#include "reacting_species.hh"
#include "reaction_pieces.hh"
#include "species_pieces.hh"
#include "e_vib_coupling.hh"

/** \brief Base abstract class for a ReactionScheme

\author Rowan J Gollan
\version 24-Feb-2006

This abstract class for a ReactionScheme encapsulates
the basic behaviour of a ReactionScheme.
An object of this type provides a direct link to 
a gas dynamic flow solver.
Given an input gas state and a flow timestep, this object
can update the gas state at the end of that timestep due
to chemical reactions.
At this level of abstraction, we do not care how the
object does the update.
Two candidate ideas are solving the ODE system with an integrator
or some form of Richardson extrapolation to compute the new
mass fractions (and concentrations).

**/

class ReactionScheme {
public:
    /// \brief Normal constructor
    ReactionScheme( const std::string name, const vector<Reaction*> reactions,
		    int nsp, double err_tol, double T_trigger );
    ReactionScheme( const std::string name, int nsp,
		    double err_tol, double T_trigger );
    
    ReactionScheme( const std::string name, const std::string input_file );
    /// \brief Copy constructor
    ReactionScheme( const ReactionScheme &r );

    /// \brief Default destructor
    virtual ~ReactionScheme();

    /// \brief clone() function
    ///
    /// \author Rowan J Gollan
    /// \version 24-Feb-2006
    ///
    virtual ReactionScheme* clone() = 0;

    /// \brief String representation
    virtual std::string str() const = 0;

    Reaction* get_reaction_pointer(int ir) { std::cout << "Trying to grab reaction: " << ir << std::endl;  return reactions_[ir]; }
    ReactingSpecies* get_reacting_species_pointer(int isp) { return reacting_species_[isp]; }
    int get_number_of_reactions() { return int( reactions_.size() ); }

    // ---- Basic behaviour of a ReactionScheme ---- //
    virtual int update_gas_state( Gas_data &Q, double dt_flow, bool include_evib_exchange=true ) = 0;

    void add_reaction( Reaction *r );
    virtual void update_reacting_species();

    int get_nu_for_reaction(int ir, int index)
    { return reactions_[ir]->get_nu(index); }

    double get_k_f_for_reaction(int ir, Gas_data &Q);

    ThermallyPerfectGasMix *get_tpgm_pointer()
    { return pgm_; }

protected:
    GasModel *g_;
    ThermallyPerfectGasMix *pgm_;
    std::string name_;
    double err_tol_;
    double T_trigger_;
    std::vector<Reaction*> reactions_;
    std::vector<ReactingSpecies*> reacting_species_;
};


///// This class uses multiple inheritance because
///// I want to share the data _reactions.
///// I want the ReactionScheme to be able to use it,
///// but the eval function of an OdeSystem
///// also needs to access _reactions.
///// I could make a special OdeSystem and have it
///// point to its parent but I don't want to manage
///// the pointers to reactions in two places.

class ReactionSchemeODE : public ReactionScheme, public OdeSystem {
public:
    /// \brief Normal constructor
    ReactionSchemeODE( const std::string name, const std::vector<Reaction*> reactions,
		       int nsp, double err_tol, double T_trigger,
		       OdeSolver *ode_solver );

    ReactionSchemeODE( const std::string name,
		       int nsp, double err_tol, double T_trigger,
		       OdeSolver *ode_solver );

    /// \brief Copy constructor
    ReactionSchemeODE( const ReactionSchemeODE &r );

    /// \brief Default destructor
    virtual ~ReactionSchemeODE();

    virtual ReactionSchemeODE* clone();

    /// \brief string representation
    virtual std::string str() const;
    
    // -------- ReactionScheme behaviour ---------- //

    /// \brief Update the gas state due to chemical reactions (using an ODE solver)
    int update_gas_state( Gas_data &Q, double dt_flow, bool include_evib_exchange=true );

    // -------- OdeSystem behaviour ----------- //
    
    /// \brief Compute the rate of species production.
    int eval( const std::vector<double> &y, std::vector<double> &ydot );
    
    /// \brief Compute the rate of species production (split)
    int eval_split( const std::vector<double> &y, 
		    std::vector<double> &q, std::vector<double> &L );

    /// \brief Select an appropriate timestep for system integration.
    double stepsize_select( const std::vector<double> &y );

    /// \brief Check the mass fractions for consistency.
    bool passes_system_test( std::vector<double> &y );
    
    void update_reacting_species();

    void set_gas_data_ptr( Gas_data *Q ) { Q_ = Q; }

    void evolution_rates(int isp, std::vector<double> &W, Gas_data &Q );

protected:
    Gas_data *Q_;
    OdeSolver *ode_solver_;
    std::vector<double> yin_;
    std::vector<double> yout_;
    std::vector<double> w_f_;
    std::vector<double> w_b_;
    std::vector<double> q_;
    std::vector<double> L_;
    std::vector<double> mf_;
    std::vector<double> c_;
		       
};

// -----------------------------------------------------------
// Functions that provide access to the "managed" finite-rate
// chemistry module.

/// \brief Selects and sets up the managed ReactionScheme
int set_reaction_scheme( const std::string name,
			 const std::string input_file );

void clear_reaction_scheme_pointer();
ReactionScheme *get_reaction_scheme_pointer();
ThermallyPerfectGasMix *get_tpgm_pointer();

int estimate_appropriate_subcycle(double dt_flow, double dt_chem, double *dt, int *no_substeps);
int perform_chemical_increment(Gas_data &Q, double t_interval);
int get_nu_for_reaction(int ir, int index);
double get_k_f_for_reaction(int ir, Gas_data &Q);

ReactionSchemeODE* initialize_ReactionSchemeODE( const std::string name,
						 const std::string input_file);


class ReactionSchemeODE_MC : public ReactionScheme, public OdeSystem {
public:
    ReactionSchemeODE_MC( const std::string name,  const std::vector<Reaction*> reactions,
			  int nsp, double err_tol, double T_trigger,
			  OdeSolver *ode_solver );

    ReactionSchemeODE_MC( const ReactionSchemeODE_MC &r );

     /// \brief Default destructor
    virtual ~ReactionSchemeODE_MC();

    virtual ReactionSchemeODE_MC* clone();

    virtual std::string str() const;

    int initialise_evib_coupling(ConfigParser &cfg);
    // ------------ ReactionScheme behaviour -------------- //
    int update_gas_state( Gas_data &Q, double dt_flow, bool include_evib_exchange=true );

    // ------------ OdeSystem behaviour -------------- //
    
    int eval( const std::vector<double> &y, std::vector<double> &ydot );
    
    double stepsize_select( const std::vector<double> &y );
    bool passes_system_test( std::vector<double> &y );

    void set_init_conc(int isp, double conc)
    { spec_[isp]->set_init_conc(conc); }

    double eval_conc(int isp, std::vector<double> &w)
    { return spec_[isp]->eval_conc(w); }
    void set_gas_data_ptr( Gas_data *Q ) { Q_ = Q; }

    int update_evib(Gas_data &Q, std::vector<double> &y,
		    std::vector<double> &c_old, std::vector<double> &c_new);

private:
    Gas_data *Q_;
    OdeSolver *ode_solver_;
    std::vector<double> yin_;
    std::vector<double> yout_;
    std::vector<double> w_;
    std::vector<double> c_;
    std::vector<double> cinit_;

    std::vector<SpeciesPieces*> spec_;

    std::vector<Evib_coupling*> evc_;
    

};
#endif
