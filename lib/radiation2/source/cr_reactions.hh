/** \file cr_reactions.hh
 *  \ingroup radiation2
 *
 *  \author Daniel F. Potter
 *  \version 03-04-10 (port from lib/radiation)
 *  \brief Declarations for nonequilibirum population model transitions
 *
 **/
 
#ifndef CR_REACTIONS_HH
#define CR_REACTIONS_HH

#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../gas/models/gas_data.hh"
#include "../../nm/source/no_fuss_linear_algebra.hh"

#include "radiator.hh"
#include "cr_rr_coeffs.hh"

class CR_Reaction {
public:
    /// \brief Constructor
    CR_Reaction( std::string type );
    
    /// \brief Deconstructor
    virtual ~CR_Reaction();
    
public:
    std::string get_type();
    std::string get_equation();
    std::string get_forward_rate_coeff_type();
    
public:
    virtual int eval_reaction_rates( double T_f, double T_b, Gas_data &Q, double &k_f, double &k_b );
    virtual double eval_equilibrium_constant( double T ) = 0;
    virtual int add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy ) = 0;
    virtual int add_eval_contributions( Gas_data &Q, std::valarray<double> &G ) = 0;
    virtual int add_source_vector_contributions( Gas_data &Q, std::valarray<double> &C ) = 0;
    virtual std::string get_latex_string() = 0;
    
protected:
    std::string type;
    std::string equation;
    CR_ReactionRateCoefficient *forward_rate_coeff;
    CR_ReactionRateCoefficient *backward_rate_coeff;
};

class HeavyParticleImpactExcitation : public CR_Reaction {
public:
    /// \Brief Constructor from lua input file (for diatomic species)
    HeavyParticleImpactExcitation( lua_State * L, Radiator * rad );
    
    ~HeavyParticleImpactExcitation() {};
    
public:
    void set_heavy_particle_pointer( Radiator * M_pointer )
    { M = M_pointer; }
    
    double eval_equilibrium_constant( double T );
    int add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy );
    int add_eval_contributions( Gas_data &Q, std::valarray<double> &G );
    int add_source_vector_contributions( Gas_data &Q, std::valarray<double> &C );
    std::string get_latex_string();
    
public:
    std::string M_name;
    
private:
    int iTv;
    int iT;
    Radiator * rad;
    NoneqElecLev * ne_elev_l;
    NoneqElecLev * ne_elev_u;
    Radiator * M;
};

class ElectronImpactExcitation : public CR_Reaction {
public:
    /// \Brief Constructor from lua input file (for diatomic species)
    ElectronImpactExcitation( lua_State * L, Radiator * rad );
    
    /// \brief Constructor from model type (for atomic species)
    ElectronImpactExcitation( lua_State * L, std::string model, Radiator * rad, NoneqElecLev * ne_elev_l, NoneqElecLev * ne_elev_u );
    
    ~ElectronImpactExcitation() {};
    
public:
    void set_electron_pointer( Radiator * elec_pointer )
    { elec = elec_pointer; }
    
    double eval_equilibrium_constant( double T );
    int add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy );
    int add_eval_contributions( Gas_data &Q, std::valarray<double> &G );
    int add_source_vector_contributions( Gas_data &Q, std::valarray<double> &C );
    std::string get_latex_string();
   
private:
    int iTe;
    int iT;
    NoneqElecLev * ne_elev_l;
    NoneqElecLev * ne_elev_u;
    Radiator * elec;
    std::string rad_name;
};

class HeavyParticleImpactDissociation : public CR_Reaction {
public:
    /// \Brief Constructor from lua input file (for diatomic species)
    HeavyParticleImpactDissociation( lua_State * L, Radiator * rad );
    
    ~HeavyParticleImpactDissociation() {};
    
public:
    void set_atom_A_pointer( Radiator * atom_pointer )
    { atom_A = atom_pointer; }
    void set_atom_B_pointer( Radiator * atom_pointer )
    { atom_B = atom_pointer; }
    void set_heavy_particle_pointer( Radiator * M_pointer )
    { M = M_pointer; }
    
    double eval_equilibrium_constant( double T );
    int add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy );
    int add_eval_contributions( Gas_data &Q, std::valarray<double> &G );
    int add_source_vector_contributions( Gas_data &Q, std::valarray<double> &C );
    std::string get_latex_string();
    
public:
    std::string M_name;
    std::string atom_A_name;
    std::string atom_B_name;
    
private:
    Radiator * rad;
    double D;
    int iTv;
    int iT;
    NoneqElecLev * ne_elev_l;
    Radiator * atom_A;
    Radiator * atom_B;
    Radiator * M;
};

class ElectronImpactDissociation : public CR_Reaction {
public:
    /// \Brief Constructor from lua input file (for diatomic species)
    ElectronImpactDissociation( lua_State * L, Radiator * rad );
    
    ~ElectronImpactDissociation() {};
    
public:
    void set_atom_A_pointer( Radiator * atom_pointer )
    { atom_A = atom_pointer; }
    void set_atom_B_pointer( Radiator * atom_pointer )
    { atom_B = atom_pointer; }
    void set_electron_pointer( Radiator * elec_pointer )
    { elec = elec_pointer; }
    
    double eval_equilibrium_constant( double T );
    int add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy );
    int add_eval_contributions( Gas_data &Q, std::valarray<double> &G );
    int add_source_vector_contributions( Gas_data &Q, std::valarray<double> &C );
    std::string get_latex_string();
    
public:
    std::string atom_A_name;
    std::string atom_B_name;
    
private:
    Radiator * rad;
    double D;
    int iTv;
    int iTe;
    int iT;
    NoneqElecLev * ne_elev_l;
    Radiator * atom_A;
    Radiator * atom_B;
    Radiator * elec;
    std::string rad_name;
};

class ElectronImpactIonization : public CR_Reaction {
public:
    /// \Brief Constructor from lua input file (for diatomic species)
    ElectronImpactIonization( lua_State * L, Radiator * rad );
    
    /// \brief Constructor from model type (for atomic species)
    ElectronImpactIonization( lua_State * L, std::string model, Radiator * rad, NoneqElecLev * ne_elev_l );
    
    ~ElectronImpactIonization() {};
    
public:
    void set_ion_pointer( Radiator * ion_pointer )
    { ion = ion_pointer; }
    void set_electron_pointer( Radiator * elec_pointer )
    { elec = elec_pointer; }
    
    double eval_equilibrium_constant( double T );
    int add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy );
    int add_eval_contributions( Gas_data &Q, std::valarray<double> &G );
    int add_source_vector_contributions( Gas_data &Q, std::valarray<double> &C );
    std::string get_latex_string();
    
private:
    Radiator * rad;
    int iTe;
    int iT;
    NoneqElecLev * ne_elev_l;
    Radiator * ion;
    Radiator * elec;
};

class RadiativeTransition : public CR_Reaction {
public:
    /// \Brief Constructor from lua input file (for diatomic species)
    RadiativeTransition( lua_State * L, Radiator * rad );
    
    /// \brief Constructor from model type (for atomic species)
    RadiativeTransition( std::string model, NoneqElecLev * ne_elev_u, NoneqElecLev * ne_elev_l, double A_ul );
    
    ~RadiativeTransition() {};
    
public:
    int eval_reaction_rates( double T_f, double T_b, Gas_data &Q, double &k_f, double &k_b );
    double eval_equilibrium_constant( double T );
    int add_jacobian_contributions( Gas_data &Q, Valmatrix &dGdy );
    int add_eval_contributions( Gas_data &Q, std::valarray<double> &G );
    int add_source_vector_contributions( Gas_data &Q, std::valarray<double> &C );
    std::string get_latex_string();
    
private:
    NoneqElecLev * ne_elev_u;
    NoneqElecLev * ne_elev_l;
};

void tokenize_equation_string( std::string equation, std::vector<std::string> &tks );

void get_latex_species_piecies( std::string name, std::string &lname, std::string &lev_prefix );

#endif

