/** \file cr_rr_coeffs.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 05-04-10: Improved port from old lib/radiation
 *  \brief Declarations for collisional-radiative reaction rate coefficients
 *
 **/
 
#ifndef CR_RR_COEFFS_HH
#define CR_RR_COEFFS_HH
 
#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "radiator.hh"
 
#include "../../gas/models/gas_data.hh"
#include "../../geometry2/source/gpath.hh"

class CR_ReactionRateCoefficient {
public:
    virtual ~CR_ReactionRateCoefficient() {};
public:
    virtual double get_rate( double T, Gas_data &Q ) = 0;
    bool get_equilibrium_flag() { return equilibrium_flag; };
    std::string get_type_str()
    { return type; }
    virtual std::string get_latex_string() = 0;

protected:
    std::string type;
    bool equilibrium_flag;
    
};

class ZeroRate : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    ZeroRate();
    
    /// \brief Destructor
    ~ZeroRate() {};
    
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
};

class FromEquilibriumConstant : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    FromEquilibriumConstant();
    
    /// \brief Destructor
    ~FromEquilibriumConstant() {};
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
};

class RadGeneralisedArrhenius : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    RadGeneralisedArrhenius( lua_State *L );
    
    /// \brief Explicit constructor
    RadGeneralisedArrhenius( double A, double n, double T_a );
    
    /// \brief Destructor
    ~RadGeneralisedArrhenius() {};
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string();
    
private:
    double A;
    double n;
    double T_a;
};

class RadGeneralisedArrheniusPark : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    RadGeneralisedArrheniusPark( lua_State *L );
    
    /// \brief Explicit constructor
    RadGeneralisedArrheniusPark( double A, double n, double T_a );
    
    /// \brief Destructor
    ~RadGeneralisedArrheniusPark() {};
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string();
    
private:
    double A;
    double n;
    double T_a;
};

class DrawinOpticallyAllowedElectronImpactExcitation : public CR_ReactionRateCoefficient {
public:
    /// \brief Simple constructor
    DrawinOpticallyAllowedElectronImpactExcitation( double E_l, double E_u );
    
    /// \brief Destructor
    ~DrawinOpticallyAllowedElectronImpactExcitation() {};
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
private:
    double E_l;
    double E_u;
};

class DrawinOpticallyForbiddenElectronImpactExcitation : public CR_ReactionRateCoefficient {
public:
    /// \brief Simple constructor
    DrawinOpticallyForbiddenElectronImpactExcitation( double E_l, double E_u );
    
    /// \brief Destructor
    ~DrawinOpticallyForbiddenElectronImpactExcitation() {};
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
    
private:
    double E_l;
    double E_u;
};

class GryzinskiElectronImpactExcitation : public CR_ReactionRateCoefficient {
public:
    /// \brief Simple constructor
    GryzinskiElectronImpactExcitation( Radiator * rad, ElecLev * elev_l, ElecLev * elev_u );
    
    /// \brief Even simpler constructor
    GryzinskiElectronImpactExcitation( double E_l, double E_u, double E_up1, double I );
    
    /// \brief Destructor
    ~GryzinskiElectronImpactExcitation() {};
    
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
    
private:
    double sigma( double E, double E_i, double E_j );
    double compute_integral( double T, double E_i, double E_j  );
    double r( double x, double E_i, double E_j );
    double drdx( double x, double E_i, double E_j );
    
private:
    double E_l;
    double E_u, E_up1;
    double I;
};

class DrawinElectronImpactIonization : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    DrawinElectronImpactIonization( double E_l, double I );
    
    /// \brief Destructor
    ~DrawinElectronImpactIonization() {};
    
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
    
private:
    double E_l;
    double I;
};

class CJDrawinElectronImpactIonization : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    CJDrawinElectronImpactIonization( double E_l, double I );
    
    /// \brief Destructor
    ~CJDrawinElectronImpactIonization() {};
    
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
    
private:
    double E_l;
    double I;
};

class OpticallyThinExponentialDecay : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor from lua file (for diatoms)
    OpticallyThinExponentialDecay( lua_State * L );
    
    /// \brief Constructor from parameters (for atoms)
    OpticallyThinExponentialDecay( double tau );
    
    /// \brief Destructor
    ~OpticallyThinExponentialDecay() {};
    
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string();
    
private:
    double tau;
    double lambda;
};

class OpticallyVariableExponentialDecay : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor from lua file (for diatoms)
    OpticallyVariableExponentialDecay( lua_State * L );

    /// \brief Constructor
    OpticallyVariableExponentialDecay( double tau, double wavel, double wavel_switch=200.0, double lambda_lower=0.0, double lambda_upper=1.0 );
    
    /// \brief Destructor
    ~OpticallyVariableExponentialDecay() {};
    
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
    
private:
    double tau;
    double lambda;
};

class CurveFitExponentialDecay : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    CurveFitExponentialDecay( double tau );
    
    /// \brief Destructor
    ~CurveFitExponentialDecay() {};
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
    
private:
    double tau;
};

class FrostNitrogenElectronImpactExcitation : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    FrostNitrogenElectronImpactExcitation( lua_State * L, ElecLev * elev_l, ElecLev * elev_u );
    
    /// \brief Constructor from file
    FrostNitrogenElectronImpactExcitation( std::string fname, ElecLev * elev_l, ElecLev * elev_u );
    
    /// \brief Destructor
    ~FrostNitrogenElectronImpactExcitation();
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
    
private:
    Spline * gamma_fit;
    double delta_E;
    int g_l;
};

class ZatsarinnyTayalOxygenElectronImpactExcitation : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    ZatsarinnyTayalOxygenElectronImpactExcitation( lua_State * L, ElecLev * elev_l, ElecLev * elev_u );
    
    /// \brief Constructor from file
    ZatsarinnyTayalOxygenElectronImpactExcitation( std::string fname, ElecLev * elev_l, ElecLev * elev_u );
    
    /// \brief Destructor
    ~ZatsarinnyTayalOxygenElectronImpactExcitation();
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
    
private:
    Spline * gamma_fit;
    double delta_E;
    int g_l;
};

class SunoKatoCarbonElectronImpactExcitation : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    SunoKatoCarbonElectronImpactExcitation( lua_State * L, ElecLev * elev_l, ElecLev * elev_u );
    
    /// \brief Constructor from file
    SunoKatoCarbonElectronImpactExcitation( std::string fname, ElecLev * elev_l, ElecLev * elev_u );
    
    /// \brief Destructor
    ~SunoKatoCarbonElectronImpactExcitation();
    
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
    
    double compute_gamma( double y );
    
    double eval_t_from_x( double x, double y );
    
    double dtdx( double t, double y );
    
    double eval_X_from_x( double x, double y );
    
    double dXdx( double X, double y );
    
private:
    double delta_E;
    int g_l;
    
    double V_if;
    int eqn;
    double A, B, C, D, E, F;
};

class SunoKatoCarbonElectronImpactIonization : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    SunoKatoCarbonElectronImpactIonization( lua_State * L, ElecLev * elev, double I );
    
    /// \brief Constructor from file
    SunoKatoCarbonElectronImpactIonization( std::string fname, ElecLev * elev, double I );
    
    /// \brief Destructor
    ~SunoKatoCarbonElectronImpactIonization() {};
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
    
    double compute_omega( double X );
    
    double compute_gamma( double y );
    
    double eval_t_from_x( double x, double y );
    
    double dtdx( double t, double y );
    
    double eval_X_from_x( double x, double y );
    
    double dXdx( double X, double y );
private:
    double delta_E;
    int g_l;
    
    double V_if;
    int eqn;
    double A1, A2, A3, A4, A5;
};

class KuncSoonElectronImpactIonization : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    KuncSoonElectronImpactIonization( lua_State * L, ElecLev * elev, double I );
    
    /// \brief Constructor from file
    // KuncSoonElectronImpactIonization( std::string fname, ElecLev * elev, double I );
    
    /// \brief Constructor from data
    KuncSoonElectronImpactIonization( double A, double chi, double Q, ElecLev * elev, double I );
    
    /// \brief Destructor
    ~KuncSoonElectronImpactIonization() {};
    
public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }
    
private:
    double E_l;
    double I;
    double A;
    double chi;
    double Q_val;
    double l;
};

class OpticallyThinPhotoRecombination : public CR_ReactionRateCoefficient {
public:
    /// \brief Constructor
    OpticallyThinPhotoRecombination( ElecLev * elev, double I );

    /// \brief Destructor
    ~OpticallyThinPhotoRecombination() {};

public:
    double get_rate( double T, Gas_data &Q );
    std::string get_latex_string() { return ""; }

private:
    double eval_integrand( double T, double eps);

    double eval_integral( double T );

private:
    ElecLev * elev;
    double E_l;
    double I_;
    int g_ion;
    std::vector<double> nus;
    int nnus;
};

CR_ReactionRateCoefficient * create_explicit_rate_coeff( lua_State * L, std::string parameter_field );

#endif

