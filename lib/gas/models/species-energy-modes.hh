// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#ifndef SPECIES_ERNEGY_MODES_HH
#define SPECIES_ERNEGY_MODES_HH

#define EVAL_ENTROPY_AT_1ATM 1

#include <string>

#include "gas_data.hh"
#include "physical_constants.hh"
#include "diatom-electronic-level.hh"

class Species_energy_mode {
public:
    Species_energy_mode( int isp=-1, double R=0.0, double min_massf_=1.0e-10,
    	                 std::string type="none", int iT=-1 );
    virtual ~Species_energy_mode() {}
    
    void set_iT( int iT ) { iT_ = iT; }
    
    int get_iT() { return iT_; }
    
    double eval_weighted_energy( const Gas_data &Q )
    { return Q.massf[isp_]*s_eval_energy( Q.T[iT_] ); }
//    { return ( Q.massf[isp_]>min_massf_ ) ? Q.massf[isp_]*s_eval_energy( Q.T[iT_] ) : 0.0; }
    
    double eval_energy( const Gas_data &Q )
    { return s_eval_energy( Q.T[iT_] ); }
    
    double eval_energy_from_T( double T )
    { return s_eval_energy( T ); }
    
    double eval_weighted_enthalpy( const Gas_data &Q )
    { return Q.massf[isp_]*s_eval_enthalpy( Q.T[iT_] ); }
    
    double eval_enthalpy( const Gas_data &Q )
    { return s_eval_enthalpy( Q.T[iT_] ); }
    
    double eval_enthalpy_from_T( double T )
    { return s_eval_enthalpy( T ); }
    
    double eval_weighted_entropy( const Gas_data &Q )
    { return Q.massf[isp_]*s_eval_entropy(Q.T[iT_],Q.massf[isp_]*Q.rho*R_*Q.T[iT_]); }
    
    double eval_entropy( const Gas_data &Q )
#   if EVAL_ENTROPY_AT_1ATM
    { return s_eval_entropy(Q.T[iT_],PC_P_atm); }
#   else
    { return s_eval_entropy(Q.T[iT_],Q.massf[isp_]*Q.rho*R_*Q.T[iT_]); }
#   endif
    
    // NOTE: assuming we want standard state entropy, therefore p=1atm
    double eval_entropy_from_T( double T )
    { return s_eval_entropy(T,PC_P_atm); }
    
    double eval_weighted_Cv( Gas_data &Q )
    { return Q.massf[isp_]*s_eval_Cv( Q.T[iT_] ); }
    
    double eval_Cv( const Gas_data &Q )
    { return s_eval_Cv( Q.T[iT_] ); }
    
    double eval_Cv_from_T( double T=0.0 )
    { return s_eval_Cv( T ); }
    
    double eval_weighted_Cp( Gas_data &Q )
    { return Q.massf[isp_]*s_eval_Cp( Q.T[iT_] ); }
    
    double eval_Cp( const Gas_data &Q )
    { return s_eval_Cp( Q.T[iT_] ); }
    
    double eval_Cp_from_T( double T=0.0 )
    { return s_eval_Cp( T ); }
    
    double eval_Q_from_T( double T=0.0, double A=-1.0 )
    { return s_eval_Q(T,A); }
    
    std::string get_type()
    { return type_; }
    
protected:
    int isp_;
    double R_;
    double min_massf_;
    std::string type_;
    int iT_;
    
    virtual double s_eval_energy( double T ) = 0;
    virtual double s_eval_enthalpy( double T ) = 0;
    virtual double s_eval_entropy( double T, double p ) = 0;
    virtual double s_eval_Cv( double T ) = 0;
    virtual double s_eval_Cp( double T ) = 0;
    virtual double s_eval_Q( double T, double A ) = 0;
    
};

class One_level_electronic : public Species_energy_mode {
public:
    One_level_electronic( int isp, double R, double min_massf, int g, double theta );
    ~One_level_electronic() {}
    
private:
    int g_;
    double theta_;
    
    double s_eval_energy( double T );
    double s_eval_enthalpy( double T ) { return s_eval_energy(T); }
    double s_eval_entropy( double T, double p );
    double s_eval_Cv( double T );
    double s_eval_Cp( double T ) { return s_eval_Cv(T); }
    double s_eval_Q( double T, double A ) { return 0.0; }
};

class Two_level_electronic : public Species_energy_mode {
public:
    Two_level_electronic( int isp, double R, double min_massf, int g0, 
    	                  double theta0, int g1, double theta1 );
    ~Two_level_electronic() {}
    
private:
    int g0_;
    int g1_;
    double theta_;
    
    double s_eval_energy( double T );
    double s_eval_enthalpy( double T ) { return s_eval_energy(T); }
    double s_eval_entropy( double T, double p );
    double s_eval_Cv( double T );
    double s_eval_Cp( double T ) { return s_eval_Cv(T); }
    double s_eval_Q( double T, double A ) { return 0.0; }
};

class Multi_level_electronic : public Species_energy_mode {
public:
    Multi_level_electronic( int isp, double R, double min_massf, 
    			    std::vector<int> &g, std::vector<double> &theta);
    ~Multi_level_electronic() {}
    
private:
    std::vector<int> g_;
    std::vector<double> theta_;
    
    double s_eval_energy( double T );
    double s_eval_enthalpy( double T ) { return s_eval_energy(T); }
    double s_eval_entropy( double T, double p );
    double s_eval_Cv( double T );
    double s_eval_Cp( double T ) { return s_eval_Cv(T); }
    double s_eval_Q( double T, double A ) { return 0.0; }
};

class Fully_excited_translation : public Species_energy_mode {
public:
    Fully_excited_translation( int isp, double R, double min_massf );
    ~Fully_excited_translation() {}
    
private:
    double Cv_;
    double Cp_;
    double entropy_constant_;
    
    double s_eval_energy( double T );
    double s_eval_enthalpy( double T );
    double s_eval_entropy( double T, double p=0.0 );
    double s_eval_Cv( double T );
    double s_eval_Cp( double T );
    double s_eval_Q( double T, double A ) { return 0.0; }
};

class Fully_excited_rotation : public Species_energy_mode {
public:
    // NOTE: sigma = 1 or 2 for hetero/homonuclear molecules
    Fully_excited_rotation( int isp, double R, double min_massf, double theta, int sigma );
    ~Fully_excited_rotation() {}
    
    double get_theta()
    { return theta_; }
    
private:
    double Cv_;
    double theta_;
    int sigma_;
    
    double s_eval_energy( double T );
    double s_eval_enthalpy( double T ) { return s_eval_energy(T); }
    double s_eval_entropy( double T, double p=0.0 );
    double s_eval_Cv( double T );
    double s_eval_Cp( double T ) { return s_eval_Cv(T); }
    double s_eval_Q( double T, double A ) { return 0.0; }
};

class Fully_excited_nonlinear_rotation : public Species_energy_mode {
public:
    Fully_excited_nonlinear_rotation( int isp, double R, double min_massf, 
    				      double theta_A0, double theta_B0, double theta_C0,
    				      int sigma );
    ~Fully_excited_nonlinear_rotation() {}
    
private:
    double Cv_;
    double theta_A0_;
    double theta_B0_;
    double theta_C0_;
    int sigma_;
    
    double s_eval_energy( double T );
    double s_eval_enthalpy( double T ) { return s_eval_energy(T); }
    double s_eval_entropy( double T, double p=0.0 );
    double s_eval_Cv( double T );
    double s_eval_Cp( double T ) { return s_eval_Cv(T); }
    double s_eval_Q( double T, double A ) { return 0.0; }
};

class Vibration : public Species_energy_mode {
public:
    Vibration( int isp, double R, double min_massf, double theta );
    ~Vibration() {}
    
    double get_theta()
    { return theta_; }
    
protected:
    double theta_;
    
    double s_eval_energy( double T );
    double s_eval_enthalpy( double T ) { return s_eval_energy(T); }
    double s_eval_entropy( double T, double p=0.0 );
    double s_eval_Cv( double T );
    double s_eval_Cp( double T ) { return s_eval_Cv(T); }
    double s_eval_Q( double T, double A ) { return 0.0; }
};

class Harmonic_vibration : public Vibration {
public:
    Harmonic_vibration( int isp, double R, double min_massf, double theta );
    ~Harmonic_vibration() {}
    
private:
    double s_eval_energy( double T );
    double s_eval_entropy( double T, double p=0.0 );
    double s_eval_Cv( double T );
    double s_eval_Q( double T, double A );
};

class Truncated_harmonic_vibration : public Vibration {
public:
    Truncated_harmonic_vibration( int isp, double R, double min_massf, double theta, double theta_D_ );
    ~Truncated_harmonic_vibration() {}
    
    double eval_HO_energy_from_T( double T )
    { return s_eval_HO_energy( T ); }
    
private:
    double theta_D_;
    
    double s_eval_energy( double T );
    double s_eval_entropy( double T, double p=0.0 );
    double s_eval_Cv( double T );
    double s_eval_HO_energy( double T );
    double s_eval_Q( double T, double A );
};

class Fully_coupled_diatom_internal : public Species_energy_mode {
public:
    Fully_coupled_diatom_internal( int isp, double R, double min_massf, int sigma_r, std::vector< std::vector<double> > &elev_data );
    ~Fully_coupled_diatom_internal();
    
private:
    int sigma_r_;
    double m_;
    std::vector<Diatom_electronic_level*> elevs_;
    
    double s_eval_energy( double T );
    double s_eval_enthalpy( double T ) { return s_eval_energy(T); }
    double s_eval_entropy( double T, double p=0.0 );
    double s_eval_Cv( double T );
    double s_eval_Cp( double T ) { return s_eval_Cv(T); }
    double s_eval_Q( double T, double A ) { return 0.0; }
};

#endif
