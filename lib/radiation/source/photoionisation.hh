/** \file photoionisation.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 2-March-10: Initial implementation
 *
 *  \brief Declarations for photo-ionisation cross-section model classes
 *
 **/

#ifndef PHOTOIONISATION_HH
#define PHOTOIONISATION_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

class PhotoIonisationCrossSectionModel {
public:
    /// \brief Constructor
    PhotoIonisationCrossSectionModel( std::string name, double E_min=0.0 );
    
    /// \brief Deconstructor
    virtual ~PhotoIonisationCrossSectionModel() = 0;
    
public:
    std::string get_name()
    { return name; }

    double get_threshold_energy()
    { return E_min; }

    virtual void spectral_distribution( std::vector<double> &nus ) = 0;

    virtual double eval( double nu ) = 0;

protected:
    std::string name;
    double E_min;
};

class NoPICSModel : public PhotoIonisationCrossSectionModel {
public:
    /// \brief Constructor
    NoPICSModel();
    
    /// \brief Deconstructor
    ~NoPICSModel();
    
public:
    void spectral_distribution( std::vector<double> &nus );

    double eval( double nu );
    
};

class HydrogenicModel : public PhotoIonisationCrossSectionModel {
public:
    /// \brief Constructor
    HydrogenicModel( lua_State * L, double n_eff, int Z, double I, double E_min );
    
    /// \brief Deconstructor
    ~HydrogenicModel();
    
public:
    void spectral_distribution( std::vector<double> &nus );

    double eval( double nu );
    
private:
    double n_eff;
    double Z;
    double I;
    double constB;
    double constC;
    double constD;
    double constE;
    double E_max;	// Maximum energy (used for spectral distribution)
    int nnus;		// Number of points to use for adaptive distribution
};

class PICS_step {
public:
    PICS_step( double nu_a, double nu_b, double sigma_bf );
    ~PICS_step();
    
public:
    double nu_a;
    double nu_b;
    double sigma_bf;
};

class JohnstonStepModel : public PhotoIonisationCrossSectionModel {
public:
    /// \brief Constructor
    JohnstonStepModel( lua_State * L, int ilev );
    
    /// \brief Deconstructor
    ~JohnstonStepModel();
    
public:
    void spectral_distribution( std::vector<double> &nus );

    double eval( double nu );
    
private:
    std::vector<PICS_step*> steps;
};

class JohnstonThresholdModel : public PhotoIonisationCrossSectionModel {
public:
    /// \brief Constructor
    JohnstonThresholdModel( lua_State * L, int ilev, double E_min );
    
    /// \brief Deconstructor
    ~JohnstonThresholdModel();
    
public:
    void spectral_distribution( std::vector<double> &nus );

    double eval( double nu );
    
private:
    double nu_t;
    double sigma_bf_t;
    double theta;
    double E_max;	// Maximum energy (used for spectral distribution)
    int nnus;	// Number of points to use for adaptive distribution
};

class TOPBaseModel : public PhotoIonisationCrossSectionModel {
public:
    /// \brief Constructor
    TOPBaseModel(  lua_State * L, int ilev );
    
    /// \brief Deconstructor
    ~TOPBaseModel();
    
public:
    void spectral_distribution( std::vector<double> &nus );

    double eval( double nu );
    
private:
    std::vector<double> nu_list;
    std::vector<double> sigma_list;
    int i_prev;
    int N_points;
};

PhotoIonisationCrossSectionModel*
create_Hydrogenic_PICS_model( lua_State * L, int Z, double I, double E );

PhotoIonisationCrossSectionModel*
create_Johnston_PICS_model( lua_State * L, int ilev, double I, double E );

PhotoIonisationCrossSectionModel*
create_TOPBase_PICS_model( lua_State * L, int ilev, int Z, double I, double E );

PhotoIonisationCrossSectionModel*
create_new_PICS_model( lua_State * L, int ilev, int Z, double I, double E );

#endif
