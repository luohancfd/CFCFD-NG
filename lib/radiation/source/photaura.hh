/** \file photaura.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 28-Jul-09 : initial framework
 *
 **/

#ifndef PHOTAURA_HH
#define PHOTAURA_HH

#include <vector>

#include "../../util/source/useful.h"
#include "spectral_model.hh"
#include "radiator.hh"
#include "diatomic_radiator.hh"
#include "atomic_radiator.hh"

/* 'Hardcoded' Photaura control parameters */
#define MIN_CONC 1.0e-20
#define DELTA_NU_MIN 3.0e9	// min frequency interval in Hz for optimization

class Photaura : public RadiationSpectralModel {
public:
    Photaura( lua_State * L );
    Photaura( std::string input_file );

    ~Photaura();
    
    std::string str() const;
    
    int get_nrad()
    { return nrad; }
    
    Radiator * get_radiator_pointer( int irad )
    { return radiators[irad]; }
    
    Radiator * get_radiator_pointer_from_name( std::string name );
    
    DiatomicRadiator * get_diatomic_radiator_pointer_from_name( std::string name );
    
    AtomicRadiator * get_atomic_radiator_pointer_from_name( std::string name );
    
    std::string get_rad_name( int irad )
    { return radiators[irad]->name; }
    
    double integrate_emission_spectrum_for_radiator( Gas_data &Q, int irad, bool write_to_file=false );

private:
    
    void set_FirstOrderLTNE_pointers( int irad, int e_index );
    
    double integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved );
    
    double variably_integrated_emission_for_gas_state( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u, bool spectrally_resolved );
    
    double integrate_emission_spectrum( Gas_data &Q );
    
    double sum_emission_coefficients( Gas_data &Q );
    
    double sum_optically_variable_emission_coefficients(Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u );
    
    void spectra_for_gas_state( Gas_data &Q, CoeffSpectra &X );
    
    void spectral_distribution_for_gas_state(Gas_data &Q, std::vector<double> &nus);
    
    void write_line_widths( Gas_data &Q );
    
    void initialise_all_radiators( Gas_data &Q );
    
    double get_rad_conc( Gas_data &Q, int irad );
    
    void prep_rad_pop_files();
    
    void append_current_rad_pops( double x );
    
    void write_QSS_analysis_files( Gas_data &Q, int index );
    
private:
    
    int nrad;
    
    std::vector<Radiator*> radiators;
};

void impose_min_interval( std::vector<double> &X, double delta_X_min );

#endif
