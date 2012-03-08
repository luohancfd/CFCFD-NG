/** \file spradian.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 11-Jan-12 : initial framework
 *
 **/

#ifndef SPRADIAN_HH
#define SPRADIAN_HH

#include <vector>

#include "../../util/source/useful.h"
#include "spectral_model.hh"
#include "spradian_radiator.hh"

class Spradian : public RadiationSpectralModel {
public:
    Spradian( lua_State * L );

    Spradian( const std::string input_file );

    void initialise( lua_State * L );

    ~Spradian();
    
    std::string str() const;
    
    int get_nrad()
    { return nrad; }

    std::string get_rad_name( int irad )
    { return rad_names[irad]; }

private:
    double integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved );
    
    double variably_integrated_emission_for_gas_state( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u, bool spectrally_resolved );
    
    void spectra_for_gas_state( Gas_data &Q, CoeffSpectra &X );
    
    void spectral_distribution_for_gas_state(Gas_data &Q, std::vector<double> &nus);
    
    void write_line_widths( Gas_data &Q );
    
    void prep_rad_pop_files();
    
    void append_current_rad_pops( double x );
    
    void write_QSS_analysis_files( Gas_data &Q, int index );

    void read_spradian_template_file( std::string spradian_template_filename );

    void create_spradian_control_file( Gas_data &Q );

private:
    int nrad;

    int e_index;

    int iT, iTe;

    std::vector<std::string> rad_names;

    std::vector<SpradianRadiator*> radiators;

    std::stringstream spradian_template_file_buffer;
};

#endif
