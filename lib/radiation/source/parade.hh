/** \file parade.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 11-Jan-12 : initial framework
 *
 **/

#ifndef PARADE_HH
#define PARADE_HH

#include <vector>

#include "../../util/source/useful.h"
#include "spectral_model.hh"
#include "parade_radiator.hh"

#define USE_FLO_INPUT_FILES 1

class Parade : public RadiationSpectralModel {
public:
    Parade( lua_State * L );

    Parade( const std::string input_file );

    void initialise( lua_State * L );

    ~Parade();
    
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

    void read_control_template_file( std::string control_template_filename );

    void read_snbopt_template_file( std::string snbcon_template_filename );

    void create_parade_control_files( Gas_data &Q );

private:
    int nrad;

    int e_index;

    int iT, iTe;

    std::vector<std::string> rad_names;

    std::vector<ParadeRadiator*> radiators;

    std::stringstream control_template_file_buffer;

    std::stringstream snbopt_template_file_buffer;
};

#endif
