/** \file spectral_model.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 10-Oct-09
 *
 **/

#ifndef SPECTRAL_MODEL_HH
#define SPECTRAL_MODEL_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../gas/models/gas_data.hh"
#include "spectra_pieces.hh"

#define ECHO_RAD_INPUT 1
#define DEBUG_RAD      0

class Radiator;

class RadiationSpectralModel {
public:
    RadiationSpectralModel();
    
    RadiationSpectralModel( lua_State *L );
    
    RadiationSpectralModel( const std::string input_file );
    
    virtual ~RadiationSpectralModel();
    
    virtual std::string str() const = 0;
    
    double radiative_integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved )
    { return integrated_emission_for_gas_state( Q, spectrally_resolved ); }
    
    double radiative_variably_integrated_emission_for_gas_state( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u, bool spectrally_resolved )
    { return variably_integrated_emission_for_gas_state( Q, wavel_switch, Lambda_l, Lambda_u, spectrally_resolved ); }
    
    void radiative_spectra_for_gas_state( Gas_data &Q, CoeffSpectra &X )
    { return spectra_for_gas_state( Q, X ); }
    
    void radiative_spectral_distribution_for_gas_state( Gas_data &Q, std::vector<double> &nus )
    { return spectral_distribution_for_gas_state( Q, nus ); }
    
    void reset_spectral_params( int isb );
    
    double get_lambda_min() { return lambda_min_star; }
    
    double get_lambda_max() { return lambda_max_star; }
    
    int get_spectral_points() { return spectral_points_star; }
    
    int get_spectral_blocks() { return spectral_blocks; }
    
    void new_spectral_params( double _lambda_min, double _lambda_max, int _spectral_points, int _spectral_blocks );

    double get_delta_nu() { return delta_nu; }

    bool get_adaptive_spectral_grid() { return adaptive_spectral_grid; }

    void write_line_widths_to_file( Gas_data &Q )
    { write_line_widths(Q); }
    
    void prep_radiator_population_files()
    { return prep_rad_pop_files(); }
    
    void append_current_radiator_populations( double x )
    { return append_current_rad_pops(x); }
    
    void write_QSS_population_analysis_files( Gas_data &Q, int index )
    { return write_QSS_analysis_files( Q, index ); }

protected:
    virtual double integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved ) = 0;
    
    virtual double variably_integrated_emission_for_gas_state( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u, bool spectrally_resolved ) = 0;
    
    virtual void spectra_for_gas_state( Gas_data &Q, CoeffSpectra &X ) = 0;
    
    virtual void spectral_distribution_for_gas_state( Gas_data &Q, std::vector<double> &nus ) = 0;
    
    virtual void write_line_widths( Gas_data &Q ) = 0;
    
    virtual void prep_rad_pop_files();
    
    virtual void append_current_rad_pops( double x );
    
    virtual void write_QSS_analysis_files( Gas_data &Q, int index );
    
protected:
    double lambda_min, lambda_min_star;
    double lambda_max, lambda_max_star;
    int spectral_points, spectral_points_star;
    int spectral_blocks, spectral_block;
    double delta_nu;
    bool adaptive_spectral_grid;
};

/* Functions for creating a RadiationSpectralModel */

RadiationSpectralModel * create_radiation_spectral_model( const std::string input_file="dummy_file" );

lua_State* initialise_radiation_lua_State();

#endif
