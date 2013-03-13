/** \file spradian.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 11-Jan-12 : initial framework
 *
 **/

#ifndef SPRADIAN_HH
#define SPRADIAN_HH

#define SPR_N_ATM_NEQS  6
#define SPR_N_ATM_RADS 18
#define SPR_N_DTM_NEQS 10
#define SPR_N_DTM_BNDS 40
#define SPR_N_TRM_BNDS 10

#include <vector>

#include "../../util/source/useful.h"
#include "spectral_model.hh"
#include "spradian_radiator.hh"

extern "C" int radipac_( double* z, int* nnode, double* depth, int* method, double* tblack, double* stand_off, double* nose_radius, 
    			 double* wavmin, double* wavmax, int* nwav, int* avg_num,
			 char (*atom_noneqs)[2], char (*atom_rads)[2][2], char (*diatom_noneqs)[4], char (*diatom_bands)[2][4], char (*triatom_bands)[2][4],
			 double* tran, double* trot, double* tvib, double* telec,
			 double* concC, double* concC2, double* concC2H, double* concC3, double* concCH, double* concCN, double* concCO, double* concCp, double* concH, double* concH2,
			 double* concHp, double* concN, double* concN2, double* concN2p, double* concNE, double* concNO, double* concNp, double* concO, double* concO2, double* concOH,
			 double* concOp, double* conchvy, double* concatm, double* concmol, double* avg_molwt, double (*norm_int)[1], double (*dibydx)[1], double (*flux)[1], double (*dqbydx)[1],
			 double* divq );

class SpradianParams {
public:
    SpradianParams( std::vector<SpradianRadiator*> &rad_vec, double lambda_min, double lambda_max, int spectral_points );

    SpradianParams( std::string fname );

    ~SpradianParams();

public:
    void set_char_arrays_to_white_space();

    void reset_variable_parameters();

    void write_to_file( std::string fname );

    void read_from_file( std::string fname );

    int call_radipac();

public:
    int method;

    double depth;
    double stand_off;
    double nose_radius;

    int nnode;
    double z[1];

    double wavmin;
    double wavmax;
    int nwav;
    int avg_num;

    double conchvy[1];
    double concC[1];
    double concC2[1];
    double concC2H[1];
    double concC3[1];
    double concCH[1];
    double concCN[1];
    double concCO[1];
    double concCp[1];
    double concH[1];
    double concH2[1];
    double concHp[1];
    double concN[1];
    double concN2[1];
    double concN2p[1];
    double concNE[1];
    double concNO[1];
    double concNp[1];
    double concO[1];
    double concO2[1];
    double concOH[1];
    double concOp[1];
    double concatm[1];
    double concmol[1];
    // NOTE: these are not used by the spradian subroutine
    double concC2p[1];
    double concCOp[1];
    double concO2p[1];

    double tran[1];
    double trot[1];
    double tvib[1];
    double telec[1];

    double avg_molwt[1];

    double tblack[2];

    char atom_noneqs[SPR_N_ATM_NEQS][2];
    char atom_rads[SPR_N_ATM_RADS][2][2];
    char diatom_noneqs[SPR_N_DTM_NEQS][4];
    char diatom_bands[SPR_N_DTM_BNDS][2][4];
    char triatom_bands[SPR_N_TRM_BNDS][2][4];

    double norm_int[2][1], dibydx[2][1], flux[2][1], dqbydx[2][1], divq[1];
};

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

    void setup_parameters( Gas_data &Q );

    SpradianParams * get_params_pointer()
    { return params; }

private:
    double integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved );
    
    double variably_integrated_emission_for_gas_state( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u, bool spectrally_resolved );
    
    void spectra_for_gas_state( Gas_data &Q, CoeffSpectra &X );
    
    void write_line_widths( Gas_data &Q );
    
    void prep_rad_pop_files();
    
    void append_current_rad_pops( double x );
    
    void write_QSS_analysis_files( Gas_data &Q, int index );

private:
    int nrad;

    int e_index;

    int iT, iTr, iTv, iTe;

    std::vector<std::string> rad_names;

    std::vector<SpradianRadiator*> radiators;

    SpradianParams * params;
};

#endif
