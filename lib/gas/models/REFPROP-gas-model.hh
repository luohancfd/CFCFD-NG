/** \file REFPROP-gas-model.hh
 *  \ingroup gas
 *
 *  \author Peter J. Blyton
 *  \version 15-May-2012
 */

/**
 *  \brief Class for using the REFPROP thermodynamic database as a gas model.
 *
 *  Pure fluid (.FLD), pseudo pure fluid (.PPF) and pre-defined mixture (.MIX)
 *  files from REFPROP may be used in this gas model. The molecular weights
 *  and mole fractions are specified by these REFPROP fluid files, and are set
 *  up and saved locally to the gas model object when creating a gas model.
 *  This could be improved in the future by updating the mole fraction vector
 *  with the gas data object within each method call.
 *
 *  REFPROP uses units of mol, L, kPa and uPa.s. All unit conversion is performed
 *  in-function so function interfaces and the gas model and data objects remain
 *  in the correct units.
 *
 *  Refer to MANUAL.TXT provided with REFPROP for documentation on the functions called.
 */

#ifndef REFPROP_GAS_MODEL_HH
#define REFPROP_GAS_MODEL_HH

#include "gas_data.hh"
#include "gas-model.hh"

using namespace std;

// Function pointers to the Fortran functions and subroutines in REFPROP.
// Fortran subroutines are all call-by-reference, functions return a value.
typedef void (*fp_SETUP_TYPE)(long &,char*,char*,char*,long &,char*);
typedef void (*fp_SETMIX_TYPE)(char*,char*,char*,long &,char*,double *,long &,char*);
typedef double (*fp_WMOL_TYPE)(double *);// Fortran function, returns a value
typedef void (*fp_DEFLSH_TYPE)(double &,double &,double *,double &,double &,double &,
                               double &,double *,double *,double &,double &,double &,
                               double &,double &,double &,long &,char*);
typedef void (*fp_TPFLSH_TYPE)(double &,double &,double *,double &,double &,double &,
                               double *,double *,double &,double &,double &,double &,
                               double &,double &,double &,long &,char*);
typedef void (*fp_TDFLSH_TYPE)(double &,double &,double *,double &,double &,double &,
                               double *,double *,double &,double &,double &,double &,
                               double &,double &,double &,long &,char*);
typedef void (*fp_TRNPRP_TYPE)(double &,double &,double *,double &,double &,long &,char*);
typedef void (*fp_DEFL1_TYPE)(double &,double &,double *,double &,long &,char*);
typedef void (*fp_CVCP_TYPE)(double &,double &,double *,double &,double &);
typedef void (*fp_PRESS_TYPE)(double &,double &,double *,double &);
typedef void (*fp_SATT_TYPE)(double &,double *,long &,double &,double &,double &,
                             double *,double *,long &,char*);

class REFPROP_gas_model : public Gas_model {
public:
    REFPROP_gas_model(std::string cfile);
    // This constructor overrides the constructor of the base Gas_model class,
    // in constrast to other gas models which append to it. We override it as
    // the base Gas_model constructor expects the lua gas config file to provide
    // a molecular weight, however this is handled by REFPROP in this gas model.
    // This constructor has a few extra jobs that were performed by the base
    // Gas_model, such as filling the s_names_ vector.
    ~REFPROP_gas_model();

private:
    string phase;
    // Constants for REFPROP
    static const long reflength = 3; // 3 character thermo reference state
    static const long stringlength = 255; // file path and error message length
    static const long ncmax = 20; // maximum number of components in mixture

    // Components for REFPROP
    void* refprop_handle;
    fp_SETUP_TYPE SETUP;
    fp_SETMIX_TYPE SETMIX;
    fp_WMOL_TYPE WMOL;
    fp_DEFLSH_TYPE DEFLSH;
    fp_TPFLSH_TYPE TPFLSH;
    fp_TDFLSH_TYPE TDFLSH;
    fp_TRNPRP_TYPE TRNPRP;
    fp_DEFL1_TYPE DEFL1;
    fp_CVCP_TYPE CVCP;
    fp_PRESS_TYPE PRESS, DPDD, ENERGY, ENTHAL, ENTRO;
    fp_SATT_TYPE SATT;

    long i,ierr;
    double x[ncmax],xliq[ncmax],xvap[ncmax];
    char hfld[stringlength*ncmax],hrf[reflength],
         herr[stringlength],hfmix[stringlength],
         hfiles[stringlength*ncmax],setpath[stringlength];
    double wm,d,e,t,p,dl,dv,q,h,s,cv,cp,w,eta,tcx,dpdrho;
    double liq_SS, vap_SS;
    long phase_flag;

    int s_eval_thermo_state_rhoe(Gas_data &Q);
    int s_eval_thermo_state_pT(Gas_data &Q);
    int s_eval_thermo_state_rhoT(Gas_data &Q);
    int s_eval_transport_coefficients(Gas_data &Q);
    int s_eval_diffusion_coefficients(Gas_data &Q);
    double s_dedT_const_v(const Gas_data &Q, int &status);
    double s_dhdT_const_p(const Gas_data &Q, int &status);
    double s_internal_energy(const Gas_data &Q, int isp);
    double s_enthalpy(const Gas_data &Q, int isp);
    double s_entropy(const Gas_data &Q, int isp);

    void is_slib_error(const string &method_name, bool close_handle);
    void is_REFPROP_error(const string &method_name);
    double two_phase_sound_speed();
};

Gas_model* create_REFPROP_gas_model(string cfile);

#endif
