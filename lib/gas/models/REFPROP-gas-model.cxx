/** \file REFPROP-gas-model.cxx
 *  \ingroup gas
 *
 *  \author Peter J. Blyton
 *  \version 15-May-2012
 */

#include <iostream>
#include <string>
#include <dlfcn.h>
#include <string.h>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "REFPROP-gas-model.hh"

using namespace std;

REFPROP_gas_model::
REFPROP_gas_model(string cfile)
{
    lua_State *L = initialise_lua_State();

    if( luaL_dofile(L, cfile.c_str()) != 0 ) {
        ostringstream ost;
        ost << "REFPROP_gas_model():\n"
            << "Error in gas model input file: " << cfile << endl;
        input_error(ost);
    }
    lua_getglobal(L, "species");
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "REFPROP_gas_model::REFPROP_gas_model():\n"
            << "Error in the declaration of species: a table is expected.\n";
        input_error(ost);
    }
    if ( lua_objlen(L, -1) != 2 ) {
        ostringstream ost;
        ost << "REFPROP_gas_model::REFPROP_gas_model():\n"
            << "Error in the declaration of species: a single filename expected,\n"
            << "with the second item specifying 'single phase' or 'two phase'.\n";
        input_error(ost);
    }

    lua_rawgeti(L, -1, 1); // Put the specific species name to TOS
    const char* sp = luaL_checkstring(L, -1);
    lua_pop(L, 1); // pop specific species name off stack
    lua_rawgeti(L, -1, 2); // Put the phase type to TOS
    phase = luaL_checkstring(L, -1);
    lua_pop(L, 1); // pop phase type off stack
    lua_pop(L, 1); // pop list of species table off stack

    // Setup the REFPROP fluid model
    refprop_handle = dlopen("refprop.so", RTLD_LAZY);
    if (!refprop_handle) is_slib_error("REFPROP_gas_model", false);

    // Get pointers into the shared object for the required functions
    *(void **)(&SETUP) = dlsym(refprop_handle, "setup");
    if (!SETUP) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&SETMIX) = dlsym(refprop_handle, "setmix");
    if (!SETMIX) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&WMOL) = dlsym(refprop_handle, "wmol");
    if (!WMOL) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&DEFLSH) = dlsym(refprop_handle, "deflsh");
    if (!DEFLSH) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&TPFLSH) = dlsym(refprop_handle, "tpflsh");
    if (!TPFLSH) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&TDFLSH) = dlsym(refprop_handle, "tdflsh");
    if (!TDFLSH) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&TRNPRP) = dlsym(refprop_handle, "trnprp");
    if (!TRNPRP) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&DEFL1) = dlsym(refprop_handle, "defl1");
    if (!DEFL1) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&CVCP) = dlsym(refprop_handle, "cvcp");
    if (!CVCP) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&PRESS) = dlsym(refprop_handle, "press");
    if (!PRESS) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&DPDD) = dlsym(refprop_handle, "dpdd");
    if (!DPDD) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&ENERGY) = dlsym(refprop_handle, "energy");
    if (!ENERGY) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&ENTHAL) = dlsym(refprop_handle, "enthal");
    if (!ENTHAL) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&ENTRO) = dlsym(refprop_handle, "entro");
    if (!ENTRO) is_slib_error("REFPROP_gas_model", true);
    *(void **)(&SATT) = dlsym(refprop_handle, "satt");
    if (!SATT) is_slib_error("REFPROP_gas_model", true);

    string home = getenv("HOME"), species_path;
    string mixfile_path = home + "/e3bin/species/refprop/fluids/HMX.BNC";
    strncpy(hfmix, mixfile_path.c_str(), stringlength);
    strncpy(hrf, "DEF", reflength); // Default thermo reference state
    strncpy(herr, "Ok", stringlength);
    ierr = 0;

    if (strstr(sp, ".FLD") || strstr(sp, ".PPF")){
        species_path = home + "/e3bin/species/refprop/fluids/" + sp;
        strncpy(hfld, species_path.c_str(), stringlength*ncmax);
        i = 1; x[0] = 1.0;
        SETUP(i,hfld,hfmix,hrf,ierr,herr);
        if (ierr != 0) is_REFPROP_error("REFPROP_gas_model");
    } else if (strstr(sp, ".MIX")) {
        species_path = home + "/e3bin/species/refprop/mixtures/" + sp;
        strncpy(hfld, species_path.c_str(), stringlength*ncmax);
        SETMIX(hfld,hfmix,hrf,i,hfiles,x,ierr,herr);
        if (ierr != 0) is_REFPROP_error("REFPROP_gas_model");
    } else {
        ostringstream ost;
        ost << "REFPROP_gas_model::REFPROP_gas_model():\n"
            << "You can only provide a .FLD, .PPF or .MIX file." << endl;
        dlclose(refprop_handle);
        input_error(ost);
    }
    // Even for REFPROP mixtures, let lib/gas think it is single species
    set_number_of_species(1);
    set_number_of_modes(1);
    wm = WMOL(x); // Mixture molecular weight, g/mol
    M_.resize(1);
    M_[0] = wm*1e-3; // kg/mol
    s_names_.resize(1);
    s_names_[0] = string(sp);
    lua_close(L);
}

REFPROP_gas_model::
~REFPROP_gas_model()
{
    dlclose(refprop_handle);
}

int
REFPROP_gas_model::
s_eval_thermo_state_rhoe(Gas_data &Q)
{   // Calculate the thermodynamic state from density, energy and composition
    d = Q.rho/wm; // kg/m^3 to mol/L
    e = Q.e[0]*wm*1e-3; // J/kg to J/mol
    if (phase == "two phase") {
        // Use the general purpose flash calculator when the phase is not known.
        // This is quite slow, but allows two-phase calculations.
        DEFLSH(d,e,x,t,p,dl,dv,xliq,xvap,q,h,s,cv,cp,w,ierr,herr);
        if (ierr != 0) is_REFPROP_error("s_eval_thermo_state_rhoe");
        Q.p = p*1e3; // kPa to Pa
        Q.T[0] = t;
        Q.quality = q;
        Q.a = two_phase_sound_speed();
        return SUCCESS;
    } else if (phase == "single phase") {
        // Use the single phase flash calculator, much faster than DEFLSH.
        DEFL1(d,e,x,t,ierr,herr);
        if (ierr != 0) is_REFPROP_error("s_eval_thermo_state_rhoe");
        Q.T[0] = t;
        // Now update the pressure from the new temperature and density
        PRESS(t,d,x,p);
        Q.p = p*1e3; // kPa to Pa
        CVCP(t,d,x,cv,cp); // J/(mol.K)
        DPDD(t,d,x,dpdrho); // kPa.L/mol
        Q.a = sqrt(1e3*dpdrho/wm*(cp/cv));
        Q.quality = 1;
        return SUCCESS;
    } else {
        cout << "REFPROP_gas_model::s_eval_thermo_state_rhoe():\n"
             << "Second species item should specify either\n"
             << "'single phase' or 'two phase' for calculation.\n";
        exit(BAD_INPUT_ERROR);
    }
}

int
REFPROP_gas_model::
s_eval_thermo_state_pT(Gas_data &Q)
{   // Calculate the thermodynamic state from pressure and temperature
    t = Q.T[0];
    p = Q.p*1e-3; // Pa to kPa
    TPFLSH(t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr);
    if (ierr != 0) is_REFPROP_error("s_eval_thermo_state_pT");
    Q.rho = d*wm; // mol/L to mol/L
    Q.a = w;
    Q.e[0] = e*1e3/wm; // J/mol to J/kg
    Q.quality = q;
    return SUCCESS;
}

int
REFPROP_gas_model::
s_eval_thermo_state_rhoT(Gas_data &Q)
{   // Calculate pressure from density and temperature
    t = Q.T[0];
    d = Q.rho/wm; // kg/m^3 to mol/L
    if (phase == "two phase") {
        // Use the general purpose flash calculator when the phase is not known.
        // This is quite slow, but allows two-phase calculations.
        TDFLSH(t,d,x,p,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr);
        if (ierr != 0) is_REFPROP_error("s_eval_thermo_state_rhoT");
        Q.p = p*1e3; // kPa to Pa
        Q.e[0] = e*1e3/wm;
        Q.quality = q;
        Q.a = two_phase_sound_speed();
        return SUCCESS;
    } else if (phase == "single phase") {
        PRESS(t,d,x,p);
        Q.p = p*1e3; // kPa to Pa
        ENERGY(t,d,x,e);
        Q.e[0] = e*1e3/wm; // J/mol to J/kg
        CVCP(t,d,x,cv,cp); // J/(mol.K)
        DPDD(t,d,x,dpdrho); // kPa.L/mol
        Q.a = sqrt(1e3*dpdrho/wm*(cp/cv));
        Q.quality = 1;
        return SUCCESS;
    } else {
        cout << "REFPROP_gas_model::s_eval_thermo_state_rhoT():\n"
             << "Second species item should specify either\n"
             << "'single phase' or 'two phase' for calculation.\n";
        exit(BAD_INPUT_ERROR);
    }
}

int
REFPROP_gas_model::
s_eval_transport_coefficients(Gas_data &Q)
{   // Transport properties from temperature, density and composition.
    t = Q.T[0];
    d = Q.rho/wm; // kg/m^3 to mol/L
    TRNPRP(t,d,x,eta,tcx,ierr,herr);
    if (ierr != 0) is_REFPROP_error("s_eval_transport_coefficients");
    Q.mu = eta*1e6; // uPa.s to Pa.s
    Q.k[0] = tcx;
    return SUCCESS;
}

int
REFPROP_gas_model::
s_eval_diffusion_coefficients(Gas_data &Q)
{
    cout << "REFPROP_gas_model::s_eval_diffusion_coefficients()\n"
         << "Diffusion coefficients have no meaning for this gas model.\n";
    exit(BAD_INPUT_ERROR);
}

double
REFPROP_gas_model::
s_dedT_const_v(const Gas_data &Q, int &status)
{   // Cv from temperature, density and composition.
    t = Q.T[0];
    d = Q.rho/wm; // kg/m^3 to mol/L
    CVCP(t,d,x,cv,cp);
    return cv*1e3/wm; // J/(mol.K) to J/(kg.K)
}

double
REFPROP_gas_model::
s_dhdT_const_p(const Gas_data &Q, int &status)
{   // Cp from temperature, density and composition.
    t = Q.T[0];
    d = Q.rho/wm; // kg/m^3 to mol/L
    CVCP(t,d,x,cv,cp);
    return cp*1e3/wm; // J/(mol.K) to J/(kg.K)
}

double
REFPROP_gas_model::
s_internal_energy(const Gas_data &Q, int isp)
{   // Internal energy from temperature, density and composition.
    t = Q.T[0];
    d = Q.rho/wm; // kg/m^3 to mol/L
    ENERGY(t,d,x,e);
    return e*1e3/wm; // J/mol to J/kg
}

double
REFPROP_gas_model::
s_enthalpy(const Gas_data &Q, int isp)
{   // Enthalpy from temperature, density and composition.
    t = Q.T[0];
    d = Q.rho/wm; // kg/m^3 to mol/L
    ENTHAL(t,d,x,h);
    return h*1e3/wm; // J/mol to J/kg
}

double
REFPROP_gas_model::
s_entropy(const Gas_data &Q, int isp)
{   // Entropy from temperature, density and composition.
    t = Q.T[0];
    d = Q.rho/wm; // kg/m^3 to mol/L
    ENTRO(t,d,x,s);
    return s*1e3/wm; // J/(mol.K) to J/(kg.K)
}

void
REFPROP_gas_model::
is_slib_error(const string &method_name, bool close_handle)
{   // Stop the program as we have an issue setting up the shared library
    ostringstream ost;
    ost << "REFPROP_gas_model::" + method_name << + "():\n"
        << dlerror() << endl;
    if (close_handle) dlclose(refprop_handle);
    input_error(ost);
}

void
REFPROP_gas_model::
is_REFPROP_error(const string &method_name)
{   // Stop the program as our REFPROP function failed
    ostringstream ost;
    ost << "REFPROP_gas_model::" + method_name << + "():\n"
        << herr << endl;
    dlclose(refprop_handle);
    input_error(ost);
}

double
REFPROP_gas_model::
two_phase_sound_speed()
{   // Calculate mass average of liquid and vapor sound speeds if two phase
    if (q <= 0.0 || q >= 1.0) { // Single phase
        return w; // Already computed sound speed
    } else { // Two phase
        phase_flag = 1;
        SATT(t,x,phase_flag,p,dl,dv,xliq,xvap,ierr,herr);
        if (ierr != 0) is_REFPROP_error("two_phase_sound_speed");
        // Liquid sound speed
        CVCP(t,dl,x,cv,cp); // J/(mol.K)
        DPDD(t,dl,x,dpdrho); // kPa.L/mol
        liq_SS = sqrt(1e3*dpdrho/wm*(cp/cv));
        // Vapour sound speed
        CVCP(t,dv,x,cv,cp); // J/(mol.K)
        DPDD(t,dv,x,dpdrho); // kPa.L/mol
        vap_SS = sqrt(1e3*dpdrho/wm*(cp/cv));
        return q*vap_SS + (1.0 - q)*liq_SS;
    }
}

Gas_model* create_REFPROP_gas_model(string cfile)
{
    return new REFPROP_gas_model(cfile);
}
