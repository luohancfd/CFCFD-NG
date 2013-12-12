// Author: Rowan J. Gollan
// Date: 07-Nov-2008
//       14-Jun-2010 (PJ) small variable change to fix bug
//                        and eliminate ambiguity
//       12-Dec-2013 (PJ) Add entropy to interpolation data.
// Place: Hampton, Virginia, USA
//        Lisbon, Portugal
// Note:
//   This is a port of PJs look-up table
//   implementation with some cosmetic changes:
//     - gzipped text files, instead of binary format
//     - vector storage instead of arrays
//     - reworked to fit in new class framework
//

#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <math.h>
#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "physical_constants.hh"
#include "look-up-table.hh"

using namespace std;

Look_up_table::Look_up_table(const string cfile)
    : Gas_model()
{
    set_number_of_species(1);
    set_number_of_modes(1);
    s_names_.resize(1);
    s_names_[0] = "LUT";

    // Read lua file and populate LUT
    lua_State *L = initialise_lua_State();
    
    if ( do_gzfile(L, cfile.c_str()) != 0 ) {
	ostringstream ost;
	ost << "Look_up_table():\n";
	ost << "    Error in look-up table input file: " << cfile << endl;
	ost << "    Something is wrong with the Lua script itself." << endl;
	input_error(ost);
    }

    int ie, ir;

    try {
	with_entropy = get_int(L, LUA_GLOBALSINDEX, "with_entropy") == 1;
    } catch ( runtime_error &e ) {
	with_entropy = false;
    }
    if ( !with_entropy ) cout << "Look_up_table(): No entropy data available." << endl;
    iesteps_ = get_positive_int(L, LUA_GLOBALSINDEX, "iesteps");
    irsteps_ = get_positive_int(L, LUA_GLOBALSINDEX, "irsteps");
    emin_ = get_number(L, LUA_GLOBALSINDEX, "emin");
    de_ = get_positive_number(L, LUA_GLOBALSINDEX, "de");
    lrmin_ = get_number(L, LUA_GLOBALSINDEX, "lrmin");
    dlr_ = get_positive_number(L, LUA_GLOBALSINDEX, "dlr");

    emax_ = emin_ + de_ * iesteps_;
    lrmax_ = lrmin_ + dlr_ * irsteps_;
    
    lua_getglobal(L, "data");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Look_up_table():\n";
	ost << "    Error in look-up table input file: " << cfile << endl;
	ost << "    A table of 'data' is expected, but not found.\n";
	input_error(ost);
    }

    int ne = lua_objlen(L, -1);
    if ( ne != iesteps_ + 1 ) {
	ostringstream ost;
	ost << "Look_up_table():\n"
	    << "    Error in look-up table input file: " << cfile << endl
	    << "    Inconsistent numbers for energy steps: "
	    << "points=" << ne << " steps=" << iesteps_ << endl;
	input_error(ost);
    }
    Cv_hat_.resize(ne);
    Cv_.resize(ne);
    R_hat_.resize(ne);
    g_hat_.resize(ne);
    mu_hat_.resize(ne);
    k_hat_.resize(ne);
    if ( with_entropy ) s_.resize(ne);

    lua_rawgeti(L, -1, 1);
    int nr = lua_objlen(L, -1);
    lua_pop(L, 1);
    if ( nr != irsteps_ + 1 ) {
	ostringstream ost;
	ost << "Look_up_table():\n"
	    << "    Error in look-up table input file: " << cfile << endl
	    << "    Inconsistent numbers for density steps: "
	    << "points=" << nr << " steps=" << irsteps_ << endl;
	input_error(ost);
    }

    for ( ie = 0; ie < ne; ++ie ) {
	Cv_hat_[ie].resize(nr);
	Cv_[ie].resize(nr);
	R_hat_[ie].resize(nr);
	g_hat_[ie].resize(nr);
	mu_hat_[ie].resize(nr);
	k_hat_[ie].resize(nr);
	if ( with_entropy ) s_[ie].resize(nr);

	lua_rawgeti(L, -1, ie+1);
	for ( ir = 0; ir < nr; ++ir ) {
	    lua_rawgeti(L, -1, ir+1);

	    lua_rawgeti(L, -1, 1);
	    Cv_hat_[ie][ir] = luaL_checknumber(L, -1);
	    lua_pop(L, 1);

	    lua_rawgeti(L, -1, 2);
	    Cv_[ie][ir] = luaL_checknumber(L, -1);
	    lua_pop(L, 1);

	    lua_rawgeti(L, -1, 3);
	    R_hat_[ie][ir] = luaL_checknumber(L, -1);
	    lua_pop(L, 1);

	    lua_rawgeti(L, -1, 4);
	    g_hat_[ie][ir] = luaL_checknumber(L, -1);
	    lua_pop(L, 1);

	    lua_rawgeti(L, -1, 5);
	    mu_hat_[ie][ir] = luaL_checknumber(L, -1);
	    lua_pop(L, 1);
	    
	    lua_rawgeti(L, -1, 6);
	    k_hat_[ie][ir] = luaL_checknumber(L, -1);
	    lua_pop(L, 1);
	    
	    if ( with_entropy ) {
		lua_rawgeti(L, -1, 7);
		s_[ie][ir] = luaL_checknumber(L, -1);
		lua_pop(L, 1);
	    }

	    lua_pop(L, 1); // pop data[ie][ir] off.
	}
	lua_pop(L, 1); // pop data[ie] off.
    }

    lua_pop(L, 1); // pop data off.

    lua_close(L);
} // end of constructor

Look_up_table::~Look_up_table() 
{
    s_names_.clear();
    size_t ne = Cv_hat_[0].size();
    for ( size_t ie = 0; ie < ne; ++ie ) {
	Cv_hat_[ie].clear();
	Cv_[ie].clear();
	R_hat_[ie].clear();
	g_hat_[ie].clear();
	mu_hat_[ie].clear();
	k_hat_[ie].clear();
	if ( with_entropy ) s_[ie].clear();
    }
    Cv_hat_.clear();
    Cv_.clear();
    R_hat_.clear();
    g_hat_.clear();
    mu_hat_.clear();
    k_hat_.clear();
    if ( with_entropy ) s_.clear();
} // end of destructor

int
Look_up_table::
determine_interpolants(const Gas_data &Q, int &ir, int &ie, double &lrfrac, double &efrac)
{
    if ( Q.rho <= 0.0 ) {
	cout << "Look_up_table::determine_interpolants(): density= " << Q.rho 
	     << " is zero or negative\n";
	cout << "   Supplied Q:" << endl;
	Q.print_values();
	return VALUE_ERROR;
    }

    // Find the enclosing cell. 
    double logrho = log10(Q.rho);
    ir = (int) ((logrho - lrmin_) / dlr_);
    ie = (int) ((Q.e[0] - emin_) / de_);
    // cout << "ir= " << ir << " ie= " << ie << endl;
    // cout << "dlr= " << dlr_ << " de= " << de_ << endl;
    // cout << "lrmin= " << lrmin_ << endl;
    // cout << "emin= " << emin_ << endl;
    
    // Make sure that we don't try to access data outside the
    // actual arrays.
    if ( ir < 0 ) ir = 0;
    if ( ir > (irsteps_ - 1) ) ir = irsteps_ - 1;
    if ( ie < 0 ) ie = 0;
    if ( ie > (iesteps_ - 1) ) ie = iesteps_ - 1;

    // Calculate bilinear interpolation(/extrapolation) fractions.
    lrfrac = (logrho - (lrmin_ + ir * dlr_)) / dlr_;
    efrac  = (Q.e[0] - (emin_ + ie * de_)) / de_;
    // cout << "desired lrfrac= " << lrfrac << " efrac= " << efrac << endl;
#   if 0
    // 2013-12-13 (PJ) removed the limitation on extrapolation because
    // it was messing with my estimates of entropy for low temperatures.
    // Limit the extrapolation to small distances.
    constexpr double EXTRAP_MARGIN = 1.0;
    lrfrac = max(lrfrac, -EXTRAP_MARGIN);
    lrfrac = min(lrfrac, 1.0+EXTRAP_MARGIN);
    efrac = max(efrac, -EXTRAP_MARGIN);
    efrac = min(efrac, 1.0+EXTRAP_MARGIN);
    // cout << "actual lrfrac= " << lrfrac << " efrac= " << efrac << endl;
#   endif

    return (SUCCESS);
}

int
Look_up_table::
s_eval_thermo_state_rhoe(Gas_data &Q)
{
    double efrac, lrfrac, Cv_eff, R_eff, g_eff;
    int    ir, ie;

    if ( determine_interpolants(Q, ir, ie, lrfrac, efrac) != SUCCESS ) {
	cout << "Bailing out!\n";
	exit(1);
    }

    Cv_eff = (1.0 - efrac) * (1.0 - lrfrac) * Cv_hat_[ie][ir] +
	efrac         * (1.0 - lrfrac) * Cv_hat_[ie+1][ir] +
	efrac         * lrfrac         * Cv_hat_[ie+1][ir+1] +
	(1.0 - efrac) * lrfrac         * Cv_hat_[ie][ir+1];
   
    R_eff  = (1.0 - efrac) * (1.0 - lrfrac) * R_hat_[ie][ir] +
	efrac         * (1.0 - lrfrac) * R_hat_[ie+1][ir] +
	efrac         * lrfrac         * R_hat_[ie+1][ir+1] +
	(1.0 - efrac) * lrfrac         * R_hat_[ie][ir+1];

    g_eff  = (1.0 - efrac) * (1.0 - lrfrac) * g_hat_[ie][ir] +
	efrac         * (1.0 - lrfrac) * g_hat_[ie+1][ir] +
	efrac         * lrfrac         * g_hat_[ie+1][ir+1] +
	(1.0 - efrac) * lrfrac         * g_hat_[ie][ir+1];

    // Reconstruct the thermodynamic properties.
    Q.T[0] = Q.e[0] / Cv_eff;
    Q.p = Q.rho*R_eff*Q.T[0];
    Q.a = sqrt(g_eff*R_eff*Q.T[0]);
    if ( Q.T[0] < gas_Tmin() ) {
	cout << "Look_up_table::eval_thermo_state_rhoe(): Low temperature: rho= " << Q.rho
	    << " e= " << Q.e[0]
	    << " T= " << Q.T[0] << endl;
	cout << "   Supplied Q:" << endl;
       Q.print_values();
       return BAD_TEMPERATURE_ERROR;
   }

   // Fix meaningless values if they arise
   if ( Q.p < 0.0 ) Q.p = 0.0;
   if ( Q.T[0] < 0.0 ) Q.T[0] = 0.0;
   if ( Q.a < 0.0 ) Q.a = 0.0;

   return (SUCCESS);
}

int
Look_up_table::
s_eval_transport_coefficients(Gas_data &Q)
{
    double efrac, lrfrac;
    double mu_eff, k_eff;
    int    ir, ie;

    if ( determine_interpolants(Q, ir, ie, lrfrac, efrac) != SUCCESS ) {
	cout << "Bailing out!\n";
	exit(1);
    }

    mu_eff = (1.0 - efrac) * (1.0 - lrfrac) * mu_hat_[ie][ir] +
	efrac         * (1.0 - lrfrac) * mu_hat_[ie+1][ir] +
	efrac         * lrfrac         * mu_hat_[ie+1][ir+1] +
	(1.0 - efrac) * lrfrac         * mu_hat_[ie][ir+1];

    k_eff = (1.0 - efrac) * (1.0 - lrfrac) * k_hat_[ie][ir] +
	efrac         * (1.0 - lrfrac) * k_hat_[ie+1][ir] +
	efrac         * lrfrac         * k_hat_[ie+1][ir+1] +
	(1.0 - efrac) * lrfrac         * k_hat_[ie][ir+1];
   

    Q.mu = mu_eff;
    Q.k[0] = k_eff;

    return (SUCCESS);
}

int
Look_up_table::
s_eval_diffusion_coefficients(Gas_data &Q)
{
    // These have no meaning for an equilibrium gas.
    Q.D_AB[0][0] = 0.0;
    return (SUCCESS);
}

double
Look_up_table::
s_molecular_weight(int isp)
{
    // This method is not very meaningful for an equilibrium
    // gas.  The molecular weight is best obtained from
    // the mixture molecular weight methods which IS a function
    // of gas composition and thermodynamic state, however,
    // there are times when a value from this function makes
    // other code simpler, in that the doesn't have to treat
    // the look-up gas specially.
    // cout << "Caution: calling s_molecular_weight for LUT species." << endl;
    if ( isp != 0 ) {
	throw runtime_error("LUT gas: should not be looking up isp != 0");
    }
    int ie = 0; // coldest
    int ir = irsteps_ - 1; // quite dense 
    double Rgas = R_hat_[ie][ir]; // J/kg/deg-K
    double M = PC_R_u_kmol / Rgas;
    return M;
}

double
Look_up_table::
s_internal_energy(const Gas_data &Q, int isp)
{
    // This method should never be called expecting quality data
    // because the LUT gas doesn not keep the species information.
    // This implementation is here to keep C++ happy that
    // all of the methods are implemented as required,
    // and there may be times when having this function
    // return something reasonable may make other code
    // simpler because it doesn't have to treat the
    // LUT gas specially.
    // cout << "Caution: calling s_internal_energy for LUT species." << endl;
    if ( isp != 0 ) {
	throw runtime_error("LUT gas: should not be looking up isp != 0");
    }
    // Finally, we assume that the thermodynamic state is current.
    return Q.e[0];
}

double
Look_up_table::
s_enthalpy(const Gas_data &Q, int isp)
{
    // This method assumes that the internal energy,
    // pressure and density are up-to-date in the
    // gas_data struct. Then enthalpy is computed
    // from definition.
    // cout << "Caution: calling s_enthalpy for LUT species." << endl;
    if ( isp != 0 ) {
	throw runtime_error("LUT gas: should not be looking up isp != 0");
    }
    double h = Q.e[0] + Q.p/Q.rho;
    return h;
}

double
Look_up_table::
s_entropy(const Gas_data &Q, int isp)
{
    if ( isp != 0 ) {
	throw runtime_error("LUT gas: should not be looking up isp != 0");
    }
    double s_eff;
    if ( with_entropy ) {
	int ir, ie;
	double lrfrac, efrac;
	if ( determine_interpolants(Q, ir, ie, lrfrac, efrac) != SUCCESS ) {
	    cout << "Bailing out!\n";
	    exit(1);
	}
	s_eff  = (1.0 - efrac) * (1.0 - lrfrac) * s_[ie][ir] +
	    efrac         * (1.0 - lrfrac) * s_[ie+1][ir] +
	    efrac         * lrfrac         * s_[ie+1][ir+1] +
	    (1.0 - efrac) * lrfrac         * s_[ie][ir+1];
    } else {
	// Without having the entropy recorded as part of the original table,
	// the next best is to use a model of an ideal gas.
	cout << "Caution: calling s_entropy for LUT species without tabular data." << endl;
	int ie = 0; // coldest
	int ir = irsteps_ - 1; // quite dense 
	double R = R_hat_[ie][ir]; // J/kg/deg-K
	double Cp = R + Cv_hat_[ie][ir];
	constexpr double T1 = 300.0; // degrees K
	constexpr double p1 = 100.0e3; // Pa
	s_eff = Cp * log(Q.T[0]/T1) - R * log(Q.p/p1);
    } 
    return s_eff;
}

double
Look_up_table::
s_dedT_const_v(const Gas_data &Q, int &status)
{
    double efrac, lrfrac;
    int    ir, ie;
    double Cv_actual;

    if ( determine_interpolants(Q, ir, ie, lrfrac, efrac) != SUCCESS ) {
	cout << "Bailing out!\n";
	exit(1);
    }

    Cv_actual = (1.0 - efrac) * (1.0 - lrfrac) * Cv_[ie][ir] +
	efrac         * (1.0 - lrfrac) * Cv_[ie+1][ir] +
	efrac         * lrfrac         * Cv_[ie+1][ir+1] +
	(1.0 - efrac) * lrfrac         * Cv_[ie][ir+1];

    status = SUCCESS;
    return Cv_actual;

}

double
Look_up_table::
s_dhdT_const_p(const Gas_data &Q, int &status)
{
    double efrac, lrfrac;
    int    ir, ie;
    double Cv_actual, R_eff;

    if ( determine_interpolants(Q, ir, ie, lrfrac, efrac) != SUCCESS ) {
	cout << "Bailing out!\n";
	exit(1);
    }

    Cv_actual = (1.0 - efrac) * (1.0 - lrfrac) * Cv_[ie][ir] +
	efrac         * (1.0 - lrfrac) * Cv_[ie+1][ir] +
	efrac         * lrfrac         * Cv_[ie+1][ir+1] +
	(1.0 - efrac) * lrfrac         * Cv_[ie][ir+1];

    R_eff  = (1.0 - efrac) * (1.0 - lrfrac) * R_hat_[ie][ir] +
	efrac         * (1.0 - lrfrac) * R_hat_[ie+1][ir] +
	efrac         * lrfrac         * R_hat_[ie+1][ir+1] +
	(1.0 - efrac) * lrfrac         * R_hat_[ie][ir+1];

    status = SUCCESS;
    return (Cv_actual + R_eff);

}

double
Look_up_table::
s_gas_constant(const Gas_data &Q, int &status)
{
    double efrac, lrfrac;
    int    ir, ie;
    double R_eff;

    if ( determine_interpolants(Q, ir, ie, lrfrac, efrac) != SUCCESS ) {
	cout << "Bailing out!\n";
	exit(1);
    }

    R_eff  = (1.0 - efrac) * (1.0 - lrfrac) * R_hat_[ie][ir] +
	efrac         * (1.0 - lrfrac) * R_hat_[ie+1][ir] +
	efrac         * lrfrac         * R_hat_[ie+1][ir+1] +
	(1.0 - efrac) * lrfrac         * R_hat_[ie][ir+1];

    status = SUCCESS;
    return R_eff;
}

Gas_model* create_look_up_table_gas_model(const string cfile)
{
    return new Look_up_table(cfile);
}
