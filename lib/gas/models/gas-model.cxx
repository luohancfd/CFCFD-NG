/** \file gas-model.cxx
 *  \ingroup gas
 *
 *  \author Rowan J Gollan
 *  \version 03-Jul-2008
 **/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <string>
#include <stdexcept>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "../../nm/source/Richardson_extrapolation.hh"
#include "gas-model.hh"

// When the object factory is implemented, the following
// headers should be removed.
#include "user-defined-gas-model.hh"
#include "composite-gas-model.hh"
#include "look-up-table.hh"
#include "LUT-plus-composite-gas-model.hh"
#include "REFPROP-gas-model.hh"

using namespace std;

Gas_model::
Gas_model() {}

Gas_model::
Gas_model(string cfile)
{
    // In general, assume that a gas model is NOT compatible with the 
    // finite-rate chemistry module
    set_reaction_compatibility(false);
    // For the cases where it is, SEE: composite-gas-model.cxx where
    // the thermal behaviour model is set.

    lua_State *L = initialise_lua_State();
    if( luaL_dofile(L, cfile.c_str()) != 0 ) {
	ostringstream ost;
	ost << "Gas_model():\n";
	ost << "Error in gas model input file: " << cfile << endl;
	input_error(ost);
    }

    // Add a reverse mapping to species table.
    lua_getglobal(L, "species");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Gas_model():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }
    int nsp = lua_objlen(L, -1);
    s_names_.resize(nsp);
    for ( int isp = 1; isp <= nsp; ++isp ) {
	lua_rawgeti(L, 1, isp);
	const char* species = luaL_checkstring(L, -1); lua_pop(L, 1);
	s_names_[isp-1] = string(species);
    }

    M_.resize(nsp);
    charge_.resize(nsp);
    atomic_constituents_.resize(nsp);
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1); // A Lua list is offset one from the C++ vector index
	const char* sp = luaL_checkstring(L, -1);
	lua_pop(L, 1);

	// Now bring the specific species table to TOS
	lua_getglobal(L, sp);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Gas_model::Gas_model()\n";
	    ost << "Error locating information table for species: " << sp << endl;
	    input_error(ost);
	}

	M_[isp] = get_positive_value(L, -1, "M");
	try {
	    charge_[isp] = get_int(L, -1, "charge");
	}
	catch ( runtime_error &e ) {
	    cout << e.what();
	    cout << "No charge information available in input file for species: " << sp << endl;
	    cout << "Applying default value of 0." << endl;
	    charge_[isp] = 0;
	}
	    
	// Bring atomic constituents table to TOS
	lua_getfield(L, -1, "atomic_constituents");
	if ( lua_istable(L, -1) ) {
	    int t = lua_gettop(L);
	    // Start table traversal by using "nil" as a key
	    lua_pushnil(L);
	    while ( lua_next(L, t) != 0 ) {
		/* uses 'key' (at index -2) and 'value' (at index -1) */
		string atom = luaL_checkstring(L, -2);
		int n_atom = luaL_checkint(L, -1);
		atomic_constituents_[isp].insert(pair<string, int>(atom, n_atom));
		/* removes 'value'; keeps 'key' for next iteration */
		lua_pop(L, 1);
	    }
	}
	else {
	    cout << "No table of atomic constituents found in input file for species: " << sp << endl;
	    cout << "The atomic constituents of this species will not be available." << endl;
	}
	lua_pop(L, 1); // pop "atomic_constituents" off stack

	lua_pop(L, 1); // pop "sp" off stack
    }
    lua_pop(L, 1); // pop "species" off stack
    lua_close(L);
}

Gas_model::
~Gas_model() {}

int
Gas_model::
number_of_values_in_gas_data_copy() const
{
    int nv = 5; // rho, p, p_e, a, mu
    nv += nsp_; // array of M
    nv += nsp_; // array of massf
    nv += nsp_*nsp_; // matrix of D_AB
    nv += 3*nmodes_; // array of e, T, k
    return nv;
}

int
Gas_model::
no_atoms_of(string atom, int isp)
{
    map<string,int>::iterator it = atomic_constituents_[isp].find(atom);
    if ( it != atomic_constituents_[isp].end() ) {
	return it->second;
    }
    // else
    return 0;
}

void
Gas_model::
atomic_constituents(int isp, map<string, int> &m)
{
    m = atomic_constituents_[isp];
}

#define MAX_STEPS 30
#define MAX_RELATIVE_STEP 0.1

int
Gas_model::
s_eval_thermo_state_pT(Gas_data &Q)
{
    double drho, rho_old, rho_new, e_old, e_new, de;
    double drho_sign, de_sign;
    double Cv_eff, R_eff, T_old;
    double fp_old, fT_old, fp_new, fT_new;
    double dfp_drho, dfT_drho, dfp_de, dfT_de, det;
    int converged, count;

    double p_given = Q.p;
    double T_given = Q.T[0];
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect 
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    double fp_tol = 1.0e-6 * p_given;
    double fT_tol = 1.0e-6 * T_given;
    double fp_tol_fail = 0.02 * p_given;
    double fT_tol_fail = 0.02 * T_given;

    // Get an idea of the gas properties by calling the original
    // equation of state with some dummy values for density
    // and internal energy.
    Q.rho = 1.0; // kg/m**3 
    Q.e[0] = 2.0e5; // J/kg 
    eval_thermo_state_rhoe(Q);
    
    T_old = Q.T[0];
    R_eff = Q.p / (Q.rho * T_old);
    de = 0.01 * Q.e[0];
    Q.e[0] += de;
    if( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	cout << "eval_thermo_state_pT():\n";
	cout << "    Duff call to eval_thermo_state_rhoe, starting guess 1.\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }

    Cv_eff = de / (Q.T[0] - T_old);
    // Now, get a better guess for the appropriate density and
    // internal energy.
    e_old = Q.e[0] + (T_given - Q.T[0]) * Cv_eff;
    rho_old = p_given / (R_eff * T_given);

    // Evaluate state variables using this guess.
    Q.rho = rho_old;
    Q.e[0] = e_old;
    if ( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	cout << "eval_thermo_state_pT:\n";
	cout << "    Duff call to eval_thermo_state_rhoe, starting guess 2.\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    fp_old = p_given - Q.p;
    fT_old = T_given - Q.T[0];
    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fp_old) < fp_tol) && (fabs(fT_old) < fT_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
	// Perturb first dimension to get derivatives.
	rho_new = rho_old * 1.001;
	e_new = e_old;
	Q.rho = rho_new;
	Q.e[0] = e_new;
	if( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	    cout << "eval_thermo_state_pT():\n";
	    cout << "    Duff call to eval_thermo_state_rhoe, iteration " << count << ", call A\n";
	    Q.print_values();
	    return DUFF_EOS_ERROR;
	}
	fp_new = p_given - Q.p;
	fT_new = T_given - Q.T[0];
	dfp_drho = (fp_new - fp_old) / (rho_new - rho_old);
	dfT_drho = (fT_new - fT_old) / (rho_new - rho_old);
	// Perturb other dimension to get derivatives.
	rho_new = rho_old;
	e_new = e_old * 1.001;
	Q.rho = rho_new;
	Q.e[0] = e_new;
	if( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	    cout << "eval_thermo_state_pT():\n";
	    cout << "    Duff call to eval_thermo_state_rhoe, iteration " << count << ", call B\n";
	    Q.print_values();
	    return DUFF_EOS_ERROR;
	}

	fp_new = p_given - Q.p;
	fT_new = T_given - Q.T[0];
	dfp_de = (fp_new - fp_old) / (e_new - e_old);
	dfT_de = (fT_new - fT_old) / (e_new - e_old);

	det = dfp_drho * dfT_de - dfT_drho * dfp_de;
	if( fabs(det) < 1.0e-12 ) {
	    cout << "eval_thermo_state_pT():\n";
	    cout << "    Nearly zero determinant, det = " << det << endl;
	    return ZERO_DETERMINANT_ERROR;
	}
	drho = (-dfT_de * fp_old + dfp_de * fT_old) / det;
	de = (dfT_drho * fp_old - dfp_drho * fT_old) / det;
	if( fabs(drho) > MAX_RELATIVE_STEP * rho_old ) {
	    // move a little toward the goal 
	    drho_sign = (drho > 0.0 ? 1.0 : -1.0);
	    drho = drho_sign * MAX_RELATIVE_STEP * rho_old;
	} 
	if( fabs(de) > MAX_RELATIVE_STEP * e_old ) {
	    // move a little toward the goal
	    de_sign = (de > 0.0 ? 1.0 : -1.0);
	    de = de_sign * MAX_RELATIVE_STEP * e_old;
	} 
	rho_old += drho;
	e_old += de;
	// Make sure of consistent thermo state.
	Q.rho = rho_old;
	Q.e[0] = e_old;
	if( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	    cout << "eval_thermo_state_pT():\n";
	    cout << "    Duff call to eval_thermo_state_rhoe, iteration " << count << ", call C\n";
	    Q.print_values();
	    return DUFF_EOS_ERROR;
	}
	// Prepare for next iteration.
	fp_old = p_given - Q.p;
	fT_old = T_given - Q.T[0];
	converged = (fabs(fp_old) < fp_tol) && (fabs(fT_old) < fT_tol);
	++count;
        // cout << "eval_thermo_state_pT(): end of iteration " << count 
	//      << " rho= " << Q.rho << " e= " << Q.e[0] << endl;
	// cout << "    fp_old= " << fp_old << " fT_old= " << fT_old << endl;
    } // end while 
    if( count >= MAX_STEPS ) {
	cout << "eval_thermo_state_pT():\n";
	cout << "    Warning, iterations did not converge.\n";
	cout << "    p_given = " << p_given << ", T_given = " << T_given << endl;
	Q.print_values();
	return ITERATION_ERROR;
    }
    if( (fabs(fp_old) > fp_tol_fail) || (fabs(fT_old) > fT_tol_fail) ) {
	cout << "eval_thermo_state_pT():\n";
	cout << "    iterations failed badly.\n";
	cout << "    p_given = " << p_given << ", T_given = " << T_given << endl;
	Q.print_values();
	return ITERATION_ERROR;
    }
    // If we get to this point, assume that all is well.
    return SUCCESS;
}

int
Gas_model::
s_eval_thermo_state_rhoT(Gas_data &Q)
{
    double e_old, e_new, de, tmp, de_sign;
    double Cv_eff, T_old;
    double dfT_de, fT_old, fT_new;
    int converged, count;

    double rho_given = Q.rho;
    double T_given = Q.T[0];
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect 
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    double fT_tol = 1.0e-6 * T_given;
    double fT_tol_fail = 0.02 * T_given;

    // Get an idea of the gas properties by calling the original
    // equation of state with some dummy values for density
    // and internal energy.
    Q.rho = rho_given; // kg/m**3 
    Q.e[0] = 2.0e5; // J/kg 
    eval_thermo_state_rhoe(Q);
    T_old = Q.T[0];
    de = 0.01 * Q.e[0];
    Q.e[0] += de;
    if ( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	cout << "eval_thermo_state_rhoT():\n";
	cout << "    Duff call to eval_thermo_state_rhoe, starting guess 0.\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    Cv_eff = de / (Q.T[0] - T_old);
    // Now, get a better guess for the appropriate density and internal energy.
    e_old = Q.e[0] + (T_given - Q.T[0]) * Cv_eff;
    // Evaluate state variables using this guess.
    Q.rho = rho_given;
    Q.e[0] = e_old;
    if ( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	cout << "eval_thermo_state_rhoT:\n";
	cout << "    Duff call to eval_thermo_state_rhoT, starting guess 1.\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    fT_old = T_given - Q.T[0];
    // Perturb to get derivative.
    e_new = e_old * 1.001;
    Q.rho = rho_given;
    Q.e[0] = e_new;
    if ( eval_thermo_state_rhoe(Q) != 0 ) {
	printf("EOS_rhoT(): Duff call to EOS, starting guess 2\n");
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    fT_new = T_given - Q.T[0];
    dfT_de = (fT_new - fT_old) / (e_new - e_old);

    // At the start of iteration, we want *_old to be the best guess.
    if ( fabs(fT_new) < fabs(fT_old) ) {
	tmp = fT_new; fT_new = fT_old; fT_old = tmp;
	tmp = e_new; e_new = e_old; e_old = tmp;
    }
    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fT_old) < fT_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
	de = -fT_old / dfT_de;
	if ( fabs(de) > MAX_RELATIVE_STEP * e_old ) {
	    // move a little toward the goal 
	    de_sign = (de > 0.0 ? 1.0 : -1.0);
	    de = de_sign * MAX_RELATIVE_STEP * fabs(e_old);
	} 
	e_new = e_old + de;
	Q.rho = rho_given;
	Q.e[0] = e_new;
	if ( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	    cout << "eval_thermo_state_rhoT():\n";
	    cout << "    Duff call to eval_thermo_state_rhoe, iteration " << count << endl;
	    Q.print_values();
	    return DUFF_EOS_ERROR;
	}
	fT_new = T_given - Q.T[0];
	dfT_de = (fT_new - fT_old) / (e_new - e_old);
	// Prepare for the next iteration.
	++count;
	fT_old = fT_new;
	e_old = e_new;
	converged = fabs(fT_old) < fT_tol;
    }   // end while 
    // Ensure that we have the current data for all EOS variables.
    Q.rho = rho_given;
    Q.e[0] = e_old;
    if ( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	cout << "eval_thermo_state_rhoT()\n";
	cout << "    Duff call to EOS, after finishing iteration\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    if ( count >= MAX_STEPS ) {
	cout << "eval_thermo_state_rhoT():\n";
	cout << "    Warning, iterations did not converge.\n";
	cout << "    rho_given = " << rho_given << ", T_given = " << T_given << endl;
	Q.print_values();
	return ITERATION_ERROR;
    }
    if ( fabs(fT_old) > fT_tol_fail ) {
	cout << "eval_thermo_state_rhoT():\n";
	cout << "    iterations failed badly.\n";
	cout << "    rho_given = " << rho_given << ", T_given = " << T_given << endl;
	Q.print_values();
	return ITERATION_ERROR;
    }
    // If we get to this point, assume that all is well.
    return SUCCESS;
}

int
Gas_model::
s_eval_thermo_state_rhop(Gas_data &Q)
{
    double e_old, e_new, de, dedp, tmp, de_sign;
    double p_old;
    double dfp_de, fp_old, fp_new;
    int converged, count;

    double rho_given = Q.rho;
    double p_given = Q.p;
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect 
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    double fp_tol = 1.0e-6 * p_given;
    double fp_tol_fail = 0.02 * p_given;

    // Get an idea of the gas properties by calling the original
    // equation of state with some dummy values for density
    // and internal energy.
    Q.rho = rho_given; // kg/m**3
    Q.e[0] = 2.0e5; // J/kg 
    eval_thermo_state_rhoe(Q);
    p_old = Q.p;
    de = 0.01 * Q.e[0];
    Q.e[0] += de;
    if ( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	cout << "eval_thermo_state_rhop():\n";
	cout << "    Duff call to eval_thermo_state_rhoe, starting guess 0.\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    dedp = de / (Q.p - p_old);
    // Now, get a better guess for the appropriate internal energy.
    e_old = Q.e[0] + (p_given - Q.p) * dedp;
//     printf( "Initial guess e_old= %g dedp= %g\n", e_old, dedp );
    // Evaluate state variables using this guess.
    Q.rho = rho_given;
    Q.e[0] = e_old;
    if ( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	cout << "eval_thermo_state_rhop():\n";
	cout << "    Duff call to eval_thermo_state_rhoe, starting guess 1.\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    fp_old = p_given - Q.p;
    // Perturb to get derivative.
    e_new = e_old * 1.001;
    Q.rho = rho_given;
    Q.e[0] = e_new;
    if ( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	cout << "eval_thermo_state_rhop():\n";
	cout << "    Duff call to eval_thermo_state_rhoe, starting guess 2.\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    fp_new = p_given - Q.p;
    dfp_de = (fp_new - fp_old) / (e_new - e_old);

    // At the start of iteration, we want *_old to be the best guess.
    if ( fabs(fp_new) < fabs(fp_old) ) {
	tmp = fp_new; fp_new = fp_old; fp_old = tmp;
	tmp = e_new; e_new = e_old; e_old = tmp;
    }
    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fp_old) < fp_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
	de = -fp_old / dfp_de;
	if ( fabs(de) > MAX_RELATIVE_STEP * e_old ) {
	    // move a little toward the goal
	    de_sign = (de > 0.0 ? 1.0 : -1.0);
	    de = de_sign * MAX_RELATIVE_STEP * fabs(e_old);
	} 
	e_new = e_old + de;
	Q.rho = rho_given;
	Q.e[0] = e_new;
	if ( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	    cout << "eval_thermo_state_rhop():\n";
	    cout << "    Duff call to eval_thermo_state_rhoe, iteration " << count << endl;
	    Q.print_values();
	    return DUFF_EOS_ERROR;
	}
	fp_new = p_given - Q.p;
	dfp_de = (fp_new - fp_old) / (e_new - e_old);
	// Prepare for next iteration.
	++count;
	fp_old = fp_new;
	e_old = e_new;
	converged = fabs(fp_old) < fp_tol;
    }   // end while 
    // Ensure that we have the current data for all EOS variables.
    Q.rho = rho_given;
    Q.e[0] = e_old;
    if ( eval_thermo_state_rhoe(Q) != SUCCESS ) {
	cout << "eval_thermo_state_rhop():\n";
	cout << "    Duff call to eval_thermo_state_rhoe, after finishing iteration\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    if ( count >= MAX_STEPS ) {
	cout << "eval_thermo_state_rhop():\n";
	cout << "    Warning, iterations did not converge.\n";
	cout << "    rho_given = " << rho_given << ", p_given = " << p_given << endl;
	Q.print_values();
	return ITERATION_ERROR;
    }   // end if 
    if ( fabs(fp_old) > fp_tol_fail ) {
	cout << "eval_thermo_state_rhop():\n";
	cout << "    iterations failed badly.\n";
	cout << "    rho_given = " << rho_given << ", p_given = " << p_given << endl;
	Q.print_values();
	return ITERATION_ERROR;
    }
    return SUCCESS;
}

int
Gas_model::
s_eval_thermo_state_ps(Gas_data &Q, double p, double s)
{
    double T_old, T_new, dT, tmp, dT_sign;
    double dfs_dT, fs_old, fs_new;
    int converged, count;

    double s_given = s;
    double p_given = p;
    Q.p = p;
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect 
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    double fs_tol = 1.0e-6 * s_given;
    double fs_tol_fail = 0.02 * s_given;

    // Guess the thermo state assuming that T is a good guess.
    T_old = Q.T[0];
    if ( eval_thermo_state_pT(Q) != SUCCESS ) {
	cout << "eval_thermo_state_ps():\n";
	cout << "    Duff call to eval_thermo_state_pT, starting guess 0.\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    double s_old = mixture_entropy(Q);
    fs_old = s_given - s_old;
    // Perturb T to get a derivative estimate
    T_new = T_old * 1.001;
    Q.T[0] = T_new;
    if ( eval_thermo_state_pT(Q) != SUCCESS ) {
	cout << "eval_thermo_state_ps():\n";
	cout << "    Duff call to eval_thermo_state_pT, starting guess 1.\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    double s_new = mixture_entropy(Q);
    fs_new = s_given - s_new;
    dfs_dT = (fs_new - fs_old)/(T_new - T_old);
    // At the start of iteration, we want *_old to be the best guess.
    if ( fabs(fs_new) < fabs(fs_old) ) {
	tmp = fs_new; fs_new = fs_old; fs_old = tmp;
	tmp = s_new; s_new = s_old; s_old = tmp;
    }

    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fs_old) < fs_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
	dT = -fs_old / dfs_dT;
	if ( fabs(dT) > MAX_RELATIVE_STEP * T_old ) {
	    // move a little toward the goal
	    dT_sign = (dT > 0.0 ? 1.0 : -1.0);
	    dT = dT_sign * MAX_RELATIVE_STEP * fabs(T_old);
	} 
	T_new = T_old + dT;
	Q.T[0] = T_new;
	if ( eval_thermo_state_pT(Q) != SUCCESS ) {
	    cout << "eval_thermo_state_ps():\n";
	    cout << "    Duff call to eval_thermo_state_pT, iteration " << count << endl;
	    Q.print_values();
	    return DUFF_EOS_ERROR;
	}
	s_new = mixture_entropy(Q);
	fs_new = s_given - s_new;
	dfs_dT = (fs_new - fs_old) / (T_new - T_old);
	// Prepare for next iteration.
	++count;
	fs_old = fs_new;
	T_old = T_new;
	converged = (fabs(fs_old) < fs_tol);
    }   // end while 
    // Ensure that we have the current data for all EOS variables.
    Q.T[0] = T_old;
    if ( eval_thermo_state_pT(Q) != SUCCESS ) {
	cout << "eval_thermo_state_ps():\n";
	cout << "    Duff call to eval_thermo_state_pT, after finishing iteration\n";
	Q.print_values();
	return DUFF_EOS_ERROR;
    }
    if ( count >= MAX_STEPS ) {
	cout << "eval_thermo_state_ps():\n";
	cout << "    Warning, iterations did not converge.\n";
	cout << "    p_given = " << p_given << ", s_given = " << s_given << endl;
	Q.print_values();
	return ITERATION_ERROR;
    }   // end if 
    if ( fabs(fs_old) > fs_tol_fail ) {
	cout << "eval_thermo_state_ps():\n";
	cout << "    iterations failed badly.\n";
	cout << "    p_given = " << p_given << ", s_given = " << s_given << endl;
	Q.print_values();
	return ITERATION_ERROR;
    }
    return SUCCESS;
}

int
Gas_model::
s_eval_thermo_state_hs(Gas_data &Q, double h, double s)
{
    double dp, p_old, p_new, T_old, T_new, dT;
    double dp_sign, dT_sign;
    double fh_old, fs_old, fh_new, fs_new;
    double dfh_dp, dfs_dp, dfh_dT, dfs_dT, det;
    int converged, count;

    double h_given = h;
    double s_given = s;
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect 
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    double fh_tol = 1.0e-6 * h_given;
    double fs_tol = 1.0e-6 * s_given;
    double fh_tol_fail = 0.02 * h_given;
    double fs_tol_fail = 0.02 * s_given;

    // Use current gas state as guess
    p_old = Q.p;
    T_old = Q.T[0];
    double h_new = mixture_enthalpy(Q);
    double s_new = mixture_entropy(Q);
    fh_old = h_given - h_new;
    fs_old = s_given - s_new;

    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fh_old) < fh_tol) && (fabs(fs_old) < fs_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
	// Perturb first dimension to get derivatives.
	p_new = p_old * 1.001;
	T_new = T_old;
	Q.p = p_new;
	Q.T[0] = T_new;
	if ( eval_thermo_state_pT(Q) != SUCCESS ) {
	    cout << "eval_thermo_state_hs():\n";
	    cout << "    Duff call to eval_thermo_state_pT, iteration " << count << ", call A\n";
	    Q.print_values();
	    return DUFF_EOS_ERROR;
	}
	h_new = mixture_enthalpy(Q);
	s_new = mixture_entropy(Q);
	fh_new = h_given - h_new;
	fs_new = s_given - s_new;
	dfh_dp = (fh_new - fh_old) / (p_new - p_old);
	dfs_dp = (fs_new - fs_old) / (p_new - p_old);
	// Perturb other dimension to get derivatives.
	p_new = p_old;
	T_new = T_old * 1.001;
	Q.p = p_new;
	Q.T[0] = T_new;
	if( eval_thermo_state_pT(Q) != SUCCESS ) {
	    cout << "eval_thermo_state_hs():\n";
	    cout << "    Duff call to eval_thermo_state_pT, iteration " << count << ", call B\n";
	    Q.print_values();
	    return DUFF_EOS_ERROR;
	}
	h_new = mixture_enthalpy(Q);
	s_new = mixture_entropy(Q);
	fh_new = h_given - h_new;
	fs_new = s_given - s_new;
	dfh_dT = (fh_new - fh_old) / (T_new - T_old);
	dfs_dT = (fs_new - fs_old) / (T_new - T_old);

	det = dfh_dp * dfs_dT - dfs_dp * dfh_dT;
	if( fabs(det) < 1.0e-12 ) {
	    cout << "eval_thermo_state_hs():\n";
	    cout << "    Nearly zero determinant, det = " << det << endl;
	    return ZERO_DETERMINANT_ERROR;
	}
	dp = (-dfs_dT * fh_old + dfh_dT * fs_old) / det;
	dT = (dfs_dp * fh_old - dfh_dp * fs_old) / det;
	if( fabs(dp) > MAX_RELATIVE_STEP * p_old ) {
	    // move a little toward the goal 
	    dp_sign = (dp > 0.0 ? 1.0 : -1.0);
	    dp = dp_sign * MAX_RELATIVE_STEP * p_old;
	} 
	if( fabs(dT) > MAX_RELATIVE_STEP * T_old ) {
	    // move a little toward the goal
	    dT_sign = (dT > 0.0 ? 1.0 : -1.0);
	    dT = dT_sign * MAX_RELATIVE_STEP * T_old;
	} 
	p_old += dp;
	T_old += dT;
	// Make sure of consistent thermo state.
	Q.p = p_old;
	Q.T[0] = T_old;
	if( eval_thermo_state_pT(Q) != SUCCESS ) {
	    cout << "eval_thermo_state_hs():\n";
	    cout << "    Duff call to eval_thermo_state_pT, iteration " << count << ", call C\n";
	    Q.print_values();
	    return DUFF_EOS_ERROR;
	}
	h_new = mixture_enthalpy(Q);
	s_new = mixture_entropy(Q);
	// Prepare for next iteration.
	fh_old = h_given - h_new;
	fs_old = s_given - s_new;
	converged = (fabs(fh_old) < fh_tol) && (fabs(fs_old) < fs_tol);
	++count;
    } // end while 
    if( count >= MAX_STEPS ) {
	cout << "eval_thermo_state_hs():\n";
	cout << "    Warning, iterations did not converge.\n";
	cout << "    h_given = " << h_given << ", s_given = " << s_given << endl;
	Q.print_values();
	return ITERATION_ERROR;
    }
    if( (fabs(fh_old) > fh_tol_fail) || (fabs(fs_old) > fs_tol_fail) ) {
	cout << "eval_thermo_state_hs():\n";
	cout << "    iterations failed badly.\n";
	cout << "    h_given = " << h_given << ", s_given = " << s_given << endl;
	Q.print_values();
	return ITERATION_ERROR;
    }
    // If we get to this point, assume that all is well.
    return SUCCESS;
}


int
Gas_model::
s_eval_sound_speed(Gas_data &Q)
{
    // Reference:
    // Cengel and Boles (1998)
    // Thermodynamics: an Engineering Approach, 3rd edition
    // McGraw Hill
    // Equation 16-10 on p. 849
    
    // "frozen" sound speed
    cout << "eval_sound_speed()\n";
    int status, status1, status2;
    Q.a = sqrt( gamma(Q, status1)*dpdrho_const_T(Q, status2) );

    status = SUCCESS;
    if ( status1 != SUCCESS )
	status = status1;
    if ( status2 != SUCCESS )
	status = status2;
    return status;
}

double
Gas_model::
s_dTdp_const_rho(const Gas_data &Q, int &status)
{
    Gas_data Q_work(this);
    Q_work.copy_values_from(Q);
	
    const double tol = 1.0e-6;
    const int max_steps = 5;
    dTdp_functor dTdp_f(*this, Q_work);
    double p = Q.p;
    double h = 0.001 * Q.p;
    return R_extrap_deriv(dTdp_f, p, h, status, tol, max_steps);
}

double
Gas_model::
s_dTdrho_const_p(const Gas_data &Q, int &status)
{
    Gas_data Q_work(this);
    Q_work.copy_values_from(Q);

    const double tol = 1.0e-6;
    const int max_steps = 5;
    dTdrho_functor dTdrho_f(*this, Q_work);
    double rho = Q.rho;
    double h = 0.001 * Q.rho;
    return R_extrap_deriv(dTdrho_f, rho, h, status, tol, max_steps);
}

double
Gas_model::
s_dTdrho_const_s( const Gas_data &Q, int &status )
{
    return Q.T[0] / s_dedT_const_v( Q, status ) / Q.rho / Q.rho / s_dTdp_const_rho( Q, status );
}

double
Gas_model::
s_dpdrho_const_T(const Gas_data &Q, int &status)
{
    Gas_data Q_work(this);
    Q_work.copy_values_from(Q);

    const double tol = 1.0e-6;
    const int max_steps = 5;
    dpdrho_functor dpdrho_f(*this, Q_work);
    double rho = Q.rho;
    double h = 0.001 * Q.rho;
    return R_extrap_deriv(dpdrho_f, rho, h, status, tol, max_steps);
}

double
Gas_model::
s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
{
    Gas_data Q_work(this);
    Q_work.copy_values_from(Q);

    const double tol = 1.0e-6;
    const int max_steps = 5;
    dpdrho_i_functor dpdrho_i_f(*this, Q_work, isp);
    double rho_i = Q.rho * Q.massf[isp];
    double h = 0.001 * rho_i;
    return R_extrap_deriv(dpdrho_i_f, rho_i, h, status, tol, max_steps);
}

double
Gas_model::
s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
{
    // Default behaviour is to revert to inverse of net dTdp function
    // This requires only one thermal mode to be correct, so lets test for that
    if ( nmodes_!=1 || itm!=0 ) {
    	cout << "Gas_model::s_dpdT_i_const_rho()" << endl
    	     << "This function is only correct when one thermal mode is present." << endl
    	     << "Exiting program." << endl;
    	exit( BAD_INPUT_ERROR );
    }
    Gas_data Q_work(this);
    Q_work.copy_values_from(Q);
	
    const double tol = 1.0e-6;
    const int max_steps = 5;
    dTdp_functor dTdp_f(*this, Q_work);
    double p = Q.p;
    double h = 0.001 * Q.p;
    return 1.0/R_extrap_deriv(dTdp_f, p, h, status, tol, max_steps);
}

double
Gas_model::
s_dedT_const_v(const Gas_data &Q, int &status)
{
    Gas_data Q_work(this);
    Q_work.copy_values_from(Q);

    const double tol = 1.0e-6;
    const int max_steps = 5;
    dedT_functor dedT_f(*this, Q_work);
    double T = Q.T[0];
    double h = 0.001 * Q.T[0];
    return R_extrap_deriv(dedT_f, T, h, status, tol, max_steps);
}

double
Gas_model::
s_dhdT_const_p(const Gas_data &Q, int &status)
{
    Gas_data Q_work(this);
    Q_work.copy_values_from(Q);

    const double tol = 1.0e-6;
    const int max_steps = 5;
    dhdT_functor dhdT_f(*this, Q_work);
    double T = Q.T[0];
    double h = 0.001 * Q.T[0];
    return R_extrap_deriv(dhdT_f, T, h, status, tol, max_steps);
}

double
Gas_model::
s_gas_constant(const Gas_data &Q, int &status)
{
    status = SUCCESS;
    return calculate_gas_constant(Q.massf, M_);
}

double
Gas_model::
s_mixture_molecular_weight(const Gas_data &Q)
{
    return calculate_molecular_weight(Q.massf, M_);  
}

double
Gas_model::
s_molecular_weight(int isp)
{
    return M_[isp];
}

double
Gas_model::
s_modal_enthalpy( const Gas_data &Q, int isp, int itm )
{
    // NOTE: this function should never be used as the multi-temperature thermal
    //       models implement their own s_modal_enthalpy() function
    cout << "Gas_model::s_modal_enthalpy()" << endl
         << "Something has gone wrong as this function is not intended for use." << endl;
    exit( FAILURE );
}

double
Gas_model::
s_modal_Cv(Gas_data &Q, int itm )
{
    if ( itm!=0 || nmodes_!=1 ) {
	cout << "Gas_model::s_modal_Cv()" << endl
	     << "Something has gone wrong as this function is not intended for use" << endl
	     << "when itm != 0 or nmodes != 1" << endl;
	exit( FAILURE );
    }
    
    int status;
    return s_dedT_const_v(Q,status);
}

double
Gas_model::
Gibbs_free_energy(const Gas_data &Q, int isp)
{
    double h = s_enthalpy(Q, isp);
    double s = s_entropy(Q, isp);
    double g = h - Q.T[0]*s;
    return g;
}

int
Gas_model::
get_isp_from_species_name( string name )
{
    for ( int isp=0; isp<nsp_; ++isp )
	if ( s_names_[isp] == name ) return isp;
    
    cout << "Species with name: " << name << " is not part of this gas-model" << endl;
    exit(VALUE_ERROR);
}

int
Gas_model::
get_imode_from_mode_name(string name)
{
    for ( int imode = 0; imode < nmodes_; ++imode ) {
	if ( m_names_[imode] == name )
	    return imode;
    }
    // If we get this far, then haven't found the name
    cout << "Mode name '" << name << "' could not be found in this gas model.\n";
    exit(VALUE_ERROR);
}

int
Gas_model::
set_mole_fractions(Gas_data &Q,
		   vector<string> &sp, // species names to set
		   const vector<double> &X)  // corresponding mole fractions
{
    // since this is likely to be called on the user level
    // use assertions to check for mistakes
    
    // ASSERT(sp.size() == X.size());
    double sum = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
	sum += X[i];
    }
    // ASSERT((sum < 1.001) && (sum > 0.999));
    
    int isp;
    vector<double> molef(nsp_, 0.0);
    vector<double> massf(nsp_, 0.0);
    
    int j = 0;
    vector<string>::iterator spname = sp.begin();
    for (; spname != sp.end(); ++spname) {
	isp = get_isp_from_species_name(*spname);
	molef[isp] = X[j];
	//printf("molef[%i] (%s) = %g\n", isp, (*spname).c_str(), molef[isp]);
	++j;
    }
    
    // convert to massf
    convert_molef2massf(M_, molef, massf);
    for (int isp = 0; isp < nsp_; ++isp) {
	Q.massf[isp] = massf[isp];
    }
    return SUCCESS;
}

double
Gas_model::
mixture_internal_energy(const Gas_data &Q, 
		      double T) 
{
    double e = 0.0;
    if (T == 0.0) {
	for (size_t isp = 0; isp < (size_t)nsp_; ++isp) {
	    e += Q.massf[isp]*s_internal_energy(Q, isp);
	}
    } else {
	// make a copy of the gas and set the temperature
	Gas_data Q0(this);
	Q0.copy_values_from(Q);
	Q0.T[0] = T;
	eval_thermo_state_pT(Q0);

	for (size_t isp = 0; isp < (size_t)nsp_; ++isp) {
	    e += Q0.massf[isp]*s_internal_energy(Q0, isp);
	}
    }

    return e;
}

double
Gas_model::
mixture_enthalpy(const Gas_data &Q, 
	       double T)
{
    double h = 0.0;
    if (T == 0.0) {
	for (size_t isp = 0; isp < (size_t)nsp_; ++isp) {
	    h += Q.massf[isp]*s_enthalpy(Q, isp);
	}
    } else {
	// make a copy of the gas and set the temperature
	Gas_data Q0(this);
	Q0.copy_values_from(Q);
	Q0.T[0] = T;
	eval_thermo_state_pT(Q0);

	for (size_t isp = 0; isp < (size_t)nsp_; ++isp) {
	    h += Q0.massf[isp]*s_enthalpy(Q0, isp);
	}
    }
    
    return h;
}

double
Gas_model::
mixture_entropy( const Gas_data &Q )
{
    double s = 0.0;
    for ( int isp=0; isp<nsp_; ++isp )
    	s += Q.massf[isp]*s_entropy(Q,isp);
    
    return s;
}

Gas_model* 
create_gas_model(string cfile)
{
    // Later this should be implemented as an object factory.
    map<string, Gas_model* (*)(string)> gas_models;
    gas_models.insert(pair<string, Gas_model* (*)(string)>("user-defined", create_user_defined_gas_model));
    gas_models.insert(pair<string, Gas_model* (*)(string)>("composite gas", create_composite_gas_model));
    gas_models.insert(pair<string, Gas_model* (*)(string)>("look-up table", create_look_up_table_gas_model));
    gas_models.insert(pair<string, Gas_model* (*)(string)>("LUT-plus-composite", create_LUT_plus_composite_gas_model));
    gas_models.insert(pair<string, Gas_model* (*)(string)>("REFPROP", create_REFPROP_gas_model));

    lua_State *L = initialise_lua_State();

    if( do_gzfile(L, cfile) != 0 ) {
	ostringstream ost;
	ost << "create_gas_model():\n";
	ost << "Error in gas model input file: " << cfile << endl;
	input_error(ost);
    }

    lua_getglobal(L, "model");
    if ( !lua_isstring(L, -1) ) {
	ostringstream ost;
	ost << "create_gas_model():\n";
	ost << "Error in gas model input file: " << cfile << endl;
	ost << "The gas model has not been declared as a string, that is,\n";
	ost << "   model = 'model_type'\n";
	input_error(ost);
    }

    string model(luaL_checkstring(L, -1));

    if ( gas_models.find(model) == gas_models.end() ) {
	ostringstream ost;
	ost << "create_gas_model():\n";
	ost << "Error in gas model input file: " << cfile << endl;
	ost << "The gas model: " << model << " is unknown.\n";
	input_error(ost);
    }

    lua_close(L);
    
    return gas_models[model](cfile);
}

void call_gas_model_deconstructor( Gas_model *gm )
{
    gm->~Gas_model();
}

/*
int declare_model(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 || !lua_istable(L, 1) ) {
	ostringstream ost;
	ost << "Only one table argument is expected after the 'model' declaration.\n";
	input_error(ost);
    }

    lua_rawgeti(L, 1, 1);
    const char* model = luaL_checkstring(L, -1); lua_pop(L, 1);
    lua_pushstring(L, model);
    lua_setglobal(L, "model");
    return 0;
}
*/

int 
declare_species(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 || !lua_istable(L, 1) ) {
	ostringstream ost;
	ost << "Only one table argument is expected after the 'species' declaration.\n";
	input_error(ost);
    }

    int nsp = lua_objlen(L, 1);

    // Lua indexes a table from 1
    for ( int isp = 1; isp <= nsp; ++isp ) {
	lua_rawgeti(L, 1, isp);
	const char* species = luaL_checkstring(L, -1); lua_pop(L, 1);
	// Construct an empty table associated with the species
	lua_newtable(L);
	// Push some items into that empty table.
	// 1. an empty viscosity table
	//	lua_newtable(L);
	//	lua_setfield(L, -2, "viscosity");
	// 2. an empty thermal conductivity table.
	//	lua_newtable(L);
	//	lua_setfield(L, -2, "thermal_conductivity");
	// And give this main table the species name in global register.
	lua_setglobal(L, species);
    }

    // Set the species function to a table now that
    // we have used the function.
    lua_insert(L, 1);
    lua_setglobal(L, "species");

    return 0;
}


lua_State* 
initialise_lua_State()
{
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    // Setup 'species' function calls
    lua_pushcfunction(L, declare_species);
    lua_setglobal(L, "species");

    return L;
}

void 
convert_massf2molef(const std::vector<double> &massf,
		    const std::vector<double> &M,
		    std::vector<double> &molef)
{
    double M_mix = calculate_molecular_weight(massf, M);
    for (size_t i = 0; i < molef.size(); ++i) {
	molef[i] = massf[i]*M_mix/M[i];
    }
}

void
convert_molef2massf(const vector<double> &molef,
		    const vector<double> &M, 
		    vector<double> &massf)
{
    double M_mix = mole_average(molef, M);
    for (int i = 0; i < (int)massf.size(); ++i) {
	massf[i] = molef[i]*M[i]/M_mix;
    }
}

void 
convert_massf2conc(double rho,
		   const std::vector<double> &massf,
		   const std::vector<double> &M,
		   std::vector<double> &c)
{
    const double min_moles = 1.0e-30;
    for( size_t isp = 0; isp < massf.size(); ++isp ) {
	c[isp] = massf[isp] * rho / M[isp];
	if( c[isp] < min_moles ) // helps with chemistry
	    c[isp] = 0.0;
    }
}

void 
convert_conc2massf(double rho,
		   const std::vector<double> &c,
		   const std::vector<double> &M,
		   std::vector<double> &massf)
{
    const double min_mass_frac = 1.0e-50;
    for ( size_t isp = 0; isp < massf.size(); ++isp ) {
	massf[isp] = c[isp] * M[isp] / rho;
	if ( massf[isp] < min_mass_frac ) {
	    massf[isp] = 0.0;
	}
    }
}

void 
convert_conc2molef(double rho_bar,
		   const vector<double> &c,
		   vector<double> &molef)
{
    const double min_mole_frac = 1.0e-50;
    for ( size_t isp = 0; isp < molef.size(); ++isp ) {
	molef[isp] = c[isp] / rho_bar;
	if ( molef[isp] < min_mole_frac ) {
	    molef[isp] = 0.0;
	}
    }
}

// Python-friendly functions.
vector<double> 
convert_massf2molef(const vector<double> &massf,
		    const vector<double> &M)
{
    vector<double> molef(massf.size(), 0.0);
    convert_massf2molef(massf, M, molef);
    return molef;
}
vector<double>
convert_molef2massf(const vector<double> &molef,
		    const vector<double> &M)
{
    vector<double> massf(molef.size(), 0.0);
    convert_molef2massf(molef, M, massf);
    return massf;
}

vector<double> 
convert_massf2conc(double rho,
		   const vector<double> &massf,
		   const vector<double> &M)
{
    vector<double> c(massf.size(), 0.0);
    convert_massf2conc(rho, massf, M, c);
    return c;
}
vector<double> 
convert_conc2massf(double rho,
                   const vector<double> &c,
                   const vector<double> &M) 
{
    vector<double> massf(c.size(), 0.0);
    convert_conc2massf(rho, c, M, massf);
    return massf;
}

double
calculate_molecular_weight(const vector<double> &massf, 
			   const vector<double> &M)
{
    // returns mixture molecular weight in kg/mol

    // Reference:
    // Cengel and Boles (2002)
    // Thermodynamics: an Engineering Approach, 3rd edition
    // 4th Ed.
    // McGraw Hill
    // Equation 11-40 on p. 634
    return 1.0/mass_average_inverse(massf, M);
}

double
calculate_gas_constant(const vector<double> &massf,
		       const vector<double> &M)
{
    // returns mixture gas constant in J/kg.K

    // Reference:
    // Cengel and Boles (2002)
    // Thermodynamics: an Engineering Approach, 3rd edition
    // 4th Ed.
    // McGraw Hill
    // Equation 11-40 on p. 634
    return PC_R_u/calculate_molecular_weight(massf, M);
}
