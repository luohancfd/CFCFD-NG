// Author: PJ, DFP, RJG (see below)
// Date: 05-Dec-2008
// Place: Hampton, Virginia, USA
// Note:
//   This is a port of PJs implementation of the 
//   LUT + ideal gas mix model.  DFP updated this
//   previously for libgas2.  This implementation
//   fits with the new Gas_model class.
//

#include <string>
#include <sstream>
#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"

#include "physical_constants.hh"
#include "LUT-plus-composite-gas-model.hh"

using namespace std;

const double insignificant_massf = 0.0001;

LUT_plus_composite::
LUT_plus_composite(string cfile)
    : Gas_model()
{
    
    // Read lua file and populate LUT
    lua_State *L = initialise_lua_State();
    
    if ( do_gzfile(L, cfile) != 0 ) {
	ostringstream ost;
	ost << "LUT_plus_composite():\n";
	ost << "Error in LUT_plus_composite input file: " << cfile << endl;
	input_error(ost);
    }

    // All we get from file, at this stage, is the name
    // of the look-up-table file.
    string lut_file = get_string(L, LUA_GLOBALSINDEX, "lut_file");
    LUT_ = new Look_up_table(lut_file);
    if ( LUT_ == 0 ) {
	ostringstream ost;
	ost << "LUT_plus_composite():\n";
	ost << "Error initialising look-up table component of the\n";
	ost << "LUT_plus_composite gas model.\n";
	input_error(ost);
    }

    // We now let the composite gas model re-open the
    // config file and parse its own information.
    CGM_ = new Composite_gas_model(cfile);
    if ( CGM_ == 0 ) {
	ostringstream ost;
	ost << "LUT_plus_composite():\n";
	ost << "Error initialising composite gas component of the\n";
	ost << "LUT_plus_composite gas model.\n";
	input_error(ost);
    }

    // plus 1 for LUT
    set_number_of_species(CGM_->get_number_of_species() + 1);
    if ( CGM_->get_number_of_modes() > 1 ) {
	ostringstream ost;
	ost << "LUT_plus_composite():\n";
	ost << "Error initialising LUT_plus_composite gas model.\n";
	ost << "This model assumes only one thermal mode, but the\n";
	ost << "supplied composite gas component has " 
	    << CGM_->get_number_of_modes() << " modes.\n";
	input_error(ost);
    }
    set_number_of_modes(1);

    // Initialise space in working gas data structures
    Q_LUT_ = new Gas_data(LUT_);
    Q_CGM_ = new Gas_data(CGM_);

    s_names_.resize(nsp_);
    s_names_[0] = "LUT";

    lua_getglobal(L, "species");
    for ( int isp = 1; isp < nsp_; ++isp ) {
	lua_rawgeti(L, 1, isp);
	const char* species = luaL_checkstring(L, -1); lua_pop(L, 1);
	s_names_[isp] = string(species);
    }
    lua_pop(L, 1);
    lua_close(L);

}

LUT_plus_composite::
~LUT_plus_composite()
{
    delete LUT_;
    delete CGM_;
    delete Q_LUT_;
    delete Q_CGM_;
}

int
LUT_plus_composite::
s_eval_thermo_state_rhoe(Gas_data &Q)
{
    double f_LUT, f_CGM, f_sum; 
    double Cv_LUT, Cp_LUT, R_LUT;
    double T_CGM, Cv_CGM, Cp_CGM, R_CGM;
    int result_flag, count;
    double factor_LUT, Cv_overall, Cv_hat;
    double x1, y1, x2, y2, slope, spare;
    double desired_increment, actual_increment;

    int nsp = get_number_of_species();

    if ( Q.rho <= 0.0 ) {
        cerr << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
	     << "    density= " << Q.rho << " is zero or negative on entry\n";
	cout << "    Original Q:" << endl;
	Q.print_values();
	return FAILURE;
    }

    // The LUT component is always the first.
    f_LUT = Q.massf[0];
    Q_LUT_->massf[0] = 1.0; // LUT pretends that there is no other species

    if ( f_LUT > (1.0 - insignificant_massf) ) {
        // The LUT gas completely dominates. 
        // Ignore the other components and delegate computation of the gas properties.
	Q_LUT_->rho = Q.rho;
	Q_LUT_->e[0] = Q.e[0];
	result_flag = LUT_->eval_thermo_state_rhoe(*Q_LUT_);
	if ( result_flag == SUCCESS ) {
	    Q.T[0] = Q_LUT_->T[0];
	    if ( Q.T[0] < gas_Tmin() ) {
		cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
		     << "    after delegating to look-up-table:" << endl;
		cout << "    Bad temperature." << endl;
		cout << "    Original Q:" << endl;
		Q.print_values();
		cout << "    Q_LUT:" << endl;
		Q_LUT_->print_values();
		return BAD_TEMPERATURE_ERROR;
	    }
	    Q.p = Q_LUT_->p;
	    Q.a = Q_LUT_->a;
	} else {
	    cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
		 << "    after after calling look-up-table:" << endl;
	    cout << "    Bad return flag." << endl;
	    cout << "    Original Q:" << endl;
	    Q.print_values();
	}
	return result_flag;
    }

    // The ideal gas components make up the rest.
    f_CGM = 1.0 - f_LUT;
    f_sum = 0.0;
    for ( int i = 1; i < nsp; ++i ) {
	Q_CGM_->massf[i-1] = Q.massf[i] / f_CGM;
	f_sum += Q_CGM_->massf[i-1];
    }
    if ( fabs(f_sum - 1.0) > 1.0e-3 ) {
	cout << "LUT_plus_composite::eval_thermo_state_rhoe(): f_sum=" << f_sum << endl; 
	cout << "    Original Q:" << endl;
	Q.print_values();
	return MASS_FRACTION_ERROR;
    }

    if ( f_LUT < insignificant_massf ) {
        // The ideal-gas components dominate, just use those values and return.
	Q_CGM_->rho = Q.rho;
	Q_CGM_->e[0] = Q.e[0];
	result_flag = CGM_->eval_thermo_state_rhoe(*Q_CGM_);
	if ( result_flag == SUCCESS ) {
	    Q.T[0] = Q_CGM_->T[0];
	    if ( Q.T[0] < gas_Tmin() ) {
		cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
		     << "    after delegating to composite gas:" << endl;
		cout << "    Bad temperature." << endl;
		cout << "    Original Q:" << endl;
		Q.print_values();
		exit(1);
		//return BAD_TEMPERATURE_ERROR;
	    }
	    Q.p = Q_CGM_->p;
	    Q.a = Q_CGM_->a;
	} else {
	    cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
		 << "    after after delegating to composite gas:" << endl;
	    cout << "    Bad return flag." << endl;
	    cout << "    Original Q:" << endl;
	    Q.print_values();
	}
	return result_flag;
    }

    // At this point, we are continuing with a mixture calculation.
    // The central idea is that we can share the internal energy between the LUT gas
    // and the composite-gas mixture, however, we must make the split such that they
    // have the same static temperatures.

    Q_CGM_->rho = Q.rho * f_CGM;
    Q_CGM_->e[0] = Q.e[0];  // Assume the same values for LUT and ideal-gas-mix
    result_flag = CGM_->eval_thermo_state_rhoe(*Q_CGM_);
    if ( result_flag != SUCCESS ) {
	cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
	     << "    after calling composite gas in mixture calc:" << endl;
	cout << "    Bad return flag." << endl;
	cout << "    Q_CGM_ is:" << endl;
	Q_CGM_->print_values();
	cout << "    Original Q:" << endl;
	Q.print_values();
	return result_flag;
    }
    Cv_CGM = CGM_->Cv(*Q_CGM_, result_flag);
    R_CGM = CGM_->R(*Q_CGM_, result_flag);

    // Take a guess at the value of the internal energy of the LUT gas.
    factor_LUT = 1.0;
    Q_LUT_->e[0] = Q.e[0] * factor_LUT;
    Q_LUT_->rho = Q.rho * f_LUT;
    result_flag = LUT_->eval_thermo_state_rhoe(*Q_LUT_);
    if ( result_flag != SUCCESS ) {
	cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
	     << "    after calling look-up-table in mixture calc:" << endl;
	cout << "    Bad return flag." << endl;
	cout << "    Q_LUT_ is:" << endl;
	Q_LUT_->print_values();
	cout << "    Original Q:" << endl;
	Q.print_values();
	return result_flag;
    }
    T_CGM = (Q.e[0] - f_LUT * Q_LUT_->e[0]) / (f_CGM * Cv_CGM);

    if ( Q_LUT_->T[0] < gas_Tmin() || T_CGM < gas_Tmin() ) {
	cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
	     << "    Low temperature at initial guess:\n";
	cout << "    rho= " << Q.rho << ", e[0]= " << Q.e[0] 
	     << ", T_LUT= " << Q_LUT_->T[0] 
	     << ", T_CGM= " << T_CGM << endl;
	cout << "    Q_LUT_ is:" << endl;
	Q_LUT_->print_values();
	cout << "    Original Q:" << endl;
	Q.print_values();
	return BAD_TEMPERATURE_ERROR;
    }
    x1 = factor_LUT;
    y1 = Q_LUT_->T[0] - T_CGM;
    // With Cv_hat (e/T), we can work out a better split of the energy.
    Cv_hat = Q_LUT_->e[0] / Q_LUT_->T[0];
    Cv_overall = f_LUT * Cv_hat + f_CGM * Cv_CGM;
    factor_LUT = Cv_hat / Cv_overall;
    Q_LUT_->e[0] = Q.e[0] * factor_LUT;
    result_flag = LUT_->eval_thermo_state_rhoe(*Q_LUT_);
    if ( result_flag != SUCCESS ) {
	cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
	     << "    after 2nd call to look-up-table in mixture calc:" << endl;
	cout << "    Bad return flag." << endl;
	cout << "    Q_LUT_ is:" << endl;
	Q_LUT_->print_values();
	cout << "    Original Q:" << endl;
	Q.print_values();
	return result_flag;
    }
    T_CGM = (Q.e[0] - f_LUT * Q_LUT_->e[0]) / (f_CGM * Cv_CGM);
    x2 = factor_LUT;
    y2 = Q_LUT_->T[0] - T_CGM;

    // Set up secant iterations.
    if ( fabs(y2) > fabs(y1) ) {
	// Make x2 the best guess so far.
	spare = x2; x2 = x1; x1 = spare;
	spare = y2; y2 = y1; y1 = spare;
    }
    count = 0;
    const int max_secant_iterations = 15;
    const int T_diff_tol = 0.01; // degrees K
    while ( fabs(y2) > T_diff_tol && count < max_secant_iterations ) {
	// Improve the estimate of the energy fraction.
	++count;
	if ( fabs(x2 - x1) < 1.0e-12 ) {
	    // cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
	    //      << "    Very small increment in energy fraction: " 
	    //      << x2 - x1 << endl;
	    break;  // This is OK, so just return with the current guess.
	} else {
	    slope = (y2 - y1) / (x2 - x1);
	    if ( fabs(slope) < 1.0e-12 ) {
		cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
		     << "    Almost zero slope: " << slope << endl;
		break;
	    }
	}
	desired_increment = -y2/slope;
	if ( fabs(desired_increment) > 0.2 ) {
	    // If we are far from the solution, assume that the direction of move
	    // is correct and move a bit that way.
	    actual_increment = (desired_increment > 0.0) ? 0.2 : -0.2;
	} else {
	    actual_increment = desired_increment;
	}
	x1 = x2 + actual_increment;
	factor_LUT = x1; 
	Q_LUT_->e[0] = Q.e[0] * factor_LUT;
	result_flag = LUT_->eval_thermo_state_rhoe(*Q_LUT_);
	if ( result_flag != SUCCESS ) {
	    cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
		 << "    after after calling look-up-table in secant method:"
		 << endl;
	    cout << "    Bad return flag." << endl;
	    cout << "    Q_LUT_ is:" << endl;
	    Q_LUT_->print_values();
	    cout << "    Original Q:" << endl;
	    Q.print_values();
	    return result_flag;
	}
	T_CGM = (Q.e[0] - f_LUT * Q_LUT_->e[0]) / (f_CGM * Cv_CGM);
	y1 = Q_LUT_->T[0] - T_CGM;
	if ( fabs(y2) > fabs(y1) ) {
	    // Make x2 the best guess so far.
	    spare = x2; x2 = x1; x1 = spare;
	    spare = y2; y2 = y1; y1 = spare;
	}
    } //end while loop
	
    factor_LUT = x2; 
    Q_LUT_->e[0] = Q.e[0] * factor_LUT;
    result_flag = LUT_->eval_thermo_state_rhoe(*Q_LUT_);
    if ( result_flag != SUCCESS ) {
	cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
	     << "    after final call to look-up table:\n";
	cout << "    Bad return flag." << endl;
	cout << "    Q_LUT_ is:" << endl;
	Q_LUT_->print_values();
	cout << "    Original Q:" << endl;
	Q.print_values();
	return result_flag;
    }
    T_CGM = (Q.e[0] - f_LUT * Q_LUT_->e[0]) / (f_CGM * Cv_CGM);

    // We settled on a suitable split of energies.
    // Now we need to set the other thermo properties of the LUT+composite mixture.
    Q.p = Q_LUT_->p + f_CGM * Q.rho * R_CGM * T_CGM;
    Q.T[0] = Q_LUT_->T[0];
    // All of the following is to get an approximate sound speed.
    int status1, status2, status3, status4, status5, status6;
    Cv_LUT = LUT_->Cv(*Q_LUT_, status1);
    Cv_CGM = CGM_->Cv(*Q_CGM_, status2);
    double Cvm = f_LUT*Cv_LUT + f_CGM*Cv_CGM;
    Cp_LUT = LUT_->Cp(*Q_LUT_, status3);
    Cp_CGM = CGM_->Cp(*Q_CGM_, status4);
    double Cpm = f_LUT*Cp_LUT + f_CGM*Cp_CGM;
    R_LUT = LUT_->R(*Q_LUT_, status5);
    R_CGM = CGM_->R(*Q_CGM_, status6);
    double Rm = f_LUT*R_LUT + f_CGM*R_CGM;
    if ( status1 == SUCCESS && status2 == SUCCESS && status3 == SUCCESS &&
	 status4 == SUCCESS && status5 == SUCCESS && status6 == SUCCESS ) {
	result_flag = SUCCESS;
    } else {
	result_flag = FAILURE;
    }
    if ( result_flag != SUCCESS ) {
	cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n";
	cout << "    There was a problem evaluating Cv, Cp, R.\n";
	cout << "    The gas data for Q_LUT_ is: \n";
	Q.print_values();
	cout << "    The gas data for Q_LUT_ is: \n";
	Q.print_values();
	return result_flag;
    }
    double g = Cpm / Cvm;
    Q.a = sqrt(g * Rm * Q.T[0]);

    const double moderately_small_T_diff = 5.0; // degrees K
    if ( count >= max_secant_iterations && fabs(y2) > moderately_small_T_diff) {
	cout << "LUT_plus_composite::eval_thermo_state_rhoe():\n"
	     << "    Secant iteration did not converge." << endl;
	cout << "    Original Q:" << endl;
	Q.print_values();
	cout << "    f_LUT= " << f_LUT << ", f_CGM=" << f_CGM << endl;
	cout << "    Q_CGM_:" << endl;
	Q_CGM_->print_values();
	cout << "    Q_LUT_:" << endl;
	Q_LUT_->print_values();
	return ITERATION_ERROR;
    }

    // If we get to this point, assume that all was OK.
    return SUCCESS;
}

int
LUT_plus_composite::
s_eval_transport_coefficients(Gas_data &Q)
{
    int result_flag;
    int nsp = get_number_of_species();
    // First ensure that Q_LUT properties are set
    double f_LUT = Q.massf[0];
    Q_LUT_->massf[0] = 1.0; // LUT pretends that there is no other species
    Q_LUT_->rho = Q.rho;
    Q_LUT_->e[0] = Q.e[0];

    result_flag = LUT_->eval_transport_coefficients(*Q_LUT_);
    if (result_flag != SUCCESS ) {
	cout << "LUT_plus_composite::s_eval_transport_coeffcients()\n";
	cout << "    There is a problem evaluating the transport\n";
	cout << "    coefficients for the LUT model.\n";
	Q.print_values();
	return result_flag;
    }

    if ( f_LUT > (1.0 - insignificant_massf) ) {
	// The LUT gas dominates.
	Q.mu = Q_LUT_->mu;
	Q.k[0] = Q_LUT_->k[0];
	return result_flag;
    }

    Q_CGM_->T[0] = Q.T[0];
    Q_CGM_->p = Q.p;
    double f_CGM = 1.0 - f_LUT;
    double f_sum = 0.0;
    for ( int i = 1; i < nsp; ++i ) {
	Q_CGM_->massf[i-1] = Q.massf[i] / f_CGM;
	f_sum += Q_CGM_->massf[i-1];
    }
    if ( fabs(f_sum - 1.0) > 1.0e-3 ) {
	cout << "LUT_plus_composite::s_eval_transport_coefficients():\n"
	     << "    f_sum=" << f_sum << endl; 
	cout << "    Original Q:" << endl;
	Q.print_values();
	return MASS_FRACTION_ERROR;
    }

    result_flag = CGM_->eval_transport_coefficients(*Q_CGM_);
    if (result_flag != SUCCESS ) {
	cout << "LUT_plus_composite::s_eval_transport_coeffcients()\n";
	cout << "    There is a problem evaluating the transport\n";
	cout << "    coefficients for the composite gas model.\n";
	Q.print_values();
	return result_flag;
    }

    if ( f_LUT < insignificant_massf ) {
	// The CGM dominates.
	Q.mu = Q_CGM_->mu;
	Q.k[0] = Q_CGM_->k[0];
	return result_flag;
    }
    
    // Presently, a simple mass-weighted average of
    // the transport coefficients.
    // -- FIX ME --
    Q.mu = f_LUT*Q_LUT_->mu + f_CGM*Q_CGM_->mu;
    Q.k[0] = f_LUT*Q_LUT_->k[0] + f_CGM*Q_CGM_->k[0];

    return SUCCESS;
}

int
LUT_plus_composite::
s_eval_diffusion_coefficients(Gas_data &Q)
{
    cout << "LUT_plus_composite::s_eval_diffusion_coefficients()\n";
    cout << "Diffusion coefficients have no meaning for\n";
    cout << "this gas model.\n";
    return SUCCESS;
}

double
LUT_plus_composite::
s_molecular_weight(int isp)
{
    int status_flag;
    double R;
    if ( isp == 0 ) {
	// Dealing with LUT
	R = LUT_->R(*Q_LUT_, status_flag);
	if ( status_flag != SUCCESS ) {
	    cout << "LUT_plus_composite::s_molecular_weight()\n";
	    cout << "There was a problem computing the gas constant\n";
	    cout << "for the LUT model.\n";
	    cout << "Bailing out!\n";
	    exit(DUFF_EOS_ERROR);
	}
	return PC_R_u/R;
    }
    else {
	// Dealing with CGM
	return CGM_->molecular_weight(isp-1);
    }
}

double
LUT_plus_composite::
s_internal_energy(const Gas_data &Q, int isp)
{
    UNUSED_VARIABLE(isp);
    // This method should never be called.
    // This implementation is here to keep C++ happy that
    // all of the methods are implemented as required.
    cout << "LUT_plus_composite::s_internal_energy()\n";
    cout << "This function does nothing, it should not be\n";
    cout << "called.\n";
    return 0.0;
}

double
LUT_plus_composite::
s_enthalpy(const Gas_data &Q, int isp)
{
    UNUSED_VARIABLE(isp);
    // This method should never be called.
    // It is a utility for the finite-rate chemistry,
    // and an equilibrium gas should not be part of the
    // the finite-rate chemistry model by definition.
    // This implementation is here to keep C++ happy that
    // all of the methods are implemented as required.
    cout << "LUT_plus_composite::s_enthalpy()\n";
    cout << "This function does nothing, it should not be\n";
    cout << "called.  It is only needed for finite-rate chemistry\n";
    cout << "simulations.\n";
    return 0.0;
}

double
LUT_plus_composite::
s_entropy(const Gas_data &Q, int isp)
{
    UNUSED_VARIABLE(isp);
    // This method should never be called.
    // It is a utility for the finite-rate chemistry,
    // and an equilibrium gas should not be part of the
    // the finite-rate chemistry model by definition.
    // This implementation is here to keep C++ happy that
    // all of the methods are implemented as required.
    cout << "LUT_plus_composite::s_entropy()\n";
    cout << "This function does nothing, it should not be\n";
    cout << "called.  It is only needed for finite-rate chemistry\n";
    cout << "simulations.\n";
    return 0.0;
}

double
LUT_plus_composite::
s_dedT_const_v(const Gas_data &Q, int &status)
{
    int status1, status2;
    double f_LUT, f_CGM, f_sum;
    int nsp = get_number_of_species();

    f_LUT = Q.massf[0];
    Q_LUT_->massf[0] = 1.0; // LUT pretends that there is no other species
    Q_LUT_->rho = Q.rho;
    Q_LUT_->e[0] = Q.e[0];

    double Cv_LUT = LUT_->Cv(*Q_LUT_, status1);
    if ( f_LUT > (1.0 - insignificant_massf) ) {
	// The LUT gas dominates.
	status = status1;
	return Cv_LUT;
    }

    Q_CGM_->T[0] = Q.T[0];
    Q_CGM_->p = Q.p;
    Q_CGM_->rho = Q.rho;
    Q_CGM_->e[0] = Q.e[0];
    f_CGM = 1.0 - f_LUT;
    f_sum = 0.0;
    for ( int i = 1; i < nsp; ++i ) {
	Q_CGM_->massf[i-1] = Q.massf[i] / f_CGM;
	f_sum += Q_CGM_->massf[i-1];
    }
    if ( fabs(f_sum - 1.0) > 1.0e-3 ) {
	cout << "LUT_plus_composite::s_dedT_const_v(): f_sum=" << f_sum << endl; 
	cout << "    Original Q:" << endl;
	Q.print_values();
	status = MASS_FRACTION_ERROR;
	return 0.0;
    }
    
    double Cv_CGM = CGM_->Cv(*Q_CGM_, status2);
    if ( f_LUT < insignificant_massf ) {
	// The CGM dominates.
	status = status2;
	return Cv_CGM;
    }

    double Cv = f_LUT*Cv_LUT + f_CGM*Cv_CGM;
    
    status = SUCCESS;
    if ( status1 != SUCCESS )
	status = status1;
    if ( status2 != SUCCESS )
	status = status2;

    return Cv;
}

double
LUT_plus_composite::
s_dhdT_const_p(const Gas_data &Q, int &status)
{
    int status1, status2;
    double f_LUT, f_CGM, f_sum;
    int nsp = get_number_of_species();

    f_LUT = Q.massf[0];
    Q_LUT_->massf[0] = 1.0; // LUT pretends that there is no other species
    Q_LUT_->rho = Q.rho;
    Q_LUT_->e[0] = Q.e[0];
    double Cp_LUT = LUT_->Cp(*Q_LUT_, status1);

    if ( f_LUT > (1.0 - insignificant_massf) ) {
        // The LUT gas dominates
	status = status1;
	return Cp_LUT;
    }

    Q_CGM_->T[0] = Q.T[0];
    Q_CGM_->p = Q.p;
    Q_CGM_->rho = Q.rho;
    Q_CGM_->e[0] = Q.e[0];
    f_CGM = 1.0 - f_LUT;
    f_sum = 0.0;
    for ( int i = 1; i < nsp; ++i ) {
	Q_CGM_->massf[i-1] = Q.massf[i] / f_CGM;
	f_sum += Q_CGM_->massf[i-1];
    }
    if ( fabs(f_sum - 1.0) > 1.0e-3 ) {
	cout << "LUT_plus_composite::s_dhdT_const_p(): f_sum=" << f_sum << endl; 
	cout << "    Original Q:" << endl;
	Q.print_values();
	status = MASS_FRACTION_ERROR;
	return 0.0;
    }
    
    double Cp_CGM = CGM_->Cp(*Q_CGM_, status2);

    if ( f_LUT < insignificant_massf ) {
	// The CGM dominates
	status = status2;
	return Cp_CGM;
    }
    
    double Cp = f_LUT*Cp_LUT + f_CGM*Cp_CGM;
    
    status = SUCCESS;
    if ( status1 != SUCCESS )
	status = status1;
    if ( status2 != SUCCESS )
	status = status2;

    return Cp;
}

double
LUT_plus_composite::
s_gas_constant(const Gas_data &Q, int &status)
{
    int status1, status2;
    double f_LUT, f_CGM, f_sum;
    int nsp = get_number_of_species();

    f_LUT = Q.massf[0];
    Q_LUT_->massf[0] = 1.0; // LUT pretends that there is no other species
    Q_LUT_->rho = Q.rho;
    Q_LUT_->e[0] = Q.e[0];

    double R_LUT = LUT_->R(*Q_LUT_, status1);
    if ( f_LUT > (1.0 - insignificant_massf) ) {
	// The LUT gas dominates.
	status = status1;
	return R_LUT;
    }
    
    Q_CGM_->T[0] = Q.T[0];
    Q_CGM_->p = Q.p;
    Q_CGM_->rho = Q.rho;
    Q_CGM_->e[0] = Q.e[0];
    f_CGM = 1.0 - f_LUT;
    f_sum = 0.0;
    for ( int i = 1; i < nsp; ++i ) {
	Q_CGM_->massf[i-1] = Q.massf[i] / f_CGM;
	f_sum += Q_CGM_->massf[i-1];
    }
    if ( fabs(f_sum - 1.0) > 1.0e-3 ) {
	cout << "LUT_plus_composite::s_dhdT_const_p(): f_sum=" << f_sum << endl; 
	cout << "    Original Q:" << endl;
	Q.print_values();
	status = MASS_FRACTION_ERROR;
	return 0.0;
    }
    
    double R_CGM = CGM_->R(*Q_CGM_, status2);
    if ( f_LUT < insignificant_massf ) {
	status = status2;
	return R_CGM;
    }

    double R = f_LUT*R_LUT + f_CGM*R_CGM;
    
    status = SUCCESS;
    if ( status1 != SUCCESS )
	status = status1;
    if ( status2 != SUCCESS )
	status = status2;

    return R;
}

Gas_model* create_LUT_plus_composite_gas_model(string cfile)
{
    return new LUT_plus_composite(cfile);
}
