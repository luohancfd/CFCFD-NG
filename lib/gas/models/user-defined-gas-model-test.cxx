// Author: Rowan J. Gollan
// Version: 08-Jul-2008

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>

#include "gas_data.hh"
#include "user-defined-gas-model.hh"

using namespace std;

int main()
{
    
    const double tol = 1.0e-9;
    // Reference answer values
    const double rho = 1.2;
    const double T = 300.0;
    const double p = 103356.0;
    const double R = 287.1;
    const double e = 215325.0;
    const double a = 347.2491900638;

    ofstream out;
    out.open("user-defined-gas-model-test.result");
    if ( out.fail() ) {
	cout << "Error opening file user-defined-gas-model-test.result\n";
	cout << "Bailing out!\n";
	exit(1);
    }

    // 1. Test a minimal implementation by the user.
    //    This relies on some of the default implementations
    //    in the base class.
    UD_gas_model gmodel("ideal-air-minimal-udm.lua");
    Gas_data Q(&gmodel);

    out << "------------------------------------------------" << endl;
    out << "Title   'user-defined-gas-model'" << endl;
    out << "File    'user-defined-gas-model.cxx'" << endl;
    out << "------------------------------------------------" << endl;


    // thermo properties given rho, e
    Q.rho = rho;
    Q.massf[0] = 1.0;
    Q.e[0] = e;
    gmodel.eval_thermo_state_rhoe(Q);
    bool result = ( fabs(Q.p - p) <= tol && 
		    fabs(Q.T[0] - T) <= tol &&
		    fabs(Q.a - a) <= tol );
    
    out << "Test    'minimal-ideal-air: eval_thermo_properties_rhoe'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;
    
    // thermo propertie given p, T
    Q.p = p;
    Q.T[0] = T;
    gmodel.eval_thermo_state_pT(Q);
    result = ( fabs(Q.rho - rho) <= tol &&
	       fabs(Q.e[0] - e) <= tol &&
	       fabs(Q.a - a) <= tol );
    out << "Test    'minimal-ideal-air: eval_thermo_properties_pT'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // thermo properties given rho, T
    Q.rho = rho;
    Q.T[0] = T;
    gmodel.eval_thermo_state_rhoT(Q);
    result = ( fabs(Q.p - p) <= tol && 
	       fabs(Q.e[0] - e) <= tol &&
	       fabs(Q.a - a) <= tol );
    out << "Test    'minimal-ideal-air: eval_thermo_properties_rhoT'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // thermo properties given rho, p
    Q.rho = rho;
    Q.p = p;
    gmodel.eval_thermo_state_rhop(Q);
    result = ( fabs(Q.T[0] - T) <= tol && 
	       fabs(Q.e[0] - e) <= tol && 
	       fabs(Q.a - a) <= tol);
    out << "Test    'minimal-ideal-air: eval_thermo_properties_rhop'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;


    // derivative dTdp at constant rho
    Q.p = p;
    Q.rho = rho;
    double dTdp = 1.0/(rho*R);
    int status;
    double deriv = gmodel.dTdp_const_rho(Q, status);
    result = (fabs(dTdp - deriv) <= tol);
    out << "Test    'minimal-ideal-air: dTdp_const_rho'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dTdrho at constant p
    Q.T[0] = T;
    Q.rho = rho;
    double dTdrho = (-1.0*p)/(R*rho*rho);
    deriv = gmodel.dTdrho_const_p(Q, status);
    result = (fabs(dTdrho - deriv) <= tol);
    out << "Test    'minimal-ideal-air: dTdrho_const_p'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dpdrho at constant T
    Q.rho = rho;
    Q.T[0] = T;
    double dpdrho = R*Q.T[0];
    deriv = gmodel.dpdrho_const_T(Q, status);
    // Use a different tolerance on this test.
    double tol2 = 1.0e-6;
    result = (fabs(dpdrho - deriv) <= tol2);
    out << "Test    'minimal-ideal-air: dpdrho_const_T'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dedT at constant v
    Q.T[0] = T;
    Q.rho = rho;
    double gamma = 1.4;
    double Cv = R / (gamma - 1.0);
    deriv = gmodel.dedT_const_v(Q, status);
    result = (fabs(Cv - deriv) <= tol);
    out << "Test    'minimal-ideal-air: dedT_const_v'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dhdT at constant p
    Q.T[0] = T;
    Q.p = p;
    double Cp = R + Cv;
    deriv = gmodel.dhdT_const_p(Q, status);
    result = (fabs(Cp - deriv) <= tol);
    out << "Test    'minimal-ideal-air: dhdT_const_p'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // gas constant
    double gc = gmodel.R(Q, status);
    result = (fabs(R - gc) <= tol);
    out << "Test    'minimal-ideal-air: R'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // tranport coefficients calculation
    Q.T[0] = T;
    // Sutherland values for air from White (2006), Table 1.2
    double mu0 = 1.716e-5; 
    double T0 = 273.0;
    double S = 111.0;
    double mu_expected = mu0 * pow(T/T0, 3.0/2.0) * (T0 + S) / (T + S);
    // Sutherland values for air from White (2006), Table 1.3
    double k0 = 0.0241;
    T0 = 273.0;
    S = 194.0;
    double k_expected = k0 * pow(T/T0, 3.0/2.0) * (T0 + S) / (T + S);

    gmodel.eval_transport_coefficients(Q);
    
    result = (fabs(mu_expected - Q.mu) <= tol &&
	      fabs(k_expected - Q.k[0]) <= tol);
    out << "Test    'minimal-ideal-air: eval_tranport_coefficients'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // 2. Test a maximal implementation by the user.
    //    This assumes that the user has provided a function for
    //    all available user-controlled methods.
    UD_gas_model gmodel2("ideal-air-maximal-udm.lua");
    // thermo properties given rho, e
    Q.rho = rho;
    Q.e[0] = e;
    gmodel2.eval_thermo_state_rhoe(Q);
    result = ( fabs(Q.p - p) <= tol && 
	       fabs(Q.T[0] - T) <= tol &&
	       fabs(Q.a - a) <= tol );
    
    out << "Test    'maximal-ideal-air: eval_thermo_properties_rhoe'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;
    
    // thermo propertie given p, T
    Q.p = p;
    Q.T[0] = T;
    gmodel2.eval_thermo_state_pT(Q);
    result = ( fabs(Q.rho - rho) <= tol &&
	       fabs(Q.e[0] - e) <= tol &&
	       fabs(Q.a - a) <= tol );
    out << "Test    'maximal-ideal-air: eval_thermo_properties_pT'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // thermo properties given rho, T
    Q.rho = rho;
    Q.T[0] = T;
    gmodel2.eval_thermo_state_rhoT(Q);
    result = ( fabs(Q.p - p) <= tol && 
	       fabs(Q.e[0] - e) <= tol &&
	       fabs(Q.a - a) <= tol );
    out << "Test    'maximal-ideal-air: eval_thermo_properties_rhoT'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // thermo properties given rho, p
    Q.rho = rho;
    Q.p = p;
    gmodel2.eval_thermo_state_rhop(Q);
    result = ( fabs(Q.T[0] - T) <= tol && 
	       fabs(Q.e[0] - e) <= tol && 
	       fabs(Q.a - a) <= tol);
    out << "Test    'maximal-ideal-air: eval_thermo_properties_rhop'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;


    // derivative dTdp at constant rho
    Q.p = p;
    Q.rho = rho;
    deriv = gmodel2.dTdp_const_rho(Q, status);
    result = (fabs(dTdp - deriv) <= tol);
    out << "Test    'maximal-ideal-air: dTdp_const_rho'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dTdrho at constant p
    Q.T[0] = T;
    Q.rho = rho;
    deriv = gmodel2.dTdrho_const_p(Q, status);
    result = (fabs(dTdrho - deriv) <= tol);
    out << "Test    'maximal-ideal-air: dTdrho_const_p'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dpdrho at constant T
    Q.rho = rho;
    Q.T[0] = T;
    deriv = gmodel2.dpdrho_const_T(Q, status);
    result = (fabs(dpdrho - deriv) <= tol);
    out << "Test    'maximal-ideal-air: dpdrho_const_T'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dedT at constant v
    Q.T[0] = T;
    Q.rho = rho;
    deriv = gmodel2.dedT_const_v(Q, status);
    result = (fabs(Cv - deriv) <= tol);
    out << "Test    'maximal-ideal-air: dedT_const_v'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dhdT at constant p
    Q.T[0] = T;
    Q.p = p;
    deriv = gmodel2.dhdT_const_p(Q, status);
    result = (fabs(Cp - deriv) <= tol);
    out << "Test    'maximal-ideal-air: dhdT_const_p'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // gas constant
    gc = gmodel2.R(Q, status);
    result = (fabs(R - gc) <= tol);
    out << "Test    'maximal-ideal-air: R'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;
    // tranport coefficients calculation
    Q.T[0] = T;

    gmodel2.eval_transport_coefficients(Q);
    
    result = (fabs(mu_expected - Q.mu) <= tol &&
	      fabs(k_expected - Q.k[0]) <= tol);
    out << "Test    'maximal-ideal-air: eval_transport_coefficients'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'UD_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;


    out.close();
    
    return 0;
}
