// Author: Rowan J. Gollan
// Date: 16-Jul-2008

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "gas_data.hh"
#include "gas-model.hh"
#include "composite-gas-model.hh"

using namespace std;

int main()
{

    const double tol = 1.0e-9;
    // Reference answer values
    const double rho = 1.2;
    const double T = 300.0;
    const double p = 103342.8946278138;
    const double R = 287.063596188371775;
    const double e = 215297.697141278855;
    const double a = 347.227174050529;

    ofstream out;
    out.open("composite-gas-model-test.result");
    if ( out.fail() ) {
 	cout << "Error opening file composite-gas-model-test.result\n";
 	cout << "Bailing out!\n";
 	exit(1);
    }

    Gas_model *ideal = create_gas_model("ideal-air.lua");
    Gas_data Q(ideal);

    out << "---------------------------------------------------" << endl;
    out << "Title    'composite-gas-model'" << endl;
    out << "File     'composite-gas-model.cxx'" << endl;
    out << "---------------------------------------------------" << endl;

    // thermo properties given rho, e
    Q.rho = rho;
    Q.e[0] = e;
    Q.massf[0] = 1.0;
    ideal->eval_thermo_state_rhoe(Q);
    bool result = ( fabs(Q.p - p) <= tol && 
 		    fabs(Q.T[0] - T) <= tol &&
 		    fabs(Q.a - a) <= tol );
    
    out << "Test    'eval_thermo_properties_rhoe'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Composite_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // thermo properties given p, T
    Q.p = p;
    Q.T[0] = T;
    ideal->eval_thermo_state_pT(Q);
    result = ( fabs(Q.rho - rho) <= tol &&
 	       fabs(Q.e[0] - e) <= tol &&
 	       fabs(Q.a - a) <= tol );

    out << "Test    'eval_thermo_properties_pT'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Composite_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // thermo properties given rho, T
    Q.rho = rho;
    Q.T[0] = T;
    ideal->eval_thermo_state_rhoT(Q);
    result = ( fabs(Q.p - p) <= tol && 
 	       fabs(Q.e[0] - e) <= tol &&
 	       fabs(Q.a - a) <= tol );
    out << "Test    'eval_thermo_properties_rhoT'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Composite_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // thermo properties given rho, p
    Q.rho = rho;
    Q.p = p;
    ideal->eval_thermo_state_rhop(Q);
    result = ( fabs(Q.T[0] - T) <= tol && 
 	       fabs(Q.e[0] - e) <= tol && 
 	       fabs(Q.a - a) <= tol);
    out << "Test    'eval_thermo_properties_rhop'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Composite_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dTdp at constant rho
    double dTdp = 1.0/(rho*R);
    int status;
    double deriv = ideal->dTdp_const_rho(Q, status);
    result = (fabs(dTdp - deriv) <= tol);
    out << "Test    'dTdp_const_rho'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Composite_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dTdrho at constant p
    Q.rho = rho;
    Q.p = p;
    double dTdrho = (-1.0*Q.p)/(R*Q.rho*Q.rho);
    deriv = ideal->dTdrho_const_p(Q, status);
    result = (fabs(dTdrho - deriv) <= tol);
    out << "Test    'dTdrho_const_p'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Composite_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dpdrho at constant T
    Q.rho = rho;
    Q.T[0] = T;
    double dpdrho = R*Q.T[0];
    deriv = ideal->dpdrho_const_T(Q, status);
    result = (fabs(dpdrho - deriv) <= tol);
    out << "Test    'dpdrho_const_T'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Composite_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dedT at constant v
    Q.T[0] = T;
    Q.rho = rho;
    double gamma = 1.4;
    double Cv = R / (gamma - 1.0);
    deriv = ideal->dedT_const_v(Q, status);
    result = (fabs(Cv - deriv) <= tol);
    out << "Test    'dedT_const_v'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Composite_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    // derivative dhdT at constant p
    Q.T[0] = T;
    Q.p = p;
    double Cp = R + Cv;
    deriv = ideal->dhdT_const_p(Q, status);
    result = (fabs(Cp - deriv) <= tol);
    out << "Test    'dhdT_const_p'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Composite_gas_model'" << endl;
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
    
    ideal->eval_transport_coefficients(Q);
    
    result = (fabs(mu_expected - Q.mu) <= tol &&
	       fabs(k_expected - Q.k[0]) <= tol);
    out << "Test    'eval_tranport_coefficients'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Composite_gas_model'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;

    out.close();
    delete ideal;
    
    return 0;
}
