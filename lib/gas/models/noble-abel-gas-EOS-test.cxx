// Author: Brendan T. O'Flaherty
//  based on perfect-gas-EOS-test code by Rowan J. Gollan
// Version: 14-May-2009
//            Initial coding.
//

#include <iostream>
#include <iomanip>

#include "gas_data.hh"
#include "noble-abel-gas-EOS.hh"

using namespace std;

int main()
{
    const double tol = 1.0e-9;
    double rho = 1.2;
    double T = 300.0;
    double p = 103360.6155;
    double R = 287.1;
    double nu_0 = 3.7212e-5;
    // ---- Test object creation and methods
    Noble_Abel_gas nag;
    Gas_data Q(nag);
    
    cout << setprecision(12) << showpoint;

    cout << "title 'noble-abel-gas-EOS-test'" << endl;
    cout << "file  'noble-abel-gas-EOS.cxx'" << endl;

    // ---- Test the functions

    double p_expected = p;
    double p_result = nag_pressure(rho, T, R, nu_0);
    cout << "" << endl;
    cout << "unit { name='nag_pressure',\n";
    cout << "       type='function',\n";
    cout << "       expected=" << p_expected << ",\n";
    cout << "       result=" << p_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    double T_expected = T;
    double T_result = nag_temperature(rho, p, R, nu_0);
    cout << "" << endl;
    cout << "unit { name='nag_temperature',\n";
    cout << "       type='function',\n";
    cout << "       expected=" << T_expected << ",\n";
    cout << "       result=" << T_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    double rho_expected = rho;
    double rho_result = nag_density(T, p, R, nu_0);
    cout << "" << endl;
    cout << "unit { name='nag_density',\n";
    cout << "       type='function',\n";
    cout << "       expected=" << rho_expected << ",\n";
    cout << "       result=" << rho_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    Q.T = T; Q.rho = rho; Q.p = p; Q.R = R;

    p_expected = p;
    p_result = nag.eval_pressure(Q);
    cout << "" << endl;
    cout << "unit { name='eval_pressure',\n";
    cout << "       type='method',\n";
    cout << "       class='Noble_Abel_gas',\n";
    cout << "       expected=" << p_expected << ",\n";
    cout << "       result=" << p_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    T_expected = T;
    T_result = nag.eval_temperature(Q);
    cout << "" << endl;
    cout << "unit { name='eval_temperature',\n";
    cout << "       type='method',\n";
    cout << "       class='Noble_Abel_gas',\n";
    cout << "       expected=" << T_expected << ",\n";
    cout << "       result=" << T_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    rho_expected = rho;
    rho_result = nag.eval_density(Q);
    cout << "" << endl;
    cout << "unit { name='eval_density',\n";
    cout << "       type='method',\n";
    cout << "       class='Noble_Abel_gas',\n";
    cout << "       expected=" << rho_expected << ",\n";
    cout << "       result=" << rho_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    return 0;
}
