// Author: Rowan J. Gollan
// Version: 25-May-2008
//            Initial coding.
//

#include <iostream>
#include <iomanip>

#include "gas_data.hh"
#include "perfect-gas-EOS.hh"

using namespace std;

int main()
{
    const double tol = 1.0e-9;
    double rho = 1.2;
    double T = 300.0;
    double p = 103356.0;
    double R = 287.1;
    // ---- Test object creation and methods
    Perfect_gas pg;
    Gas_data Q(pg);
    
    cout << setprecision(12) << showpoint;

    cout << "title 'perfect-gas-EOS-test'" << endl;
    cout << "file  'perfect-gas-EOS.cxx'" << endl;

    // ---- Test the functions

    double p_expected = p;
    double p_result = pg_pressure(rho, T, R);
    cout << "" << endl;
    cout << "unit { name='pg_pressure',\n";
    cout << "       type='function',\n";
    cout << "       expected=" << p_expected << ",\n";
    cout << "       result=" << p_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    double T_expected = T;
    double T_result = pg_temperature(rho, p, R);
    cout << "" << endl;
    cout << "unit { name='pg_temperature',\n";
    cout << "       type='function',\n";
    cout << "       expected=" << T_expected << ",\n";
    cout << "       result=" << T_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    double rho_expected = rho;
    double rho_result = pg_density(T, p, R);
    cout << "" << endl;
    cout << "unit { name='pg_density',\n";
    cout << "       type='function',\n";
    cout << "       expected=" << rho_expected << ",\n";
    cout << "       result=" << rho_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    Q.T = T; Q.rho = rho; Q.p = p; Q.R = R;

    p_expected = p;
    p_result = pg.eval_pressure(Q);
    cout << "" << endl;
    cout << "unit { name='eval_pressure',\n";
    cout << "       type='method',\n";
    cout << "       class='Perfect_gas',\n";
    cout << "       expected=" << p_expected << ",\n";
    cout << "       result=" << p_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    T_expected = T;
    T_result = pg.eval_temperature(Q);
    cout << "" << endl;
    cout << "unit { name='eval_temperature',\n";
    cout << "       type='method',\n";
    cout << "       class='Perfect_gas',\n";
    cout << "       expected=" << T_expected << ",\n";
    cout << "       result=" << T_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    rho_expected = rho;
    rho_result = pg.eval_density(Q);
    cout << "" << endl;
    cout << "unit { name='eval_density',\n";
    cout << "       type='method',\n";
    cout << "       class='Perfect_gas',\n";
    cout << "       expected=" << rho_expected << ",\n";
    cout << "       result=" << rho_result << ",\n";
    cout << "       test=equals(" << tol << ") }\n";

    return 0;
}
