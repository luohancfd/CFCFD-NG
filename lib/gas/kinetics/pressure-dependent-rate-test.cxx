// Author: Rowan J. Gollan, Brendan O'Flaherty
// Date: 17-Apr-2009
// Place: NIA, Hampton, Virginia, USA
//        UQ, St Lucia, Queensland

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "../model/gas-model.hh"
#include "pressure-dependent-rate.hh"

using namespace std;

const char *inp = "pressure-dependent-rate-test.inp";

int main()
{
    const double tol = 1.0e-9;
    // Reference calculation for test.
    // const double T = ....;
    // const double rho = ...;
    // vector<double> mf; mf.resize(4);
    // mf[0] = ...; mf[1] = ...; mf[2] = ...; mf[3] = ...;
    // const double k = ....;

    // Initialise gas model
    Gas_model *g = create_gas_model("methane-mixture.lua");

    // Initialise presure depdendent rate coefficient model
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    if ( luaL_dofile(L, inp) != 0 ) {
	ostringstream ost;
	ost << "Error in input file for testing pressure dependent reaction\n";
	ost << "rate coefficient model: " << inp << endl;
    }

    lua_getglobal(L, "test_rate");
    Pressure_dependent pd(L, *g);
    lua_close(L);

    ofstream out;
    out.open("pressure-dependent-rate-test.result");
    if ( out.fail() ) {
	cout << "Error opening file pressure-dependent-rate-test.result\n";
	cout << "Bailing out!\n";
	exit(1);
    }

    out << "-------------------------------------------------" << endl;
    out << "Title   'generalised-Arrhenius'" << endl;
    out << "File    'generalised-Arrhenius.cxx'" << endl;
    out << "-------------------------------------------------" << endl;
    
    Gas_data Q(g);
    Q.T[0] = T;
    Q.rho = rho;
    // Set mass fractions...
    for ( size_t isp = 0; isp < mf.size(); ++isp ) Q.massf[isp] = mf[isp];

    // Evaluate reaction rate coefficient...
    pd.eval(Q);

    bool result = ( fabs(pd.k() - k) <= tol );

    out << "Test    'eval'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Pressure_dependent'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;
    
    out.close();

    return 0;
}
