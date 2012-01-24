// Author: Rowan J. Gollan
// Version: 10-Oct-2008


#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "generalised-Arrhenius.hh"

using namespace std;

const char *inp = "generalised-Arrhenius-test.inp";

int main()
{
    const double tol = 1.0e-9;
    const double T = 1000.0;
    const double k = 7.640388414e-39;

    // Initialise model
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    if ( luaL_dofile(L, inp) != 0 ) {
	ostringstream ost;
	ost << "Error in input file for testing generalised Arrhenius reaction\n";
	ost << "rate coefficient model: " << inp << endl;
    }

    Gas_model *g = create_gas_model("ideal-air.lua"); // dummy gas model

    lua_getglobal(L, "test_rate");
    Generalised_Arrhenius ga(L, *g);
    lua_close(L);

    ofstream out;
    out.open("generalised-Arrhenius-test.result");
    if ( out.fail() ) {
	cout << "Error opening file generalised-Arrhenius-test.result\n";
	cout << "Bailing out!\n";
	exit(1);
    }

    out << "-------------------------------------------------" << endl;
    out << "Title   'generalised-Arrhenius'" << endl;
    out << "File    'generalised-Arrhenius.cxx'" << endl;
    out << "-------------------------------------------------" << endl;
    
    Gas_data Q(g);
    Q.T[0] = T;

    ga.eval(Q);

    bool result = ( fabs(ga.k() - k) <= tol );

    out << "Test    'eval'" << endl;
    out << "Type    'method'" << endl;
    out << "Class   'Generalised_Arrhenius'" << endl;
    out << "Result  " << (result ? "'passed'" : "'failed'") << endl;
    out << endl;
    
    out.close();

    return 0;
}
