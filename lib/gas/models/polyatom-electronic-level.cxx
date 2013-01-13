// Author: Daniel F. Potter
// Version: 24-Mar-2010
//          Ported from lib/radiation/source/diatomic_radiator.cxx

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "polyatom-electronic-level.hh"
#include "physical_constants.hh"

#include "../../util/source/useful.h"

using namespace std;

Polyatom_electronic_level::Polyatom_electronic_level( vector<double> lev_data )
{
    E = lev_data[0] * PC_c * PC_h_SI;
    r_e = lev_data[1];
    g = int(lev_data[2]);
    D = lev_data[3] * PC_c * PC_h_SI;

    A0 = lev_data[4] * PC_c * PC_h_SI;
    B0 = lev_data[5] * PC_c * PC_h_SI;
    C0 = lev_data[6] * PC_c * PC_h_SI;

    sigma = lev_data[7];
    sigma_r = lev_data[8];

    for ( size_t i=9; i<lev_data.size(); ++i )
         omega_e_vec.push_back( lev_data[i] * PC_c * PC_h_SI );
}
