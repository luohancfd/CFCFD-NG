// Author: Rowan J. Gollan
// Date: 15-May-2009
// Place: NASA Langley, Hampton, Virginia, USA
//

#include <cstdlib>
#include <iostream>
#include <string>

#include "../../../lib/util/source/useful.h"
#include "thermo-interpolator.hh"
#include "rhoe-interpolator.hh"
#include "rhop-interpolator.hh"
#include "rhoT-interpolator.hh"
#include "pT-interpolator.hh"

using namespace std;

Thermo_interpolator::
Thermo_interpolator() {}

Thermo_interpolator::
~Thermo_interpolator() {}

Thermo_interpolator* create_Thermo_interpolator(string name)
{
    if ( name == "rhoe" ) {
	return new Rhoe_interpolator();
    }
    else if ( name == "rhop" ) {
	return new Rhop_interpolator();
    }
    else if ( name == "rhoT" ) {
	return new RhoT_interpolator();
    }
    else if ( name == "pT" ) {
	return new PT_interpolator();
    }
    else {
	cout << "create_Thermo_interpolator():\n";
	cout << "The chosen interpolator: " << name << " is unknown.\n";
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
}
