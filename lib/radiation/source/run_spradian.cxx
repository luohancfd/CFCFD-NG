/** \file run_spradian.cxx
 *  \ingroup radiation2
 *  \brief Run the spradian 'radipac' routine for a given input file
 *  \author Daniel F Potter        
 *  \version 16-Apr-2012
 *
 */
#include <cstdlib>
#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <time.h>
#include <sstream>

#include "../../gas/models/gas-model.hh"
#include "../../util/source/useful.h"

#include "spectral_model.hh"
#include "spradian.hh"
#include "LOS_pieces.hh"

using namespace std;

int main( int argc, char *argv[] )
{
    if ( argc!=2 ) {
        cout << "run_spradian::main()" << endl
             << "Usage: run_spradian.x <param.fname>" << endl;
        exit( BAD_INPUT_ERROR );
    }
    
    SpradianParams params(argv[1]);
    params.write_to_file("spradian.params.check");

    params.call_radipac();
    
    return 0;
}
