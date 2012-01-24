// Author: Daniel F. Potter
// Date: 13-Apr-2010

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "gas_data.hh"
#include "gas-model.hh"

using namespace std;

int main ( int argc, char *argv[] )
{
    if ( argc!=2 ) {
    	cout << "usage: noneq-test.x GASFILE" << endl;
    	exit(0);
    }
    Gas_model * noneq = create_gas_model(argv[1]);
    Gas_data Q(noneq);

    for ( size_t iT=0; iT<Q.T.size(); ++iT )
    	Q.T[iT] = 5000.0;
    Q.massf[0] = 1.0;
    Q.p = 1.0e5;
    noneq->eval_thermo_state_pT(Q);
    noneq->eval_transport_coefficients(Q);
    Q.print_values();
    delete noneq;
    return 0;
}
