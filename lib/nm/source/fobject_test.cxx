/** \file fobject_test.cxx
 *  \ingroup nm
 *  \brief Implementation of the function-object classes. 
 *  \author PJ
 *  \version 11-Jan-2006 initial coding
 *  \version 19-Oct-2008 Added LinearFunction2 and BilinearFunction.
 *
 * For use in the libgeom2 classes so that we can pass around
 * simple function objects rather than functions themselves.
 *
 */
#include "iostream"
#include "stdio.h"
#include "fobject.hh"
using namespace std;

int main() 
{
    cout << "Begin fobject_test..." << endl;
    LinearFunction f1 = LinearFunction();
    cout << "f1: " << f1 << endl;
    LinearFunction f2 = LinearFunction(-1.0, 1.0);
    cout << "f2: " << f2 << endl;
    RobertsClusterFunction f3 = RobertsClusterFunction(0, 1, 1.1);
    cout << "f3: " << f3 << endl;
    LinearFunction2 f4 = LinearFunction2(0.0, 1.0);
    cout << "f4: " << f4 << endl;
    int n = 10;
    double dt = 1.0/n;
    cout << "         t      f1(t)      f2(t)      f3(t)      f4(t)" << endl;
    for ( int i = 0; i <= n; ++i ) {
	double t = i * dt;
	printf("%10.4f %10.4f %10.4f %10.4f %10.4f\n", t, 
	       f1.eval(t), f2.eval(t), f3.eval(t), f4.eval(t) );
    }

    cout << "Sampled points for f3..." << endl;
    vector<double>* tvalues = f3.distribute_parameter_values(n+1, 0.0, 1.0);
    cout << "tvalues.size()=" << tvalues->size() << endl;
    for ( size_t i = 0; i < tvalues->size(); ++i ) {
	cout << (*tvalues)[i] << " ";
    }
    cout << endl;
    delete tvalues;

    cout << "Sampled points for reversed clustering of f3..." << endl;
    f3.reverse_clustering();
    tvalues = f3.distribute_parameter_values(n+1, 0.0, 1.0);
    cout << "tvalues.size()=" << tvalues->size() << endl;
    for ( size_t i = 0; i < tvalues->size(); ++i ) {
	cout << (*tvalues)[i] << " ";
    }
    cout << endl;
    delete tvalues;
    
    cout << "Making valliammai instance f5" << endl;
    vector<double> *tvalue;
    ValliammaiFunction f5 = ValliammaiFunction(1e-5, 1e-5, 0.01,100);
    tvalue = f5.distribute_parameter_values(100, 0.0, 1.0);

    cout << "Sampled points for f5 (ValliammaiFunction)..." << endl;
    for ( size_t i = 0; i < tvalue->size(); ++i ) {
	cout << (*tvalue)[i] << " ";
    }
    cout << endl;

    cout << "Try out BilinearFunction" << endl;
    BilinearFunction f6 = BilinearFunction(0.0, 1.0, 2.0, 3.0);
    cout << "f6: " << f6 << endl;
    int nr = 3;
    double dr = 1.0/nr;
    int ns = 4;
    double ds = 1.0/ns;
    cout << "         r          s      f6(t)" << endl;
    for ( int i = 0; i <= nr; ++i ) {
	double r = i * dr;
	for ( int j = 0; j <= ns; ++j ) {
	    double s = j * ds;
	    printf("%10.4f %10.4f %10.4f\n", r, s, f6.eval(r,s) );
	}
    }

    cout << "Done." << endl;
    delete tvalue;
    return 0;
}

