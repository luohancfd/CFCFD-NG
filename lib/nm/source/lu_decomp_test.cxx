// Author: Rowan J. Gollan
// Date: 03-Mar-2010
// Place: Poquoson, Virginia, USA
//

#include <iostream>

#include "no_fuss_linear_algebra.hh"
#include "lu_decomp.hh"

using namespace std;

int main()
{
    cout << "Begin test: LU decomposition.\n";
    cout << "--------------------------------------\n";
    cout << endl;

    cout << "1. Initialise LUdcmp object from a Valmatrix.\n";
    Valmatrix A(2,2);
    A.set(0, 0, 0.00300); A.set(0, 1, 59.14);
    A.set(1, 0, 5.291); A.set(1, 1, -6.130);
    cout << "A = \n";
    cout << "| " << A.get(0,0) << "   " << A.get(0,1) << " |\n";
    cout << "| " << A.get(1,0) << "   " << A.get(1,1) << " |\n";
    cout << endl;

    LUdcmp Alu(A);
    cout << "Object initialised.\n";

    cout << "2. Solve the 2x2 system Ax = b where...\n";
    vector<double> b(2, 0.0);
    b[0] = 59.17;
    b[1] = 46.78;
    cout << "b = [ " << b[0] << ", " << b[1] <<  " ] \n";
    
    vector<double> x(2, 0.0);

    Alu.solve(b, x);

    cout << "After solving...\n";
    cout << "x = [ " << x[0] << ", " << x[1] << " ] \n";
    cout << "The answer should be x= [10.0, 1.0]\n";
    cout << "See Example 2. on p. 340 of Burden & Faires\n";
    cout << endl;

    cout << "3. Solve the matrix system Ax = b where...\n";
    Valmatrix A2( 3, 3);
    A2.set(0,0,1.0);  A2.set(0,1,-1.0);  A2.set(0,2,0.0);
    A2.set(1,0,-2.0);  A2.set(1,1,2.0);  A2.set(1,2,-1.0);
    A2.set(2,0,0.0);  A2.set(2,1,1.0);  A2.set(2,2,-2.0);

    LUdcmp A2lu(A2);

    cout << "A= \n";
    cout << A2;
    vector<double> B2(3);
    B2[0] = 2.0; B2[1] = -1.0; B2[2] = 6.0;
    cout << "b= [ " << B2[0] << ", " << B2[1] << ", " << B2[2] << " ]\n";

    vector<double> x2(3);
    A2lu.solve(B2, x2);

    cout << "After call to solve():\n";
    cout << "x= [ " << x2[0] << ", " << x2[1] << ", " << x2[2] << " ]\n";

    cout << "The answer should be x= [2.0, 0.0, -3.0]\n";
    cout << endl;
    
    cout << "4. Example from p.130 of Gerald and Wheatley...\n";
    cout << "   Solve Ax = b where...\n";
    Valmatrix AA(4,4);
    AA.set(0,0,0.0); AA.set(0,1,2.0); AA.set(0,2,0.0); AA.set(0,3,1.0);
    AA.set(1,0,2.0); AA.set(1,1,2.0); AA.set(1,2,3.0); AA.set(1,3,2.0);
    AA.set(2,0,4.0); AA.set(2,1,-3.0); AA.set(2,2,0.0); AA.set(2,3,1.0);
    AA.set(3,0,6.0); AA.set(3,1,1.0); AA.set(3,2,-6.0); AA.set(3,3,-5.0);
    cout << "A= \n";
    cout << AA;

    LUdcmp AAlu(AA);

    vector<double> BB(4);
    BB[0] = 0.0; BB[1] = -2.0; BB[2] = -7.0; BB[3] = 6.0;
    cout << "BB= [ " << BB[0] << ", " << BB[1] << ", " << BB[2] << ", " << BB[3] << " ]\n";

    vector<double> xx(4);

    AAlu.solve(BB, xx);

    cout << "After call to solve():\n";
    cout << "xx= [ " << xx[0] << ", " << xx[1] << ", " << xx[2] << ", " << xx[3] << " ]\n";

    cout << "The answer should be x= [ -0.5, 1, 0.333, -2 ]\n";
    cout << endl;


    cout << "Test done.\n";

    return 0;
}
