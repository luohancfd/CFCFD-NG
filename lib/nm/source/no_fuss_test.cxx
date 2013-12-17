/** \file no_fuss_test.cxx
 * \brief Testing of the no_fuss_lineare_algebra module.
 *
 * \author Rowan J Gollan
 * \date 23-Apr-2006
 *
 **/

#include <iostream>
#include <vector>

#include "no_fuss_linear_algebra.hh"

using namespace std;

int main()
{
    cout << "Begin no_fuss_test.x...\n";

    cout << "Initialize matrices using various constructors.\n";
    cout << "-----------------------------------------------\n\n";

    cout << "1. Valmatrix z( 4, 3 );\n";
    Valmatrix z( 4, 3 );
    cout << z;

    cout << "Test some matrix operations.\n";
    cout << "----------------------------\n\n";

    cout << "1. Solve the 2x2 system Ax = b where ... \n";
    Valmatrix A(2,2);
    A.set(0, 0, 0.00300); A.set(0, 1, 59.14);
    A.set(1, 0, 5.291); A.set(1, 1, -6.130);
    vector<double> B(2);
    B[0] = 59.17;
    B[1] = 46.78;

    cout << "A= \n";
    cout << A;
    cout << "B= ";
    print_vector(B);

    vector<double> x(2);

    gaussian_elimination(A, x, B);

    cout << "x= ";
    print_vector(x);
    cout << "The answer should be x= [10.0, 1.0]\n";
    cout << "See Example 2. on p. 340 of Burden & Faires\n";


    cout << "\n2. Solve the matrix system Ax = b where...\n";
    Valmatrix A2( 3, 3);
    A2.set(0,0,1.0);  A2.set(0,1,-1.0);  A2.set(0,2,0.0);
    A2.set(1,0,-2.0);  A2.set(1,1,2.0);  A2.set(1,2,-1.0);
    A2.set(2,0,0.0);  A2.set(2,1,1.0);  A2.set(2,2,-2.0);

    cout << "A= \n";
    cout << A2;
    vector<double> B2(3);
    B2[0] = 2.0; B2[1] = -1.0; B2[2] = 6.0;
    cout << "B= ";
    print_vector(B2);

    vector<double> x2(3);
    gaussian_elimination( A2, x2, B2 );
    cout << "After call to gaussian_elimination():\n";
    cout << "A= \n";
    cout << A2;
    cout << "B= ";
    print_vector(B2);
    cout << "x= ";
    print_vector(x2);

    cout << "The answer should be x= [2.0, 0.0, -3.0]\n";


    cout << "3. Repeat the same exercise using virtual scaling.\n";
    
    A2.set(0,0,1.0);  A2.set(0,1,-1.0);  A2.set(0,2,0.0);
    A2.set(1,0,-2.0);  A2.set(1,1,2.0);  A2.set(1,2,-1.0);
    A2.set(2,0,0.0);  A2.set(2,1,1.0);  A2.set(2,2,-2.0);

    B2[0] = 2.0; B2[1] = -1.0; B2[2] = 6.0;

    gaussian_elimination( A2, x2, B2, true );
   
    cout << "x= ";
    print_vector(x2);

    cout << "4. Example from p. 130 of Gerald and Wheatley...\n";
    cout << "Solve Ax = b where:\n";
    
    Valmatrix AA( 4, 4);
    AA.set(0,0,0.0); AA.set(0,1,2.0); AA.set(0,2,0.0); AA.set(0,3,1.0);
    AA.set(1,0,2.0); AA.set(1,1,2.0); AA.set(1,2,3.0); AA.set(1,3,2.0);
    AA.set(2,0,4.0); AA.set(2,1,-3.0); AA.set(2,2,0.0); AA.set(2,3,1.0);
    AA.set(3,0,6.0); AA.set(3,1,1.0); AA.set(3,2,-6.0); AA.set(3,3,-5.0);
    cout << "A= \n";
    cout << AA;
    
    vector<double> BB(4);
    BB[0] = 0.0; BB[1] = -2.0; BB[2] = -7.0; BB[3] = 6.0;
    cout << "BB= ";
    print_vector(BB);

    vector<double> xx(4);
    gaussian_elimination(AA, xx, BB);
    cout << "After call to gaussian_elimination():\n";
    cout << "A= \n";
    cout << AA;
    cout << "BB= ";
    print_vector(BB);
    cout << "xx= ";
    print_vector(xx);
    cout << "The answer should be x= [ -0.5, 1, 0.333, -2 ]\n";


    Valmatrix T(4,3);
    vector<double> tau1(3);
    T.set(0, 0, 4.0); T.set(0, 1, 3.0); T.set(0, 2, 5.0);
    T.set(1, 0, 19.0); T.set(1, 1, 47.9); T.set(1, 2, 7.0);
    T.set(2, 0, -2.0); T.set(2, 1, 107.0); T.set(2, 2, -89.0);
    T.set(3, 0, 4.5); T.set(3, 1, 9.1); T.set(3, 2, -23.0);
    cout << "\n5. Testing least squarse solver with QR algorithm\n";
    cout << "with matrix T=\n";
    cout << T;
    cout << "and solution vector, b= \n";
    vector<double> p(3);
    vector<double> b(4);
    b[0] = 33.0;
    b[1] = -17.5;
    b[2] = -88.3;
    b[3] = 5.4;

    print_vector(b);

    least_squares_solve(T, p, b);

    cout << "Answer= " << endl;
    print_vector(p);

    cout << "The answer given by octave is:\n";
    cout << "x = \n\n";
    cout << "   2.19390\n";
    cout << "  -1.07473\n";
    cout << "  -0.32910\n";
    
    cout << "\n6. Try the QR factorization to solve a normal, square system.\n";
    cout << "Repeat Test 4.\n";
    
    AA.set(0,0,0.0); AA.set(0,1,2.0); AA.set(0,2,0.0); AA.set(0,3,1.0);
    AA.set(1,0,2.0); AA.set(1,1,2.0); AA.set(1,2,3.0); AA.set(1,3,2.0);
    AA.set(2,0,4.0); AA.set(2,1,-3.0); AA.set(2,2,0.0); AA.set(2,3,1.0);
    AA.set(3,0,6.0); AA.set(3,1,1.0); AA.set(3,2,-6.0); AA.set(3,3,-5.0);
    cout << "A= \n";
    cout << AA;
    
    BB[0] = 0.0; BB[1] = -2.0; BB[2] = -7.0; BB[3] = 6.0;
    cout << "BB= ";
    print_vector(BB);

    QR_solve(AA, xx, BB);
    cout << "After call to QR_solve():\n";
    cout << "A= \n";
    cout << AA;
    cout << "BB= ";
    print_vector(BB);
    cout << "xx= ";
    print_vector(xx);
    cout << "The answer should be x= [ -0.5, 1, 0.333, -2 ]\n";

    return 0;

}
