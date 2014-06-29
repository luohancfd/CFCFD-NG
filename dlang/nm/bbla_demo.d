// bbla_demo.d

import std.stdio;
import bbla;

void main()
{
    writeln("Bare-bones linear algebra demo...");
    Matrix a = eye(3);
    writeln("a=", a);
    Matrix b = zeros(2,3);
    Matrix c = transpose(b);
    b[1,2] = 99.0;
    c[0,0] = 1.0; c[1,1] = 1.0;
    writeln("b=", b);
    writeln("c=", c);

    Matrix e = new Matrix([1.0, 2.0, 3.0]);
    writeln("e= ", e);
    Matrix e2 = new Matrix([1, 2, 3], "row");
    writeln("e2= ", e2);
    Matrix f = new Matrix([[1.0,2.0,3.0],[4.0,5.0,6.0]]);
    writeln("f= ", f);
    writeln("2*f= ", 2*f);
    writeln("f+2= ", f+2);

    Matrix g = dot(f,c);
    writeln("g= ", g);

    Matrix Ab = new Matrix([[0.0,  2.0,  0.0,  1.0,  0.0],
			    [2.0,  2.0,  3.0,  2.0, -2.0],
			    [4.0, -3.0,  0.0,  1.0, -7.0],
			    [6.0,  1.0, -6.0, -5.0,  6.0]]);
    Matrix Aonly = Ab.sliceDup(0, 4, 0, 4);
    Matrix bonly = Ab.sliceDup(0, 4, 4, 5);
    gaussJordanElimination(Ab);
    writeln("Ab= ", Ab);
    double[] x = Ab.getColumn(4);
    writeln("x= ", x);
    Matrix rhs = dot(Aonly, new Matrix(x));
    writeln("rhs= ", rhs);
    Matrix residual = rhs - bonly;
    writeln("residual= ", residual);
    writeln("Done.");
}
