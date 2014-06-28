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
    c[0,0] = 1.0;
    writeln("b=", b);
    writeln("c=", c);

    Matrix e = new Matrix([1.0, 2.0, 3.0]);
    writeln("e= ", e);
    Matrix f = new Matrix([[1.0,2.0,3.0],[4.0,5.0,6.0]]);
    writeln("f= ", f);

    Matrix g = dot(f,c);
    writeln("g= ", g);

    Matrix Ab = new Matrix([[0.0,  2.0,  0.0,  1.0,  0.0],
			    [2.0,  2.0,  3.0,  2.0, -2.0],
			    [4.0, -3.0,  0.0,  1.0, -7.0],
			    [6.0,  1.0, -6.0, -5.0,  6.0]]);
    gauss_jordan_elimination(Ab);
    writeln("Ab= ", Ab);
    writeln("Done.");
}
