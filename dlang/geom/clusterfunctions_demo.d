// clusterfunctions_demo.d

import std.stdio;
import clusterfunctions;

void main()
{
    writeln("Begin demonstration of the clustering functions.");
    auto cf = robertsFunction(false, true, 1.1);
    auto cf2 = linearFunction(1.0, 0.0);
    foreach(i; 0 .. 11) {
	double x = 0.1*i;
	writeln("x= ", x, " roberts_y= ", cf(x), " linear_y= ", cf2(x));
    }
} // end main
