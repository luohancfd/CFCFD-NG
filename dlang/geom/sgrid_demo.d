// sgrid_demo.d

import std.stdio;
import geom;
import gpath;
import surface;
import univariatefunctions;
import sgrid;

void main()
{
    writeln("Begin structured-grid demo...");
    auto p00 = Vector3(0.0, 0.1);
    auto p10 = Vector3(1.0, 0.1);
    auto p11 = Vector3(1.0, 1.1);
    auto p01 = Vector3(0.0, 1.1);
    auto my_patch = new CoonsPatch(p00, p10, p11, p01);
    auto cf = [new LinearFunction(), new LinearFunction(), 
	       new LinearFunction(), new LinearFunction()];
    auto my_grid = new StructuredGrid(my_patch, 11, 21, cf);
    writeln("grid point 5 5 at x=", my_grid.grid[5][5][0].x, 
	    " y=", my_grid.grid[5][5][0].y);
    my_grid.write_to_text_file("test_grid.vtk", true);
    writeln("Done.");
}
