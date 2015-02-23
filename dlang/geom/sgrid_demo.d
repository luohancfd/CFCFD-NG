// sgrid_demo.d

import std.stdio;
import std.conv;
import geom;
import gpath;
import surface;
import univariatefunctions;
import sgrid;

void main()
{
    writeln("Begin structured-grid demo...");
    auto p00 = Vector3(0.0, 0.1);
    auto p10 = Vector3(1.0, 0.4);
    auto p11 = Vector3(1.0, 1.1);
    auto p01 = Vector3(0.0, 1.1);
    auto my_patch = new AOPatch(p00, p10, p11, p01);
    auto cf = [new LinearFunction(), new LinearFunction(), 
	       new LinearFunction(), new LinearFunction()];
    auto my_grid = new StructuredGrid(my_patch, 11, 21, cf);
    writeln("grid point 5 5 at x=", my_grid[5,5].x, " y=", my_grid[5,5].y);
    my_grid.write_to_text_file("test_grid.vtk", true);

    writeln("Import GridPro grid...");
    auto gpgrid = import_gridpro_grid("../../examples/eilmer3/3D/gridpro-import/blk.tmp");
    foreach (i; 0 .. gpgrid.length) {
	gpgrid[i].write_to_text_file("gpgrid-"~to!string(i)~".vtk", true);
    }
    writeln("Done.");
}
