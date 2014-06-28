// globaldata_demo.d

import std.stdio;
import globaldata;

void main()
{
    writeln("test globaldata module");
    GlobalData.nblocks = 21;
    GlobalData.times ~= 1.0;
    GlobalData.times ~= 2.0;
    writeln("dimensions= ", GlobalData.dimensions);
    writeln("nblocks= ", GlobalData.nblocks);
    writeln("times(array)= ", GlobalData.times);
    writeln("done.");
}

/+ Transcript:
peterj@helmholtz ~/cfcfd3/dlang/eilmer4/play $ dmd globaldata_demo.d globaldata.d
peterj@helmholtz ~/cfcfd3/dlang/eilmer4/play $ ./globaldata_demo 
test globaldata module
dimensions= 2
nblocks= 21
times(array)= [1, 2]
done.
+/
