// globalconfig_demo.d

import std.stdio;
import globalconfig;

void main()
{
    writeln("test globalconfig module");
    GlobalConfig.nblocks = 21;
    GlobalConfig.times ~= 1.0;
    GlobalConfig.times ~= 2.0;
    writeln("dimensions= ", GlobalConfig.dimensions);
    writeln("nblocks= ", GlobalConfig.nblocks);
    writeln("times(array)= ", GlobalConfig.times);
    writeln("done.");
}

/+ Transcript:
peterj@helmholtz ~/cfcfd3/dlang/eilmer4/play $ dmd globalconfig_demo.d globalconfig.d
peterj@helmholtz ~/cfcfd3/dlang/eilmer4/play $ ./globalconfig_demo 
test globalconfig module
dimensions= 2
nblocks= 21
times(array)= [1, 2]
done.
+/
