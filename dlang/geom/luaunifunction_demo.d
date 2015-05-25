/**
 * luaunifunction_demo.d Demonstrate some of the behaviour of the UnivariateFunctions.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-02-26
 */

import std.stdio;
import util.lua;
import univariatefunctions;
import luaunifunction;

void main()
{
    writeln("Begin demonstration of LuaD connection to UnivariateFunctions.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerUnivariateFunctions(L);
    luaL_dostring(L, `
-- Make a function and sample it.
myf = LinearFunction:new{t0=0.0, t1=3.0}
print("Try evaluating a point midway on the line.")
pt = myf(0.5)
print("pt=", pt)
print("Or with an eval.")
pt1 = myf:eval(0.5)
print("pt1=", pt1)
myf2 = myf:copy()
pt2 = myf2(0.5)
print("pt2=", pt2)
myrf = RobertsFunction:new{end0=true, end1=false, beta=1.01}
print("myrf(0.5)=", myrf(0.5))
    `);
    writeln("Done with demo.");
}
