/** e4prep.d
 * Eilmer4, perpare grids and initial flow state.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-02-24
 */

import std.stdio;
import std.getopt;
import luad.all;
import luageom;
import luaflowstate;
import geom;
import gas;


void main(string[] args)
{
    writeln("Eilmer4, prepare grids and initial flow state.");
    writeln("Revision: PUT_REVISION_STRING_HERE");

    string msg = "Usage:                            Comment:\n";
    msg       ~= "e4prep [--job=<string>]           file names built from this string\n";
    msg       ~= "       [--verbosity=<int>]        defaults to 0\n";
    msg       ~= "       [--no-flow]                write grids only\n";
    msg       ~= "       [--help]                   writes this message\n";
    if ( args.length < 2 ) {
	writeln("Too few arguments.");
	write(msg);
	exit(1);
    }
    string jobName = "";
    int verbosityLevel = 0;
    bool noFlowFlag = false;
    bool helpWanted = false;
    try {
	getopt(args,
	       "job", &jobName,
	       "verbosity", &verbosityLevel,
	       "no-flow", &noFlowFlag,
	       "help", &helpWanted
	       );
    } catch (Exception e) {
	writeln("Problem parsing command-line options.");
	writeln("Arguments not processed: ");
	args = args[1 .. $]; // Dispose of program name in first argument.
	foreach (myarg; args) writeln("    arg: ", myarg);
	write(msg);
	exit(1);
    }
    if (helpWanted) {
	write(msg);
	exit(0);
    }
    if (jobName.length == 0) {
	writeln("Need to specify a job name.");
	write(msg);
	exit(1);
    }
    if (verbosityLevel > 0) {
	writeln("Begin preparation with command-line arguments.");
	writeln("  jobName: ", jobName);
	writeln("  noFlowFlag: ", noFlowFlag);
	writeln("  verbosityLevel: ", verbosityLevel);
    }

    writeln("Start up LuaD connection.");
    auto lua = new LuaState;
    lua.openLibs();
    registerVector3(lua);
    registerFlowState(lua);
    lua.doFile("e4prep.lua");
    lua.doFile(jobName~".lua");
    writeln("[TODO] Actually use the configuration data to make the grids, etc");
    if (!noFlowFlag) {
	writeln("[TODO] Write the flow data files.");
    }
    writeln("Done.");
} // end main()


