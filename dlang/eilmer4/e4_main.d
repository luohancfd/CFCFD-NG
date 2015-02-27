/** e4_main.d
 * Eilmer4 compressible-flow simulation code, top-level function.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-02-05
 */

import std.stdio;
import std.getopt;

import geom;
import gas;
import globalconfig;
import e4_core;

import luad.all;
import luaflowstate;
import luageom;
import luagpath;
import luasurface;
import luaunifunction;
import luasgrid;


void main(string[] args)
{
    writeln("Eilmer4 compressible-flow simulation code.");
    writeln("Revision: PUT_REVISION_STRING_HERE");

    string msg = "Usage:                               Comment:\n";
    msg       ~= "e4shared [--job=<string>]            file names built from this string\n";
    msg       ~= "         [--prep]                    prepare config, grid and flow files\n";
    msg       ~= "         [--run]                     run the simulation over time\n";
    msg       ~= "         [--tindx=<int>]             defaults to 0\n";
    msg       ~= "         [--verbosity=<int>]         defaults to 0\n";
    msg       ~= "         [--max-wall-clock=<int>]    in seconds\n";
    msg       ~= "         [--help]                    writes this message\n";
    if ( args.length < 2 ) {
	writeln("Too few arguments.");
	write(msg);
	exit(1);
    }
    string jobName = "";
    int verbosityLevel = 0;
    bool prepFlag = false;
    bool runFlag = false;
    int tindxStart = 0;
    int maxWallClock = 5*24*3600; // 5 days default
    bool helpWanted = false;
    try {
	getopt(args,
	       "job", &jobName,
	       "verbosity", &verbosityLevel,
	       "prep", &prepFlag,
	       "run", &runFlag,
	       "tindx", &tindxStart,
	       "max-wall-clock", &maxWallClock,
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
    if (prepFlag) {
	writeln("Start LuaD connection.");
	auto lua = new LuaState;
	lua.openLibs();
	registerVector3(lua);
	registerFlowState(lua);
	registerPaths(lua);
	registerSurfaces(lua);
	registerUnivariateFunctions(lua);
	registerStructuredGrid(lua);
	lua.doFile("e4prep.lua");
	lua.doFile(jobName~".lua");
	lua.doString("build_job_files(\""~jobName~"\")");
	writeln("Done.");
    }
    if (runFlag) {
	GlobalConfig.base_file_name = jobName;
	if (verbosityLevel > 0) {
	    writeln("Begin simulation with command-line arguments.");
	    writeln("  jobName: ", jobName);
	    writeln("  tindxStart: ", tindxStart);
	    writeln("  maxWallClock: ", maxWallClock);
	    writeln("  verbosityLevel: ", verbosityLevel);
	}
	GlobalConfig.verbosity_level = verbosityLevel;
	
	double sim_time = init_simulation(tindxStart);
	writeln("starting simulation time= ", sim_time);
	sim_time = integrate_in_time(GlobalConfig.max_time, maxWallClock);
	finalize_simulation(sim_time);
    }
    writeln("Done.");
} // end main()


