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


void main(string[] args)
{
    writeln("Eilmer4 compressible-flow simulation code.");
    writeln("Revision: PUT_REVISION_STRING_HERE");

    string msg = "Usage:                                 Comment:\n";
    msg       ~= "e4shared [--job=<string>] [--run]      file names built from this string\n";
    msg       ~= "         [--tindx=<int>]               defaults to 0\n";
    msg       ~= "         [--verbosity=<int>]           defaults to 0\n";
    msg       ~= "         [--max-wall-clock=<int>]      in seconds\n";
    msg       ~= "         [--help]                      writes this message\n";
    if ( args.length < 2 ) {
	writeln("Too few arguments.");
	write(msg);
	exit(1);
    }
    string jobName = "";
    int verbosityLevel = 0;
    bool runFlag = false;
    int tindxStart = 0;
    int maxWallClock = 5*24*3600; // 5 days default
    bool helpWanted = false;
    try {
	getopt(args,
	       "job", &jobName,
	       "verbosity", &verbosityLevel,
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
    if (runFlag) {
	if ( jobName.length == 0 ) {
	    writeln("Need to specify a job name.");
	    write(msg);
	    exit(1);
	}
	GlobalConfig.base_file_name = jobName;
	if ( verbosityLevel > 0 ) {
	    writeln("jobName: ", jobName);
	    writeln("tindxStart: ", tindxStart);
	    writeln("maxWallClock: ", maxWallClock);
	}
	GlobalConfig.verbosity_level = verbosityLevel;
	
	writeln("Begin simulation...");
	double sim_time = init_simulation();
	writeln("starting simulation time= ", sim_time);
	sim_time = integrate_in_time(sim_time);
	finalize_simulation(sim_time);
    }
    writeln("Done.");
} // end main()


