import std.stdio;
import std.string;
import std.conv;
import std.math;
import std.algorithm;
import luad.all;
import util.lua_service;

import fileutil;
import readconfig;
import globalconfig;
import globaldata;

void main(string[] args)
{
    writeln("args.length= ", args.length);
    string jobName = args[1];
    string luafname = args[2];
    writeln("jobName= ", jobName, " luafname= ", luafname);
    int tindx = 30;
    
    GlobalConfig.base_file_name = jobName;
    read_config_file();
    
    auto blk = mySolidBlocks[0];
    blk.assembleArrays();
    blk.bindFacesAndVerticesToCells();
    blk.readGrid(make_file_name!"solid-grid"(jobName, blk.id, 0));
    blk.readSolution(make_file_name!"solid"(jobName, blk.id, tindx));

    auto lua = initLuaState(luafname);
    auto lfunc = lua.get!LuaFunction("Ts");

    double LInf_norm = 0.0;
    double sum = 0.0;
    int N = 0;

    foreach (cell; blk.activeCells) {
	double T_ex = T_analytical(lfunc, cell.pos.x, cell.pos.y);
	double abs_diff = fabs(cell.T[0] - T_ex);
	LInf_norm = max(LInf_norm, abs_diff);
	sum += abs_diff^^2;
	N += 1;
    }
    double L2_norm = sqrt(sum/N);

    writeln(format("L2-norm = %20.12e", L2_norm));
    writeln(format("L-inf-norm = %20.12e", LInf_norm));
}

double T_analytical(LuaFunction lfunc, double x, double y)
{
    LuaObject[] ret = lfunc(x, y);
    return ret[0].to!double();
}
