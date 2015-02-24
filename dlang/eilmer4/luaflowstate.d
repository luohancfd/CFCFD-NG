/**
 * luaflowstate.d 
 * Lua interface to the FlowState class for users of the prep program.
 *
 * Authors: Rowan G. and Peter J.
 * Version: Initial cut.
 */

module luaflowstate;

import std.stdio;
import std.string;
import std.conv;
import luad.all;
import luad.c.lua;
import luad.c.lauxlib;
import util.lua_service;
import gas;
import flowstate;
import geom;
import luageom;

static GasModel managedGasModel;

/// name for FlowState object in Lua scripts.
immutable string FlowStateMT = "FlowState"; 

// Makes it a little more consistent to make this
// available under this name.
FlowState checkFlowState(lua_State* L, int index)
{
    return checkObj!(FlowState, FlowStateMT)(L, 1);
}

/** 
 * This function implements our constructor for the Lua interface.
 *
 * Construction of a FlowState object from in Lua will accept:
 * -----------------
 * fs = FlowState:new{p=1.0e5, T=300.0, u=1000.0, v=200.0, massf={1.0}}
 * fs = FlowState:new{p=1.0e7, T={300.0}}
 * fs = FlowState:new()
 * fs = FlowState:new{}
 * -----------------
 * Missing velocity components are set to 0.0.
 * Missing mass fraction list is set to {1.0}.
 * For one-temperature gas models, single value for T is OK.
 * Temperature will always be accepted as an array.
 * For all other missing components, the values
 * are the defaults as given by the first constructor
 * in flowstate.d
 * The empty constructors forward through to PJ's
 * constructor that accepts a GasModel argument only.
 */
extern(C) int newFlowState(lua_State* L)
{
    if ( managedGasModel is null ) {
	string errMsg = `Error in call to FlowState:new.
It appears that you have not yet set the GasModel.
Be sure to call setGasModel(fname) before using a FlowState object.`;
	luaL_error(L, errMsg.toStringz);
    }

    lua_remove(L, 1); // Remove first argument "this".
    FlowState fs;

    int narg = lua_gettop(L);
    if ( narg == 0 ) {
	fs = new FlowState(managedGasModel);
	return pushObj!(FlowState, FlowStateMT)(L, fs);
    }
    // else narg >= 1
    if ( !lua_istable(L, 1) ) {
	string errMsg = "Error in call to FlowState:new. A table is expected as first argument.";
	luaL_error(L, errMsg.toStringz);
    }
    // Now we are committed to using the first constructor
    // in class FlowState. So we have to find at least
    // a pressure and temperature(s).
    string errMsg = `Error in call to FlowState:new.
A valid pressure value 'p' is not found in arguments.
The value should be a number.`;
    double p = getNumberFromTable(L, 1, "p", true, double.init, true, errMsg);

    // Next test for T and if it's a single value or array.
    double T[];
    lua_getfield(L, 1, "T");
    if ( lua_isnumber(L, -1 ) ) {
	double Tval = lua_tonumber(L, -1);
	foreach ( i; 0..managedGasModel.n_modes ) T ~= Tval;
    }
    else if ( lua_istable(L, -1 ) ) {
	getArrayOfNumbers(L, -1, T);
	if ( T.length != managedGasModel.n_modes ) {
	    errMsg = "Error in call to FlowState:new.";
	    errMsg ~= "Length of T vector does not match number of modes in gas model.";
	    errMsg ~= format("T.length= %d; n_modes= %d\n", T.length, managedGasModel.n_modes);
	    throw new LuaInputException(errMsg);
	}
    }
    else  {
	errMsg = "Error in call to FlowState:new.";
	errMsg ~= "A valid temperature value 'T' is not found in arguments.";
	errMsg ~= "It should be listed as a single value, or list of values.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);
    // Now everything else is optional. If it has been set, then we will 
    // ensure that it can be retrieved correctly, or signal the user.
    // Values related to velocity.
    double u = 0.0;
    double v = 0.0;
    double w = 0.0;
    string errMsgTmplt = "Error in call to FlowState:new.\n";
    errMsgTmplt ~= "A valid value for '%s' is not found in arguments.\n";
    errMsgTmplt ~= "The value, if present, should be a number.";
    u = getNumberFromTable(L, 1, "u", false, 0.0, true, format(errMsgTmplt, "u"));
    v = getNumberFromTable(L, 1, "v", false, 0.0, true, format(errMsgTmplt, "v"));
    w = getNumberFromTable(L, 1, "w", false, 0.0, true, format(errMsgTmplt, "w"));
    auto vel = Vector3(u, v, w);

    // Values related to mass fractions.
    double massf[];
    lua_getfield(L, 1, "massf");
    if ( lua_isnil(L, -1) ) {
	massf ~= 1.0; // Nothing set, so massf = [1.0]
    }
    else if ( lua_istable(L, -1) ) {
	getArrayOfNumbers(L, -1, massf);
	if ( massf.length != managedGasModel.n_species ) {
	    errMsg = "Error in call to FlowState:new.";
	    errMsg ~= "Length of massf vector does not match number of species in gas model.";
	    errMsg ~= format("massf.length= %d; n_species= %d\n", massf.length, managedGasModel.n_species);
	    throw new LuaInputException(errMsg);
	}
    }
    else  {
	errMsg = "Error in call to FlowState:new.";
	errMsg ~= "A field for mass fractions was found, but the contents are not valid.";
	errMsg ~= "The mass fraction should be given as a list of values.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    // Value for quality
    double quality = getNumberFromTable(L, 1, "quality", false, 1.0, true, format(errMsgTmplt, "quality"));
    
    // Values for B (magnetic field?)
    double Bx = 0.0;
    double By = 0.0;
    double Bz = 0.0;
    Bx = getNumberFromTable(L, 1, "Bx", false, 0.0, true, format(errMsgTmplt, "Bx"));
    By = getNumberFromTable(L, 1, "By", false, 0.0, true, format(errMsgTmplt, "By"));
    Bz = getNumberFromTable(L, 1, "Bz", false, 0.0, true, format(errMsgTmplt, "Bz"));
    auto B = Vector3(Bx, By, Bz);

    // Values related to k-omega model.
    double tke = getNumberFromTable(L, 1, "tke", false, 0.0, true, format(errMsgTmplt, "tke"));
    double omega = getNumberFromTable(L, 1, "omega", false, 0.0, true, format(errMsgTmplt, "omega"));
    double mu_t = getNumberFromTable(L, 1, "mu_t", false, 0.0, true, format(errMsgTmplt, "mu_t"));
    double k_t = getNumberFromTable(L, 1, "k_t", false, 0.0, true, format(errMsgTmplt, "k_t"));

    // We won't let user set 'S' -- shock detector value.
    int S = 0;

    fs = new FlowState(managedGasModel, p, T, vel, massf, quality, B,
		       tke, omega, mu_t, k_t);
    return pushObj!(FlowState, FlowStateMT)(L, fs);
}

/**
 * Provide a peek into the FlowState data as a Lua table.
 *
 * Basically, this gives the user a table to look at the values
 * in a FlowState in a read-only manner. (Well, in truth, the
 * table values can be changed, but they won't be reflected
 * in the FlowState object. This is consistent with the methods
 * of the FlowState object. Presently, there is no automatic
 * update if one fiddles with the gas properties in FlowState
 * object.
 */

string pushGasVar(string var)
{
    return `lua_pushnumber(L, fs.gas.` ~ var ~ `);
lua_setfield(L, -2, "` ~ var ~`");`;
}

string pushGasVarArray(string var)
{
    return `lua_newtable(L);
foreach ( i, val; fs.gas.` ~ var ~ `) {
    lua_pushnumber(L, val); lua_rawseti(L, -2,to!int(i+1));
}
lua_setfield(L, -2, "` ~ var ~`");`;
}

string pushFSVar(string var)
{
return `lua_pushnumber(L, fs.` ~ var ~ `);
lua_setfield(L, -2, "` ~ var ~`");`;
}

string pushFSVecVar(string var)
{
return `pushVector3(L, fs.` ~ var ~ `);
lua_setfield(L, -2, "` ~ var ~`");`;
}


extern(C) int toTable(lua_State* L)
{
    auto fs = checkFlowState(L, 1);
    lua_newtable(L); // anonymous table { }

    lua_newtable(L); // This one will hold gas values.
    mixin(pushGasVar("p"));
    mixin(pushGasVar("a"));
    mixin(pushGasVar("rho"));
    mixin(pushGasVar("mu"));
    mixin(pushGasVar("quality"));
    mixin(pushGasVarArray("massf"));
    mixin(pushGasVarArray("T"));
    mixin(pushGasVarArray("e"));
    mixin(pushGasVarArray("k"));
    lua_setfield(L, -2, "gas");

    mixin(pushFSVar("tke"));
    mixin(pushFSVar("omega"));
    mixin(pushFSVar("mu_t"));
    mixin(pushFSVar("k_t"));
    mixin(pushFSVecVar("vel"));
    mixin(pushFSVecVar("B"));

    return 1;
}

extern(C) int setGasModel(lua_State* L)
{
    string fname = to!string(luaL_checkstring(L, 1));
    managedGasModel = init_gas_model(fname);
    lua_pushinteger(L, managedGasModel.n_species);
    lua_pushinteger(L, managedGasModel.n_modes);
    return 2;
    
}

void registerFlowState(LuaState lua)
{
    auto L = lua.state;
    luaL_newmetatable(L, FlowStateMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newFlowState);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &toStringObj!(FlowState, FlowStateMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &toTable);
    lua_setfield(L, -2, "toTable");

    lua_setglobal(L, FlowStateMT.toStringz);

    lua_pushcfunction(L, &setGasModel);
    lua_setglobal(L, "setGasModel");
}

