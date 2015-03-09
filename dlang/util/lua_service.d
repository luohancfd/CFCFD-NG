/**
 * lua_service.d
 * Home to some commonly used idioms (collected as functions)
 * when using interfacing with Lua.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-01-14
 */

module util.lua_service;

import std.c.stdlib : exit;
import std.stdio;
import std.string;
import std.conv;
import luad.all;
import luad.c.lua;
import luad.c.lauxlib;

LuaState initLuaState(string fname)
{
    auto lua = new LuaState();
    lua.openLibs();
    lua.doFile(fname);
    return lua;
}

/++
 + Grab a set of values of the same type based on a set of (string) keys.
 +/
void getValues(T)(LuaTable t, in string[] keys, out T[string] values, string tabName)
{
    foreach ( k; keys ) {
	try {
	    values[k] = t.get!T(k);
	}
	catch ( Exception e ) {
	    string msg = format("The key: %s could not be found with a useful value of type: %s in table: %s", k, typeid(T), tabName);
	    throw new Exception(msg);
	}
    }
}

/++
 + Grab an array of values of the same type out of a table in array format.
 +
 + Note 'tabName' is just used to report a useful error message.
 + The table is already passed as the first argument.
 +/
void getArray(T)(LuaTable t, out T[] values, string tabName)
{
    auto len = t.length;
    values.length = len;
    try {
	foreach ( i; 1..len+1 ) {
	    // Placed in i-1 because of 0-offset, 1-offset difference
	    // between D and Lua.
	    values[i-1] = t.get!T(i);
	}
    }
    catch ( Exception e ) {
	string msg = format("There was a problem looping over table: %s. An array of values was expected.\n", tabName);
	throw new Exception(msg);
    }
}

/**
 * Get array of numbers from index in Lua stack.
 */

void getArrayOfNumbers(lua_State* L, int index, out double[] values)
{
    auto n = to!int(lua_objlen(L, index));
    foreach ( i; 1..n+1 ) {
	lua_rawgeti(L, index, i);
	if ( lua_isnumber(L, -1) ) values ~= lua_tonumber(L, -1);
	// Silently ignore anything that isn't a value.
	lua_pop(L, 1);
    }
}

/**
 * This creates a new userdata spot on the Lua stack and
 * populates it with an object of type T.
 *
 * Note: we explicitly set the metatable name as
 * a string and do NOT try to get the type name
 * using T.stringof (or a more elaborate __traits function).
 * The reason is that we want control of the name that
 * appears in the Lua script, no matter what name the various D
 * compilers decide to give your class type.
 */
int pushObj(T, string metatableName)(lua_State* L, in T obj)
{
    auto ptr = cast(T*) lua_newuserdata(L, obj.sizeof);
    *ptr = obj.dup();
    luaL_getmetatable(L, metatableName.toStringz);
    lua_setmetatable(L, -2);
    return 1;
}

/**
 * Retrieve reference to an object from lua_State.
 *
 * This function looks for an object of type T
 * at the specified index in lua_State. A error is
 * raised by luaL_checkudata is the object is not of type
 * FlowState.
 */
T checkObj(T, string metatableName)(lua_State* L, int index)
{
    auto ptr = cast(T*) luaL_checkudata(L, index, metatableName.toStringz);
    return *ptr;
}

/**
 * Test if the object type from index position in stack matches tname.
 *
 * Assuming the object is stored as Lua userdata, this returns
 * the metatable name, which is essentially the object's type.
 * If there is no valid object, the string "nil" is returned.
 *
 * This code basically duplicates the function luaL_testudata 
 * in the Lua source file: lauxlib.c. IN LUA VERSION 5.2
 * (You won't find it in our code collection.)
 * The difference is that it only looks to see if the
 * metatable matches.
 */
bool isObjType(lua_State* L, int index, string tname)
{
    bool result;
    void *p = lua_touserdata(L, index);
    if ( p ) {  /* value is a userdata? */
	if (lua_getmetatable(L, index)) {  /* does it have a metatable? */
	    luaL_getmetatable(L, tname.toStringz);  /* get correct metatable */
	    if ( lua_rawequal(L, -1, -2) )  /* the same? */
		result = true;
	    else
		result = false;
	    lua_pop(L, 2);  /* remove both metatables */
	    return result;
	}
    }
    return false;  /* value is not a userdata with a metatable */
}


/**
 * Call an object's toString method and push result on Lua stack.
 */
extern(C) int toStringObj(T, string metatableName)(lua_State* L)
{
    auto path = checkObj!(T, metatableName)(L, 1);
    lua_pushstring(L, toStringz(path.toString));
    return 1;
}

/**
 * Attempt to retrieve a double from a field in a table.
 *
 * We can configure what happens in the event the value is missing or invalid.
 * The boolean errorIfNotFound controls what happens if the field is missing
 * (or nil -- we can't actually tell the difference in Lua). If this is set
 * to true, then we raise a Lua error and print the supplied error message.
 * If we don't care if it's missing, then we simply return the valIfError value.
 *
 * Another common case is that we don't care if the value is missing (so set
 * errorIfNotFound to false) BUT, it is set, then we want to make sure it
 * is valid. In other words, we don't want to ignore an incorrectly set value.
 * In this case, set errorIfNotValid to true. When an invalid value is encountered,
 * we will raise a Lua error using the supplied error message.
 *
 * If we really don't care if the attempt fails, then both bools should be set
 * to false. In this case, the value valIfError is returned and no errors are raised.
 */
double getNumberFromTable(lua_State* L, int index, string field,
			  bool errorIfNotFound=false, double valIfError=double.init,
			  bool errorIfNotValid=false, string errMsg="")
{
    lua_getfield(L, index, field.toStringz);
    if ( lua_isnil(L, -1) ) {
	if ( errorIfNotFound ) {
	    luaL_error(L, errMsg.toStringz);
	}
	else { // We didn't really care
	    return valIfError;
	}
    }
    // Presumably then we have something to look at.
    if ( lua_isnumber(L, -1) ) {
	auto val = lua_tonumber(L, -1);
	lua_pop(L, 1);
	return val;
    }
    // else, failed to find a number value.
    if ( errorIfNotValid ) {
	luaL_error(L, errMsg.toStringz);
    }
    // We didn't want to fail, so give back double.init
    return valIfError;
} // end getNumberFromTable()

int getIntegerFromTable(lua_State* L, int index, string field,
			bool errorIfNotFound=false, int valIfError=int.init,
			bool errorIfNotValid=false, string errMsg="")
{
    lua_getfield(L, index, field.toStringz);
    if ( lua_isnil(L, -1) ) {
	if ( errorIfNotFound ) {
	    luaL_error(L, errMsg.toStringz);
	}
	else { // We didn't really care
	    return valIfError;
	}
    }
    // Presumably then we have something to look at.
    if ( lua_isnumber(L, -1) ) {
	auto val = to!int(lua_tointeger(L, -1));
	lua_pop(L, 1);
	return val;
    }
    // else, failed to find an integer value.
    if ( errorIfNotValid ) {
	luaL_error(L, errMsg.toStringz);
    }
    // We didn't want to fail, so give back int.init
    return valIfError;
} // end getIntegerFromTable()

bool getBooleanFromTable(lua_State* L, int index, string field,
			 bool errorIfNotFound=false, bool valIfError=bool.init,
			 bool errorIfNotValid=false, string errMsg="")
{
    lua_getfield(L, index, field.toStringz);
    if ( lua_isnil(L, -1) ) {
	if ( errorIfNotFound ) {
	    luaL_error(L, errMsg.toStringz);
	}
	else { // We didn't really care
	    return valIfError;
	}
    }
    // Presumably then we have something to look at.
    if ( lua_isboolean(L, -1) ) {
	auto val = lua_toboolean(L, -1);
	lua_pop(L, 1);
	return val;
    }
    // else, failed to find a boolean value.
    if ( errorIfNotValid ) {
	luaL_error(L, errMsg.toStringz);
    }
    // We didn't want to fail, so give back bool.init
    return valIfError;
} // end getBooleanFromTable()

/**
 * Custom exception type for signalling Lua input errors.
 *
 * Recipe for custom execptions used from D Cookbook, p. 21.
 *
 * Reference:
 * Ruppe, A.D. (2014)
 * D Cookbook
 * Packt Publishing, Birmingham, UK
 */
class LuaInputException : Exception {
    this(string message, string file=__FILE__, size_t line=__LINE__,
	 Throwable next=null)
    {
	super(message, file, line, next);
    }
}



unittest
{
    auto lua = new LuaState;
    auto t = lua.newTable();
    t.set!(string,double)("A", 9.0);
    t.set!(string,double)("B", -15.8);
    t.set!(string,int)("C", 2);

    string[2] keys = ["A", "B"];
    string[1] keys2 = ["C"];
    string[3] keys3 = ["A", "B", "C"];

    double[string] vals;

    /// Test 1. Grab A and B as doubles.
    getValues(t, keys, vals, "test");
    assert(vals["A"] == 9.0);
    assert(vals["B"] == -15.8);

    /// Test 2. Grab C as int.
    int[string] vals2;
    getValues(t, keys2, vals2, "test2");
    assert(vals2["C"] == 2);

    /// Test 3. Grab A and B as doubles using slice from keys3
    getValues(t, keys3[0..2], vals, "test3");
    assert(vals["A"] == 9.0);
    assert(vals["B"] == -15.8);
    
    /// Test 4. Grab all values as doubles
    getValues(t, keys3, vals, "test4");
    assert(vals["A"] == 9.0);
    assert(vals["B"] == -15.8);
    assert(vals["C"] == 2.0);

    /// Test 5. Expect an exit exception when we go for an invalid key.
    keys3[0] = "AA";
    try {
	getValues(t, keys3, vals, "test6");
    }
    catch (Exception e) {
	assert(e);
    }
}
