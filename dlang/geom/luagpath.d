/**
 * A Lua interface for the D gpath module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-02-22
 */

module luagpath;

// We cheat to get the C Lua headers by using LuaD.
import luad.all;
import luad.c.lua;
import luad.c.lauxlib;
import std.stdio;
import std.string;
import util.lua_service;
import geom;
import gpath;
import luageom;

immutable string LineMT = "Line"; // Name of Line metatable
immutable string ArcMT = "Arc";
immutable string BezierMT = "Bezier";

// A place to hang on to references to objects that are pushed into the Lua domain.
// We don't want the D garbage collector to prematurely dispose of said objects.
static const(Path)[] pathStore; 

Path checkPath(lua_State* L, int index) {
    if ( isObjType(L, index, LineMT) ) {
	return checkObj!(Line, LineMT)(L, index);
    }
    if ( isObjType(L, index, ArcMT) ) {
	return checkObj!(Arc, ArcMT)(L, index);
    }
    if ( isObjType(L, index, BezierMT) ) {
	return checkObj!(Bezier, BezierMT)(L, index);
    }
    // if all else fails
    return null;
}

extern(C) int opCallPath(T, string MTname)(lua_State* L)
{
    auto path = checkObj!(T, MTname)(L, 1);
    auto t = luaL_checknumber(L, 2);
    auto pt = path(t);
    return pushVector3(L, pt);
}

extern(C) int t0Path(T, string MTname)(lua_State* L)
{
    // Not much error checking here because we assume
    // users are knowing what they are doing if
    // they are messing with the getter/setter functions.
    int narg = lua_gettop(L);
    auto path = checkObj!(T, MTname)(L, 1);
    if ( narg == 1 ) { // This is a getter
	lua_pushnumber(L, path.t0);
	return 1;
    }
    // else: treat as a setter
    path.t0 = luaL_checknumber(L, 2);
    return 0;
}

extern(C) int t1Path(T, string MTname)(lua_State* L)
{
    // Not much error checking here because we assume
    // users are knowing what they are doing if
    // they are messing with the getter/setter functions.
    int narg = lua_gettop(L);
    auto path = checkObj!(T, MTname)(L, 1);
    if ( narg == 1 ) { // This is a getter
	lua_pushnumber(L, path.t1);
	return 1;
    }
    // else: treat as a setter
    path.t1 = luaL_checknumber(L, 2);
    return 0;
}

extern(C) int copyPath(T, string MTname)(lua_State* L)
{
    // Sometimes it's convenient to get a copy of a path.
    auto path = checkObj!(T, MTname)(L, 1);
    pathStore ~= pushObj!(T, MTname)(L, path);
    return 1;
}

/* ----------------- Path-specific functions --------------- */

/**
 * The Lua constructor for a Line.
 *
 * Example construction in Lua:
 * ---------------------------------
 * a = Vector3:new{}
 * b = Vector3:new{1, 1}
 * ab = Line:new{a, b}
 * ab = Line:new{a, b, t0=0.0, t1=1.0}
 * ---------------------------------
 */
extern(C) int newLine(lua_State* L)
{
    lua_remove(L, 1); // remove first agurment "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to Line:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 at position 1.
    lua_rawgeti(L, 1, 1);
    auto a = checkVector3(L, -1);
    if ( a is null ) {
	string errMsg = `Error in call to Line:new{}.
A Vector3 object is expected as the first argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    // Expect Vector3 at position 2.
    lua_rawgeti(L, 1, 2);
    auto b = checkVector3(L, -1);
    if ( b is null ) {
	string errMsg = `Error in call to Line:new{}.
A Vector3 object is expected as the second argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    string errMsgTmplt = `Error in call to Line:new{}.
A valid value for '%s' was not found in list of arguments.
The value, if present, should be a number.`;
    double t0 = getNumberFromTable(L, 1, "t0", false, 0.0, true, format(errMsgTmplt, "t0"));
    double t1 = getNumberFromTable(L, 1, "t1", false, 1.0, true, format(errMsgTmplt, "t1"));
    auto ab = new Line(*a, *b, t0, t1);
    pathStore ~= pushObj!(Line, LineMT)(L, ab);
    return 1;
} // end newLine()


/**
 * The Lua constructor for an Arc.
 *
 * Example construction in Lua:
 * ---------------------------------
 * a = Vector3:new{1, 0}
 * b = Vector3:new{0, 1}
 * c = Vector3:new{0, 0}
 * abc = Arc:new{a, b, c}
 * abc = Arc:new{a, b, c, t0=0.0, t1=1.0}
 * ---------------------------------
 */
extern(C) int newArc(lua_State* L)
{
    lua_remove(L, 1); // remove first agurment "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to Arc:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 at position 1.
    lua_rawgeti(L, 1, 1);
    auto a = checkVector3(L, -1);
    if ( a is null ) {
	string errMsg = `Error in call to Arc:new{}.
A Vector3 object is expected as the first argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    // Expect Vector3 at position 2.
    lua_rawgeti(L, 1, 2);
    auto b = checkVector3(L, -1);
    if ( b is null ) {
	string errMsg = `Error in call to Arc:new{}.
A Vector3 object is expected as the second argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    // Expect Vector3 at position 3.
    lua_rawgeti(L, 1, 3);
    auto c = checkVector3(L, -1);
    if ( c is null ) {
	string errMsg = `Error in call to Arc:new{}.
A Vector3 object is expected as the third argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    string errMsgTmplt = `Error in call to Arc:new{}.
A valid value for '%s' was not found in list of arguments.
The value, if present, should be a number.`;
    double t0 = getNumberFromTable(L, 1, "t0", false, 0.0, true, format(errMsgTmplt, "t0"));
    double t1 = getNumberFromTable(L, 1, "t1", false, 1.0, true, format(errMsgTmplt, "t1"));
    auto abc = new Arc(*a, *b, *c, t0, t1);
    pathStore ~= pushObj!(Arc, ArcMT)(L, abc);
    return 1;
} // end newArc()


/**
 * The Lua constructor for a Bezier.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{1, 0}
 * p1 = Vector3:new{0.7071, 0.7071}
 * p2 = Vector3:new{0, 1}
 * -- For an arbitrary number of points in the table.
 * bez = Bezier:new{p0, p1, p2}
 * bez = Bezier:new{p0, p1, p2, t0=0.0, t1=1.0}
 * ---------------------------------
 */
extern(C) int newBezier(lua_State* L)
{
    lua_remove(L, 1); // remove first agurment "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to Bezier:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 objects at array positions.
    Vector3[] B;
    int position = 1;
    while ( true ) {
	lua_rawgeti(L, 1, position);
	if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
	auto a = checkVector3(L, -1);
	lua_pop(L, 1);
	if ( a is null ) break;
	B ~= *a;
	++position;
    }
    if ( B.length == 0 ) {
	string errMsg = `Error in call to Bezier:new{}. No valid Vector3 objects found.`;
	luaL_error(L, errMsg.toStringz());
    }
    string errMsgTmplt = `Error in call to Bezier:new{}.
A valid value for '%s' was not found in list of arguments.
The value, if present, should be a number.`;
    double t0 = getNumberFromTable(L, 1, "t0", false, 0.0, true, format(errMsgTmplt, "t0"));
    double t1 = getNumberFromTable(L, 1, "t1", false, 1.0, true, format(errMsgTmplt, "t1"));
    auto bez = new Bezier(B, t0, t1);
    pathStore ~= pushObj!(Bezier, BezierMT)(L, bez);
    return 1;
} // end newBezier()


void registerPaths(LuaState lua)
{
    auto L = lua.state;

    // Register the Line object
    luaL_newmetatable(L, LineMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newLine);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Line, LineMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Line, LineMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(Line, LineMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(Line, LineMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &t0Path!(Line, LineMT));
    lua_setfield(L, -2, "t0");
    lua_pushcfunction(L, &t1Path!(Line, LineMT));
    lua_setfield(L, -2, "t1");

    lua_setglobal(L, LineMT.toStringz);

    // Register the Arc object
    luaL_newmetatable(L, ArcMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newArc);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Arc, ArcMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Arc, ArcMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(Arc, ArcMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(Arc, ArcMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &t0Path!(Arc, ArcMT));
    lua_setfield(L, -2, "t0");
    lua_pushcfunction(L, &t1Path!(Arc, ArcMT));
    lua_setfield(L, -2, "t1");

    lua_setglobal(L, ArcMT.toStringz);

    // Register the Bezier object
    luaL_newmetatable(L, BezierMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newBezier);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Bezier, BezierMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Bezier, BezierMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(Bezier, BezierMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(Bezier, BezierMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &t0Path!(Bezier, BezierMT));
    lua_setfield(L, -2, "t0");
    lua_pushcfunction(L, &t1Path!(Bezier, BezierMT));
    lua_setfield(L, -2, "t1");

    lua_setglobal(L, BezierMT.toStringz);
} // end registerPaths()
    






