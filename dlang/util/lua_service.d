/**
 * lua_service.d
 * Home to some commonly used idioms (collected as functions)
 * when using the luaD module.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-01-14
 */

import std.c.stdlib : exit;
import std.stdio;
import std.string;
import luad.all;

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

    /// Test 5. Grab all values as ints
    getValues(t, keys3, vals2, "test5");
    assert(vals2["A"] == 9);
    assert(vals2["B"] == -15);
    assert(vals2["C"] == 2);

    /// Test 6. Expect an exit ception when we go for an invalid key.
    keys3[0] = "AA";
    try {
	getValues(t, keys3, vals, "test6");
    }
    catch (Exception e) {
	assert(e);
    }
}
