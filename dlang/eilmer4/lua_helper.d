// lua_helper.d
//
// A place to put some frequently used functions
// for interacting with the Lua stack. Some of the functions
// are for interacting in the D code. Other functions
// (marked extern(C)) are functions made available in the
// Lua script.
//
// RG & PJ 2015-03-17 -- First hack (with Guiness in hand)

import luad.all;
import luad.stack;
import luad.c.lua;
import fvcell;
import luaflowstate;
import globaldata;

// -----------------------------------------------------
// Convenience functions for user's Lua script

// [TODO] [FIXME] the following function won't work in parallel loops
// because the gasBlocks array probably won't be initialized correctly
// for any thread other than the main thread.
extern(C) int luafn_sampleFlow(lua_State *L)
{
    // Get arguments from lua_stack
    auto blkId = lua_tointeger(L, 1);
    auto i = lua_tointeger(L, 2);
    auto j = lua_tointeger(L, 3);
    auto k = lua_tointeger(L, 4);
    
    // Grab the appropriate cell
    auto cell = gasBlocks[blkId].get_cell(i, j, k);
    
    // Return the interesting bits as a table.
    lua_newtable(L);
    pushCellToTable(cell, 0, L);
    return 1;
}

// -----------------------------------------------------
// D code functions

/**
 * Push the interesting data from a cell to a LuaD LuaTable.
 */
void pushCellToLuaTable(in FVCell cell, size_t gtl, LuaTable tab)
{
    tab["x"] = cell.pos[gtl].x;
    tab["y"] = cell.pos[gtl].y;
    tab["z"] = cell.pos[gtl].z;
    tab["vol"] = cell.volume[gtl];
    pushFlowStateToLuaTable(cell.fs, tab);
}

/**
 * Push the interesting data from a cell to Lua table at TOS.
 */
void pushCellToTable(in FVCell cell, size_t gtl, lua_State* L)
{
    auto tab = getValue!LuaTable(L, -1);
    pushCellToLuaTable(cell, gtl, tab);
}




