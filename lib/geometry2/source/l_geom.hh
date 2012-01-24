#ifndef L_GEOM_HH
#define L_GEOM_HH

#include <string>

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/lunar.hh"
#include "geom.hh"

class luaVector3 {
public:
    // Lunar requirements:
    static const char className[];
    static Lunar<luaVector3>::RegType member_data[];
    static Lunar<luaVector3>::RegType methods[];
    static Lunar<luaVector3>::MetaType metamethods[];

    std::string str() const;

    // Constructor based on passed in lua_State
    luaVector3(lua_State *L);

    // Constructor based on a real Vector3
    luaVector3(const Vector3 &v);

    // Destructor
    ~luaVector3();

    Vector3& r_object()
    { return v_; }

    // Member data access
    int x(lua_State *L);
    int y(lua_State *L);
    int z(lua_State *L);

    // Member methods
    int norm(lua_State *L);

    // Copy, because Lua only makes references
    int copy(lua_State *L);
private:
    Vector3 v_;
};

bool istype(lua_State *L, int index, const char *tname);
int l_equals(lua_State *L);
int l_add(lua_State *L);
int l_subtract(lua_State *L);
int l_mul(lua_State *L);
int l_unm(lua_State *L);
int l_vabs(lua_State *L);
int l_unit(lua_State *L);
int l_dot(lua_State *L);
int l_cross(lua_State *L);

int open_geom(lua_State *L);

extern "C" int luaopen_geometry(lua_State *L);

#endif
