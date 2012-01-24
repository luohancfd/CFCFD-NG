// Author: Rowan J. Gollan
// Date: 15-Jun-2009
// Place: NASA Langley, Hampton, Virginia, USA
//

#ifndef L_SURFACE_HH
#define L_SURFACE_HH

#include <string>

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/lunar.hh"
#include "surface.hh"

class luaSurface {
public:
    luaSurface();
    virtual ~luaSurface();
    virtual std::string str() const = 0;

    virtual ParametricSurface* r_pointer()
    { return surf_; }

    virtual luaSurface* clone() const = 0;

    // member data access
    virtual int r0(lua_State *L);
    virtual int r1(lua_State *L);
    virtual int s0(lua_State *L);
    virtual int s1(lua_State *L);

    // member methods
    virtual int eval(lua_State *L);
    virtual int dpdr(lua_State *L);
    virtual int dpds(lua_State *L);
    virtual int translate(lua_State *L);
    virtual int mirror_image(lua_State *L);

protected:
    ParametricSurface *surf_;

};

class luaNurbsSurface : public luaSurface {
public:
    static const char className[];
    static Lunar<luaNurbsSurface>::RegType member_data[];
    static Lunar<luaNurbsSurface>::RegType methods[];
    static Lunar<luaNurbsSurface>::MetaType metamethods[];

    std::string str() const;

    luaNurbsSurface* clone() const;

    luaNurbsSurface(lua_State *L);
    luaNurbsSurface(const luaNurbsSurface &n);
    luaNurbsSurface(const NurbsSurface &n);

    ~luaNurbsSurface();

//    int write_IGES_file(lua_State *L);

};

int open_surface(lua_State *L, int table);

#endif
