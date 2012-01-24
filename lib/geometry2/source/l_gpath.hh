#ifndef L_GPATH_HH
#define L_GPATH_HH

#include <string>

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/lunar.hh"
#include "geom.hh"
#include "gpath.hh"


class luaPath {
public:
    luaPath();
    virtual ~luaPath();
    virtual std::string str() const = 0;

    virtual Path* r_pointer()
    { return path_; }

    virtual luaPath* clone() const = 0;

    // Member data access
    virtual int t0(lua_State *L);
    virtual int t1(lua_State *L);

    virtual int eval(lua_State *L);
    virtual int locate(lua_State *L);
    virtual int dpdt(lua_State *L);
    virtual int length(lua_State *L);
    virtual int partial_length(lua_State *L);
    virtual int point_from_length(lua_State *L);
    virtual int translate(lua_State *L);
    virtual int reverse(lua_State *L);
    virtual int mirror_image(lua_State *L);

protected:
    Path* path_;
};

class luaLine : public luaPath {
public:
    // Lunar requirements:
    static const char className[];
    static Lunar<luaLine>::RegType member_data[];
    static Lunar<luaLine>::RegType methods[];
    static Lunar<luaLine>::MetaType metamethods[];

    std::string str() const;

    luaLine(const luaLine &l);
    luaLine* clone() const;
    
    // Constructor based on passed in lua_State
    luaLine(lua_State *L);

    // Destructor
    ~luaLine();

    // Member data access
    int a(lua_State *L);
    int b(lua_State *L);

    // Special method of string representation
    int e_str(lua_State *L);
    int s_str(lua_State *L);

};

class luaArc : public luaPath {
public:
    // Lunar requirements:
    static const char className[];
    static Lunar<luaArc>::RegType member_data[];
    static Lunar<luaArc>::RegType methods[];
    static Lunar<luaArc>::MetaType metamethods[];
    virtual std::string str() const;
    
    luaArc(const luaArc &l);
    luaArc* clone() const;

    luaArc();
    // Constructor based on passed in lua_State
    luaArc(lua_State *L);

    // Destructor
    virtual ~luaArc();

    // Member data access
    int a(lua_State *L);
    int b(lua_State *L);
    int c(lua_State *L);

};

class luaArc3 : public luaArc {
public:
    // Lunar requirements:
    static const char className[];
    static Lunar<luaArc3>::RegType member_data[];
    static Lunar<luaArc3>::RegType methods[];
    static Lunar<luaArc3>::MetaType metamethods[];
    // Inherit from luaArc:    std::string str() const;
    
    luaArc3(const luaArc3 &l);
    luaArc3* clone() const;

    // Constructor based on passed in lua_State
    luaArc3(lua_State *L);

    // Destructor
    ~luaArc3();
};

class luaBezier : public luaPath {
public:
    // Lunar requirements
    static const char className[];
    static Lunar<luaBezier>::RegType member_data[];
    static Lunar<luaBezier>::RegType methods[];
    static Lunar<luaBezier>::MetaType metamethods[];
    std::string str() const;

    luaBezier(const luaBezier &l);
    luaBezier* clone() const;

    // Constructor based on passed in lua_State
    luaBezier(lua_State *L);

    // Constructor based on a real Vector3
    luaBezier(const Bezier &b);

    // Destructor
    ~luaBezier();

    // Member data access
    int B(lua_State *L);
    int nc(lua_State *L);

    // Member method
    int add_point(lua_State *L);
    int e_str(lua_State *L);
    int s_str(lua_State *L);
};

class luaNurbs : public luaPath {
public:
    // Lunar requirements
    static const char className[];
    static Lunar<luaNurbs>::RegType member_data[];
    static Lunar<luaNurbs>::RegType methods[];
    static Lunar<luaNurbs>::MetaType metamethods[];
    std::string str() const;

    luaNurbs(const luaNurbs &l);
    luaNurbs* clone() const;

    // Constructor based on passed in lua_State
    luaNurbs(lua_State *L);

    // Constructor based on a real Nurbs object
    luaNurbs(const Nurbs &n);

    // Destructor
    ~luaNurbs();

    // Methods
    int knot_insertion(lua_State *L);
    int knot_refinement(lua_State *L);

    // Member data access
    int P(lua_State *L);
    int nc(lua_State *L);
    int w(lua_State *L);
    int nw(lua_State *L);
    int p(lua_State *L);
    int U(lua_State *L);
    int nu(lua_State *L);

};


class luaPolyline : public luaPath {
public:
    // Lunar requirements
    static const char className[];
    static Lunar<luaPolyline>::RegType member_data[];
    static Lunar<luaPolyline>::RegType methods[];
    static Lunar<luaPolyline>::MetaType metamethods[];
    virtual std::string str() const;

    luaPolyline(const luaPolyline &l);
    luaPolyline* clone() const;

    luaPolyline();
    // Constructor based on passed in lua_State
    luaPolyline(lua_State *L);

    // Destructor
    ~luaPolyline();

    // Member data access (for internal [C++] usage)
    //    size_t n_segments();
    //    luaPath* segment(size_t i);
    
    // Member data access (for Lua user access)
    int t_seg(lua_State *L);
    int n_seg(lua_State *L);

    // Member methods
    int add_segment(lua_State *L);
protected:
    std::vector<luaPath*> l_path_;
};

class luaSpline : public luaPolyline {
public:
    // Lunar requirements
    static const char className[];
    static Lunar<luaSpline>::RegType member_data[];
    static Lunar<luaSpline>::RegType methods[];
    static Lunar<luaSpline>::MetaType metamethods[];
    // inherited from luaPolyline: std::string str() const;

    luaSpline(const luaSpline &l);
    luaSpline* clone() const;

    // Constructor based on passed in lua_State
    luaSpline(lua_State *L);

    // Destructor
    ~luaSpline();
};

class luaXPath : public luaPath {
public:
    luaXPath();
    virtual ~luaXPath();
    virtual std::string str() const = 0;

    virtual luaXPath* clone() const = 0;

    virtual XPath* xpath_pointer() const { return dynamic_cast<XPath*>(path_); }

    // Member data access
    virtual int x0(lua_State *L);
    virtual int x1(lua_State *L);

    // Member methods
    virtual int xeval(lua_State *L);
    virtual int dydx(lua_State *L);
    virtual int d2ydx2(lua_State *L);
    virtual int locate_x(lua_State *L);
};

class luaXPoly : public luaXPath {
public:
    // Lunar requirements
    static const char className[];
    static Lunar<luaXPoly>::RegType member_data[];
    static Lunar<luaXPoly>::RegType methods[];
    static Lunar<luaXPoly>::MetaType metamethods[];
    std::string str() const;

    // Member data access
    int B(lua_State *L);
    int n(lua_State *L);

    // Constructor based on passed in lua_State
    luaXPoly(lua_State *L);
    luaXPoly(const luaXPoly &l);
    luaXPoly(const XPoly &x);
    luaXPoly* clone() const;
    // Destructor
    ~luaXPoly();
};

class luaXSpline : public luaXPath {
public:
    // Lunar requirements
    static const char className[];
    static Lunar<luaXSpline>::RegType member_data[];
    static Lunar<luaXSpline>::RegType methods[];
    static Lunar<luaXSpline>::MetaType metamethods[];
    std::string str() const;

    // Constructor based on passed in lua_State
    luaXSpline(lua_State *L);
    luaXSpline(const luaXSpline &l);
    luaXSpline(const XSpline &x);
    luaXSpline* clone() const;
    // Destructor
    ~luaXSpline();

    // Member method
    int locate_dydx(lua_State *L);

private:
    std::vector<luaXPoly*> l_xpoly_;
};

class luaXBezier : public luaXPath {
public:
    // Lunar requirements
    static const char className[];
    static Lunar<luaXBezier>::RegType member_data[];
    static Lunar<luaXBezier>::RegType methods[];
    static Lunar<luaXBezier>::MetaType metamethods[];
    std::string str() const;

    // Constructor based on passed in lua_State
    luaXBezier(lua_State *L);
    luaXBezier(const luaXBezier &l);
    luaXBezier(const XBezier &x);
    luaXBezier* clone() const;
    // Destructor
    ~luaXBezier();
};

class luaXPolyline : public luaXPath {
public:
    // Lunar requirements
    static const char className[];
    static Lunar<luaXPolyline>::RegType member_data[];
    static Lunar<luaXPolyline>::RegType methods[];
    static Lunar<luaXPolyline>::MetaType metamethods[];
    std::string str() const;

    // Constructor based on passed in lua_State
    luaXPolyline(lua_State *L);
    luaXPolyline(const luaXPolyline &l);
    luaXPolyline* clone() const;
    
    // Destructor
    ~luaXPolyline();
private:
    std::vector<luaXPath*> l_xpath_;
};

Path* check_path(lua_State *L, int index);
luaPath* new_l_path(lua_State *L, int index);
XPath* check_xpath(lua_State *L, int index);
luaXPath* new_l_xpath(lua_State *L, int index);
void handle_path_defaults(Path *p, lua_State *L, int index);
void handle_x_params_xpath(lua_State *L, double &x0, double &x1);
void handle_xpath_defaults(Path *p, lua_State *L, int index);
int open_gpath(lua_State *L, int table=LUA_GLOBALSINDEX);

#endif
