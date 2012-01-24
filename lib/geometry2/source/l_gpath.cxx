#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/useful.h"
#include "../../util/source/lunar.hh"
#include "geom.hh"
#include "l_geom.hh"
#include "gpath.hh"
#include "l_gpath.hh"

using namespace std;

luaPath::
luaPath()
    : path_(0) {}

luaPath::
~luaPath()
{
    delete path_;
}

int
luaPath::
t0(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg == 0 ) {
	// This is a getter.
	lua_pushnumber(L, path_->t0);
	return 1;
    }
    // else
    // Treat as a setter.
    path_->t0 = luaL_checknumber(L, 1);
    return 0;
}

int
luaPath::
t1(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg == 0 ) {
	// This is a getter.
	lua_pushnumber(L, path_->t1);
	return 1;
    }
    // else
    // Treat as a setter.
    path_->t1 = luaL_checknumber(L, 1);
    return 0;
}

int
luaPath::
eval(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg != 1 ) {
	ostringstream ost;
	ost << "error in call eval():\n"
	    << "   1 argument expected. " << narg << " arguments received.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
    Vector3 tmp = path_->eval(luaL_checknumber(L, 1));
    Lunar<luaVector3>::push(L, new luaVector3(tmp), true);
    return 1;
}


int
luaPath::
locate(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg != 1 ) {
	ostringstream ost;
	ost << "error in call eval():\n"
	    << "   1 argument expected. " << narg << " arguments received.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
    luaVector3* point = Lunar<luaVector3>::check(L, 1);
    int result_flag;
    double t = path_->locate(point->r_object(), result_flag);
    lua_pushnumber(L, t);
    if( result_flag == SUCCESS )
	lua_pushboolean(L, 1);
    else
	lua_pushboolean(L, 0);
    return 2;
}



int
luaPath::
dpdt(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg != 1 ) {
	ostringstream ost;
	ost << "error in call dpdt():\n"
	    << "   1 argument expected. " << narg << " arguments received.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
    Vector3 tmp = path_->dpdt(luaL_checknumber(L, 1));
    Lunar<luaVector3>::push(L, new luaVector3(tmp), true);
    return 1;
}

int
luaPath::
length(lua_State *L)
{
    lua_pushnumber(L, path_->length());
    return 1;
}

int
luaPath::
partial_length(lua_State *L)
{
    int narg = lua_gettop(L);
    
    if( narg != 2 ) {
	ostringstream ost;
	ost << "error in call partial_length():\n"
	    << "   2 arguments expected. " << narg << " argument(s) received.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    double ta = luaL_checknumber(L, 1);
    double tb = luaL_checknumber(L, 2);

    lua_pushnumber(L, path_->partial_length(ta, tb));
    return 1;
}

int
luaPath::
point_from_length(lua_State *L)
{
    int narg = lua_gettop(L);
    
    if( narg != 1 ) {
	ostringstream ost;
	ost << "error in call point_from_length():\n"
	    << "   1 argument expected. " << narg << " arguments received.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    double len = luaL_checknumber(L, 1);
    double t = 0.0;
    Vector3 tmp = path_->point_from_length(len, t);
    Lunar<luaVector3>::push(L, new luaVector3(tmp), true);
    lua_pushnumber(L, t);
    return 2;
}

int
luaPath::
translate(lua_State *L)
{
    int narg = lua_gettop(L);
    
    if( narg == 1 ) {
	// equivalent to: Line::translate(const Vector3 &v)
	luaVector3* tmp = Lunar<luaVector3>::check(L, 1);
	path_->translate(tmp->r_object());
    }
    else if( narg == 3 ) {
	// equivalent to: Line::translate(double vx, double vy, double vz)
	double vx = luaL_checknumber(L, 1);
	double vy = luaL_checknumber(L, 2);
	double vz = luaL_checknumber(L, 3);
	path_->translate(vx, vy, vz);
    }
    else {
	ostringstream ost;
	ost << "error in call to translate():\n"
	    << "   expect 1 arg:  of type Vector3  <OR>\n"
	    << "   expect 3 args: each of type double.\n"
	    << "  " << narg << " argument(s) received.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
    return 0;
}

int
luaPath::
reverse(lua_State *L)
{
    path_->reverse();
    return 0;
}

int
luaPath::
mirror_image(lua_State *L)
{
    int narg = lua_gettop(L);
    
    if( narg == 2 ) {
	luaVector3* point = Lunar<luaVector3>::check(L, 1);
	luaVector3* normal = Lunar<luaVector3>::check(L, 2);
	path_->mirror_image(point->r_object(), normal->r_object());
    }
    else {
	ostringstream ost;
	ost << "error in call to mirror_image():\n"
	    << "   2 arguments expected. " << narg << " argument(s) received.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
    return 0;
}


string
luaLine::
str() const
{
    Line *line = dynamic_cast<Line *>(path_);
    ostringstream ost;
    ost << "Line(" << line->a.str() << ", " << line->b.str() << ", "
	<< "{label=\"" << line->label << "\", t0=" << line->t0 << ", t1=" << line->t1 << "})";
    return ost.str();
}

luaLine::
luaLine(lua_State *L)
    : luaPath()
{
    // Accept a constructor from Lua of the form:
    // > Line(a, b, {label="label", t0=0.0, t1=1.0})
    //
    
    luaVector3* pa = Lunar<luaVector3>::check(L, 1);
    luaVector3* pb = Lunar<luaVector3>::check(L, 2);

    path_ = new Line(pa->r_object(), pb->r_object());

    handle_path_defaults(path_, L, 3);
}

luaLine::
luaLine(const luaLine &l)
    : luaPath()
{
    Line* line = dynamic_cast<Line*>(l.path_);
    path_ = new Line(*line);
}

luaLine*
luaLine::
clone() const
{
    return new luaLine(*this);
}

luaLine::
~luaLine() {}

int
luaLine::
a(lua_State *L)
{
    int narg = lua_gettop(L);

    Line *line = dynamic_cast<Line *>(path_);

    if( narg == 0 ) {
	// This is a getter.
	Lunar<luaVector3>::push(L, new luaVector3(line->a), true);
	return 1;
    }
    // else
    // Treat as a setter.
    luaVector3* pa = Lunar<luaVector3>::check(L, 1);
    line->a = pa->r_object();
    return 0;
}

int
luaLine::
b(lua_State *L)
{
    int narg = lua_gettop(L);

    Line *line = dynamic_cast<Line *>(path_);

    if( narg == 0 ) {
	// This is a getter.
	Lunar<luaVector3>::push(L, new luaVector3(line->b), true);
	return 1;
    }
    // else
    // Treat as a setter.
    luaVector3* pb = Lunar<luaVector3>::check(L, 1);
    line->b = pb->r_object();
    return 0;
}

int
luaLine::
e_str(lua_State *L)
{
    Line *line = dynamic_cast<Line *>(path_);
    ostringstream ost;
    ost << setprecision(3);
    ost << "  1  " << line->length() << endl;
    lua_pushstring(L, ost.str().c_str());
    return 1;
}

int
luaLine::
s_str(lua_State *L)
{
    Line *line = dynamic_cast<Line *>(path_);
    ostringstream ost;
    Vector3 s = line->eval(0.0);
    Vector3 dpdt = line->dpdt(0.0);
    double dydx = dpdt.y/dpdt.x;
    double dzdx = dpdt.z/dpdt.x;
    
    ost << s.x << " " << s.y << " " << s.z << " " << dydx << " " << dzdx << endl;
    lua_pushstring(L, ost.str().c_str());
    return 1;
}

const char luaLine::className[] = "Line";

#define member_data(class, name) {#name, &class::name}

Lunar<luaLine>::RegType luaLine::member_data[] = {
    member_data(luaLine, t0),
    member_data(luaLine, t1),
    member_data(luaLine, a),
    member_data(luaLine, b),
    {0, 0}
};

#define method(class, name) {#name, &class::name}

Lunar<luaLine>::RegType luaLine::methods[] = {
    method(luaLine, eval),
    method(luaLine, locate),
    method(luaLine, dpdt),
    method(luaLine, length),
    method(luaLine, partial_length),
    method(luaLine, point_from_length),
    method(luaLine, translate),
    method(luaLine, reverse),
    method(luaLine, mirror_image),
    method(luaLine, e_str),
    method(luaLine, s_str),
    {0, 0}
};

Lunar<luaLine>::MetaType luaLine::metamethods[] = {
    {0, 0}
};

string
luaArc::
str() const
{
    Arc *arc = dynamic_cast<Arc*>(path_);
    ostringstream ost;
    ost << "Arc(" << arc->a.str() << ", " << arc->b.str() << ", " << arc->c.str() << ", "
	<< "{label=\"" << arc->label << "\", t0=" << arc->t0 << ", t1=" << arc->t1 << "})";\
    return ost.str();
}

luaArc::
luaArc()
    : luaPath() {}

luaArc::
luaArc(lua_State *L)
    : luaPath()
{
    // Accept a constructor from Lua of the form:
    // > Arc(a, b, c, {label="label", t0=0.0, t1=1.0})
    
    luaVector3* pa = Lunar<luaVector3>::check(L, 1);
    luaVector3* pb = Lunar<luaVector3>::check(L, 2);
    luaVector3* pc = Lunar<luaVector3>::check(L, 3);

    path_ = new Arc(pa->r_object(), pb->r_object(), pc->r_object());

    handle_path_defaults(path_, L, 4);
}

luaArc::
luaArc(const luaArc &l)
    : luaPath()
{
    Arc* arc = dynamic_cast<Arc*>(l.path_);
    path_ = new Arc(*arc);
}

luaArc*
luaArc::
clone() const
{
    return new luaArc(*this);
}

luaArc::
~luaArc() {}

int
luaArc::
a(lua_State *L)
{
    int narg = lua_gettop(L);

    Arc *arc = dynamic_cast<Arc *>(path_);

    if( narg == 0 ) {
	// This is a getter.
	Lunar<luaVector3>::push(L, new luaVector3(arc->a), true);
	return 1;
    }
    // else
    // Treat as a setter.
    luaVector3* pa = Lunar<luaVector3>::check(L, 1);
    arc->a = pa->r_object();
    return 0;
}

int
luaArc::
b(lua_State *L)
{
    int narg = lua_gettop(L);

    Arc *arc = dynamic_cast<Arc *>(path_);

    if( narg == 0 ) {
	// This is a getter.
	Lunar<luaVector3>::push(L, new luaVector3(arc->b), true);
	return 1;
    }
    // else
    // Treat as a setter.
    luaVector3* pb = Lunar<luaVector3>::check(L, 1);
    arc->b = pb->r_object();
    return 0;
}


int
luaArc::
c(lua_State *L)
{
    int narg = lua_gettop(L);

    Arc *arc = dynamic_cast<Arc *>(path_);

    if( narg == 0 ) {
	// This is a getter.
	Lunar<luaVector3>::push(L, new luaVector3(arc->c), true);
	return 1;
    }
    // else
    // Treat as a setter.
    luaVector3* pc = Lunar<luaVector3>::check(L, 1);
    arc->c = pc->r_object();
    return 0;
}

const char luaArc::className[] = "Arc";

Lunar<luaArc>::RegType luaArc::member_data[] = {
    member_data(luaArc, t0),
    member_data(luaArc, t1),
    member_data(luaArc, a),
    member_data(luaArc, b),
    member_data(luaArc, c),
    {0, 0}
};

Lunar<luaArc>::RegType luaArc::methods[] = {
    method(luaArc, eval),
    method(luaArc, locate),
    method(luaArc, dpdt),
    method(luaArc, length),
    method(luaArc, partial_length),
    method(luaArc, point_from_length),
    method(luaArc, translate),
    method(luaArc, reverse),
    method(luaArc, mirror_image),
    {0, 0}
};

Lunar<luaArc>::MetaType luaArc::metamethods[] = {
    {0, 0}
};

luaArc3::
luaArc3(lua_State *L)
    : luaArc()
{ 
    // Accept a constructor from Lua of the form:
    // > Arc3(a, b, c, {label="label", t0=0.0, t1=1.0})
    //
    
    luaVector3* pa = Lunar<luaVector3>::check(L, 1);
    luaVector3* pb = Lunar<luaVector3>::check(L, 2);
    luaVector3* pc = Lunar<luaVector3>::check(L, 3);

    path_ = new Arc3(pa->r_object(), pb->r_object(), pc->r_object());

    handle_path_defaults(path_, L, 4);
}

luaArc3::
luaArc3(const luaArc3 &l)
    : luaArc()
{
    Arc3* arc3 = dynamic_cast<Arc3*>(l.path_);
    path_ = new Arc3(*arc3);
}

luaArc3*
luaArc3::
clone() const
{
    return new luaArc3(*this);
}

luaArc3::
~luaArc3() {}
   
const char luaArc3::className[] = "Arc3";

Lunar<luaArc3>::RegType luaArc3::member_data[] = {
    member_data(luaArc3, t0),
    member_data(luaArc3, t1),
    member_data(luaArc3, a),
    member_data(luaArc3, b),
    member_data(luaArc3, c),
    {0, 0}
};

Lunar<luaArc3>::RegType luaArc3::methods[] = {
    method(luaArc3, eval),
    method(luaArc3, locate),
    method(luaArc3, dpdt),
    method(luaArc3, length),
    method(luaArc3, partial_length),
    method(luaArc3, point_from_length),
    method(luaArc3, translate),
    method(luaArc3, reverse),
    method(luaArc3, mirror_image),
    {0, 0}
};

Lunar<luaArc3>::MetaType luaArc3::metamethods[] = {
    {0, 0}
};

string
luaBezier::
str() const
{
    Bezier *bez = dynamic_cast<Bezier*>(path_);
    ostringstream ost;
    ost << "Bezier({";
    for( size_t i =0; i < bez->B.size(); ++i ) ost << bez->B[i].str() << ", ";
    ost << "},\n{label=\"" << bez->label << "\", t0=" << bez->t0 << ", t1=" << bez->t1 << "})";
    return ost.str();
}

luaBezier::
luaBezier(lua_State *L)
    : luaPath()
{ 
    // Accept a constructor from Lua of the form:
    // > Bezier( {b1, b2, b3}, {label="label", t0=0.0, t1=1.0})
    //

    vector<Vector3> tmp;
    if( !lua_istable(L, 1) ) {
	cout << "Table of Vector3 objects expected.\n";
	exit(1);
    }

    int n = lua_objlen(L, 1);
    for( int i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 1, i);
	luaVector3* v = Lunar<luaVector3>::check(L, -1);
	tmp.push_back(v->r_object());
	lua_pop(L, 1);
    }

    path_ = new Bezier(tmp);

    handle_path_defaults(path_, L, 2);
}

luaBezier::
luaBezier(const Bezier &b)
    : luaPath()
{
    path_ = new Bezier(b);
}

luaBezier::
luaBezier(const luaBezier &l)
    : luaPath()
{
    Bezier* bez = dynamic_cast<Bezier*>(l.path_);
    path_ = new Bezier(*bez);
}

luaBezier*
luaBezier::
clone() const
{
    return new luaBezier(*this);
}

luaBezier::
~luaBezier() {}

int
luaBezier::
B(lua_State *L)
{
    int narg = lua_gettop(L);

    Bezier *bez = dynamic_cast<Bezier*>(path_);

    if( narg == 1 ) {
	// This is a getter.
	int i = luaL_checkint(L, 1);
	Lunar<luaVector3>::push(L, new luaVector3(bez->B[i]), true);
	return 1;
    }
    // else
    // Treat as a setter.
    int i = luaL_checkint(L, 1);
    luaVector3* pb = Lunar<luaVector3>::check(L, 2);
    bez->B[i] = pb->r_object();
    return 0;
}

int
luaBezier::
nc(lua_State *L)
{
    // Treat as getter only.
    Bezier *bez = dynamic_cast<Bezier*>(path_);
    lua_pushinteger(L, bez->B.size());
    return 1;
}

int
luaBezier::
add_point(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg != 1 ) {
	cout << "Error to call add_point():\n";
	cout << "1 argument expected.\n";
	cout << narg << " arguments received.\n";
	cout << "Bailing out!\n";
	exit(1);
    }
    luaVector3* pb = Lunar<luaVector3>::check(L, 1);
    Bezier *bez = dynamic_cast<Bezier*>(path_);
    bez->B.push_back(pb->r_object());
    return 0;
}

int
luaBezier::
e_str(lua_State *L)
{
    Bezier *bez = dynamic_cast<Bezier*>(path_);
    
    Vector3 s = bez->B[0];
    Vector3 e = bez->B[3];
    Vector3 dpdt = bez->dpdt(1.0);
    double dydxB = dpdt.y/dpdt.x;
    double dzdxB = dpdt.z/dpdt.x;
    double p1 = (bez->B[1].x - s.x)/(e.x - s.x);
    double p2 = (bez->B[2].x - s.x)/(e.x - s.x);

    ostringstream ost;
    ost << "2  " << e.x << " " << e.y << " " << e.z << " ";
    ost << dydxB << " " << dzdxB << " " << p1 << " " << p2 << endl;
    lua_pushstring(L, ost.str().c_str());
    return 1;
}

int
luaBezier::
s_str(lua_State *L)
{
    Bezier *bez = dynamic_cast<Bezier*>(path_);
    
    Vector3 s = bez->B[0];
    Vector3 dpdt = bez->dpdt(0.0);
    double dydx = dpdt.y/dpdt.x;
    double dzdx = dpdt.z/dpdt.x;

    ostringstream ost;
    ost << s.x << " " << s.y << " " << s.z << " " << dydx << " " << dzdx << endl;
    lua_pushstring(L, ost.str().c_str());
    return 1;
}
const char luaBezier::className[] = "Bezier";

Lunar<luaBezier>::RegType luaBezier::member_data[] = {
    member_data(luaBezier, t0),
    member_data(luaBezier, t1),
    {0, 0}
};

Lunar<luaBezier>::RegType luaBezier::methods[] = {
    method(luaBezier, B),
    method(luaBezier, nc),
    method(luaBezier, add_point),
    method(luaBezier, eval),
    method(luaBezier, locate),
    method(luaBezier, dpdt),
    method(luaBezier, length),
    method(luaBezier, partial_length),
    method(luaBezier, point_from_length),
    method(luaBezier, translate),
    method(luaBezier, reverse),
    method(luaBezier, mirror_image),
    method(luaBezier, e_str),
    method(luaBezier, s_str),
    {0, 0}
};

Lunar<luaBezier>::MetaType luaBezier::metamethods[] = {
    {0, 0}
};

string
luaNurbs::
str() const
{
    Nurbs *n = dynamic_cast<Nurbs*>(path_);
    ostringstream ost;
    ost << "Nurbs({";
    for( size_t i = 0; i < n->P.size(); ++i ) ost << n->P[i].str() << ", ";
    ost << "},\n";
    ost << "{";
    for( size_t i = 0; i < n->w.size(); ++i ) ost << n->w[i] << ", ";
    ost << "},\n";
    ost << n->p << ",\n";
    ost << "{";
    for( size_t i = 0; i < n->U.size(); ++i ) ost << n->U[i] << ", ";
    ost << "},\n";
    ost << "{label=\"" << n->label << "\", t0=" << n->t0 << ", t1=" << n->t1 << "}";
    if ( !n->inf_list.empty() ) {
	ost << ",\ninf_list={";
	for( size_t i = 0; i < n->inf_list.size(); ++i ) ost << n->inf_list[i] << ", ";
	ost << "})\n";
    }
    else {
	ost << ")\n";
    }

    return ost.str();
}


luaNurbs::
luaNurbs(const luaNurbs &l)
    : luaPath()
{
    Nurbs* n = dynamic_cast<Nurbs*>(l.path_);
    path_ = new Nurbs(*n);
}

luaNurbs*
luaNurbs::
clone() const
{
    return new luaNurbs(*this);
}

luaNurbs::
luaNurbs(lua_State *L)
    : luaPath()
{ 
    // Accept a constructor from Lua of the form:
    // > Nurbs( {b1, b2, b3}, {w1, w2, w3}, p, {U1, .., Um},
    //          {label="label", t0=0.0, t1=1.0, inf_list={2, 3, 7}})
    //

    vector<Vector3> P;
    if( !lua_istable(L, 1) ) {
	ostringstream ost;
	ost << "error in Nurbs constructor:\n"
	    << "   Table of Vector3 objects expected as first argument: control points.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    size_t n = lua_objlen(L, 1);
    for ( size_t i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 1, i);
	luaVector3* v = Lunar<luaVector3>::check(L, -1);
	P.push_back(v->r_object());
	lua_pop(L, 1);
    }

    if ( !lua_istable(L, 2) ) {
	ostringstream ost;
	ost << "error in Nurbs constructor:\n"
	    << "   Table of numbers expected as second argument: weights.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    vector<double> w;
    n = lua_objlen(L, 2);
    for ( size_t i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 2, i);
	w.push_back(luaL_checknumber(L, -1));
	lua_pop(L, 1);
    }

    int p = luaL_checkint(L, 3);

    if ( !lua_istable(L, 4) ) {
	ostringstream ost;
	ost << "error in Nurbs constructor:\n"
	    << "   Table of numbers expected as fourth argument: knots.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    vector<double> U;
    n = lua_objlen(L, 4);
    for ( size_t i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 4, i);
	U.push_back(luaL_checknumber(L, -1));
	lua_pop(L, 1);
    }


    if ( lua_istable(L, 5) ) {
	lua_getfield(L, 5, "inf_list");
	vector<int> inf_list;
	if ( lua_istable(L, -1) ) {
	    n = lua_objlen(L, -1);
	    vector<int> inf_list;
	    for ( size_t i = 1; i <= n; ++i ) {
		lua_rawgeti(L, -1, i);
		inf_list.push_back(luaL_checkint(L, -1));
		lua_pop(L, 1);
	    }
	    path_ = new Nurbs(P, w, p, U, inf_list);
	}
	lua_pop(L, 1);
    }

    if ( path_ == 0 ) { // hasn't been set, no infinite control points were found
	path_ = new Nurbs(P, w, p, U);
    }

    handle_path_defaults(path_, L, 5);
}

luaNurbs::
luaNurbs(const Nurbs &n)
    : luaPath()
{
    path_ = new Nurbs(n);
}

luaNurbs::
~luaNurbs() {}

int
luaNurbs::
knot_insertion(lua_State *L)
{
    int narg = lua_gettop(L);

    if ( narg < 4 ) {
	ostringstream ost;
	ost << "error in Nurbs:knot_insertion():\n"
	    << "   4 arguments expected. " << narg << " recevied.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
    
    double u = luaL_checknumber(L, 1);
    int k = luaL_checkint(L, 2);
    int s = luaL_checkint(L, 3);
    int r = luaL_checkint(L, 4);

    Nurbs *n = dynamic_cast<Nurbs*>(path_);

    Nurbs n_new = n->knot_insertion(u, k, s, r);

    Lunar<luaNurbs>::push(L, new luaNurbs(n_new), true);

    return 1;
}

int
luaNurbs::
knot_refinement(lua_State *L)
{
    if ( ! lua_istable(L, 1) ) {
	ostringstream ost;
	ost << "error in Nurbs:knot_refinement():\n"
	    << "   A table is expected as the first argument, but wasn't received.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }


    size_t n = lua_objlen(L, 1);

    vector<double> X;
    
    for ( size_t i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 1, i);
	X.push_back(luaL_checknumber(L, -1));
	lua_pop(L, 1);
    }

    Nurbs *n1 = dynamic_cast<Nurbs*>(path_);
    Nurbs n_new = n1->knot_refinement(X);

    Lunar<luaNurbs>::push(L, new luaNurbs(n_new), true);

    return 1;
}

int
luaNurbs::
P(lua_State *L)
{
    Nurbs *n = dynamic_cast<Nurbs*>(path_);

    // Only support as getter.
    int i = luaL_checkint(L, 1);
    Lunar<luaVector3>::push(L, new luaVector3(n->P[i]), true);
    return 1;
}

int
luaNurbs::
nc(lua_State *L)
{
    // Treat as getter only.
    Nurbs *n = dynamic_cast<Nurbs*>(path_);
    lua_pushinteger(L, n->P.size());
    return 1;
}

int
luaNurbs::
w(lua_State *L)
{
    Nurbs *n = dynamic_cast<Nurbs*>(path_);

    // Only support as getter.
    int i = luaL_checkint(L, 1);
    lua_pushnumber(L, n->w[i]);
    return 1;
}

int
luaNurbs::
nw(lua_State *L)
{
    // Treat as getter only.
    Nurbs *n = dynamic_cast<Nurbs*>(path_);
    lua_pushinteger(L, n->w.size());
    return 1;
}

int
luaNurbs::
p(lua_State *L)
{
    // Treat as getter only.
    Nurbs *n = dynamic_cast<Nurbs*>(path_);
    lua_pushinteger(L, n->p);
    return 1;
}

int
luaNurbs::
U(lua_State *L)
{
    Nurbs *n = dynamic_cast<Nurbs*>(path_);

    // Only support as getter.
    int i = luaL_checkint(L, 1);
    lua_pushnumber(L, n->U[i]);
    return 1;
}

int
luaNurbs::
nu(lua_State *L)
{
    // Treat as getter only.
    Nurbs *n = dynamic_cast<Nurbs*>(path_);
    lua_pushinteger(L, n->U.size());
    return 1;
}

const char luaNurbs::className[] = "Nurbs";

Lunar<luaNurbs>::RegType luaNurbs::member_data[] = {
    member_data(luaNurbs, t0),
    member_data(luaNurbs, t1),
    {0, 0}
};

Lunar<luaNurbs>::RegType luaNurbs::methods[] = {
    method(luaNurbs, knot_insertion),
    method(luaNurbs, knot_refinement),
    method(luaNurbs, P),
    method(luaNurbs, nc),
    method(luaNurbs, w),
    method(luaNurbs, nw),
    method(luaNurbs, p),
    method(luaNurbs, U),
    method(luaNurbs, nu),
    method(luaNurbs, eval),
    method(luaNurbs, locate),
    method(luaNurbs, dpdt),
    method(luaNurbs, length),
    method(luaNurbs, partial_length),
    method(luaNurbs, point_from_length),
    method(luaNurbs, translate),
    method(luaNurbs, reverse),
    method(luaNurbs, mirror_image),
    {0, 0}
};

Lunar<luaNurbs>::MetaType luaNurbs::metamethods[] = {
    {0, 0}
};


string
luaPolyline::
str() const
{
    Polyline *pline = dynamic_cast<Polyline*>(path_);
    ostringstream ost;
    ost << "Polyline({";
    for( size_t i = 0; i < l_path_.size(); ++i ) ost << endl << l_path_[i]->str() << ", ";
    ost << "},\n{label=\"" << pline->label << "\", t0=" << pline->t0 << ", t1=" << pline->t1 << "})";
    return ost.str();
}

luaPolyline::
luaPolyline()
    : luaPath() {}

luaPolyline::
luaPolyline(lua_State *L)
    : luaPath()
{
    // Accept two constructors from Lua of the form:
    // > Polyline({path1, path2}, {label="label", t0=0.0, t1=1.0})
    // > Polyline({b1, b2, b3}, {label="label", t0=0.0, t1=1.0})

    if( !lua_istable(L, 1) ) {
	cout << "Table of objects expected.\n";
	exit(1);
    }

    int n = lua_objlen(L, 1);

    if( n >= 1 ) {
	lua_rawgeti(L, 1, 1);
	// Check if we have array of paths or
	// array of points.
		
	if( istype(L, -1, "Vector3") ) {
	    // Array of points
	    vector<Vector3*> tmp;
	    for( int i = 1; i <= n; ++i ) {
		lua_rawgeti(L, 1, i);
		luaVector3* v = Lunar<luaVector3>::check(L, -1);
		tmp.push_back(&(v->r_object()));
		lua_pop(L, 1);
	    }
	    path_ = new Polyline(tmp, 0);
	}
	else {
	    // Array of paths
	    vector<Path*> tmp;
	    for( int i = 1; i <= n; ++i ) {
		lua_rawgeti(L, 1, i);
		Path* p = check_path(L, -1);
		luaPath* lp = new_l_path(L, -1);
		if( p != 0 )
		    tmp.push_back(p);
		if( lp != 0 )
		    l_path_.push_back(lp);
		lua_pop(L, 1);
	    }
	    path_ = new Polyline(tmp);
	}
    }
    else{
	// Pass in empty paths
	vector<Path*> tmp;
	path_ = new Polyline(tmp);
    }

    handle_path_defaults(path_, L, 2);
}

luaPolyline::
luaPolyline(const luaPolyline &l)
    : luaPath()
{
    Polyline* pline = dynamic_cast<Polyline*>(l.path_);
    path_ = new Polyline(*pline);
    
    for( size_t i = 0; i < l.l_path_.size(); ++i ) {
	l_path_.push_back(l.l_path_[i]->clone());
    }
}

luaPolyline*
luaPolyline::
clone() const
{
    return new luaPolyline(*this);
}

luaPolyline::
~luaPolyline()
{
    for( size_t i = 0; i < l_path_.size(); ++i ) delete l_path_[i];
}

int
luaPolyline::
t_seg(lua_State *L)
{
    int narg = lua_gettop(L);

    Polyline *pline = dynamic_cast<Polyline*>(path_);

    // Only treat as a getter.
    if( narg != 1 ) {
	cout << "Error wrong number of arguments passed to\n";
	cout << "luaPolyline::t(). Only 1 argument, an int, expected.\n";
	cout << narg << " received.\n";
	cout << "Bailing out!\n";
	exit(1);
    }

    int i = luaL_checkint(L, 1);
    if( i >= int(pline->seg.size()) )
	lua_pushnumber(L, -1.0);
    else
	lua_pushnumber(L, pline->t_seg[i]);
    return 1;
}

int
luaPolyline::
n_seg(lua_State *L)
{
    int narg = lua_gettop(L);

    Polyline *pline = dynamic_cast<Polyline*>(path_);

    // Only treat as a getter.
    if( narg != 0 ) {
	cout << "Error wrong number of arguments passed to\n";
	cout << "luaPolyline::nseg(). Zero arguments are expected.\n";
	cout << narg << " received.\n";
	cout << "Bailing out!\n";
	exit(1);
    }

    lua_pushnumber(L, pline->seg.size());
    return 1;
}

int
luaPolyline::
add_segment(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg == 0 || narg >= 3 ) {
	cout << "Error wrong number of arguments passed to\n";
	cout << "add_segment. 1 or 2 arguments expected.\n";
	cout << narg << " received.\n";
	cout << "Bailing out!\n";
	exit(1);
    }
    
    Path* p = check_path(L, 1);
    Polyline *pline = dynamic_cast<Polyline*>(p);

    luaPath* lp = new_l_path(L, 1);

    if( narg == 2 ) {
	int dirn = luaL_checkint(L, 2);
	pline->add_segment(p, dirn);
	if( dirn == -1 ) {
	    lp->r_pointer()->reverse();
	}
	
    }
    else {
	pline->add_segment(p);
    }

    l_path_.push_back(lp);

    return 0;
}

const char luaPolyline::className[] = "Polyline";

Lunar<luaPolyline>::RegType luaPolyline::member_data[] = {
    member_data(luaPolyline, t0),
    member_data(luaPolyline, t1),
    {0, 0}
};

Lunar<luaPolyline>::RegType luaPolyline::methods[] = {
    method(luaPolyline, add_segment),
    method(luaPolyline, eval),
    method(luaPolyline, locate),
    method(luaPolyline, dpdt),
    method(luaPolyline, length),
    method(luaPolyline, partial_length),
    method(luaPolyline, point_from_length),
    method(luaPolyline, translate),
    method(luaPolyline, reverse),
    method(luaPolyline, mirror_image),
    method(luaPolyline, t_seg),
    method(luaPolyline, n_seg),
    {0, 0}
};

Lunar<luaPolyline>::MetaType luaPolyline::metamethods[] = {
    {0, 0}
};

luaSpline::
luaSpline(lua_State *L)
    : luaPolyline()
{
    // Accept a constructor from Lua of the form:
    // > Spline( {point1, point2, point3}, {label="label", t0=0.0, t1=1.0})
    //

    vector<Vector3> tmp;
    if( !lua_istable(L, 1) ) {
	cout << "Table of Vector3 objects expected.\n";
	exit(1);
    }

    int n = lua_objlen(L, 1);
    for( int i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 1, i);
	luaVector3* v = Lunar<luaVector3>::check(L, -1);
	tmp.push_back(v->r_object());
	lua_pop(L, 1);
    }

    path_ = new Spline(tmp);
    Spline* spline = dynamic_cast<Spline*>(path_);
 
    for( size_t i = 0; i < spline->seg.size(); ++i ) {
	Bezier* bez = dynamic_cast<Bezier*>(spline->seg[i]);
	l_path_.push_back(new luaBezier(*bez));
    }

    handle_path_defaults(path_, L, 2);

}

luaSpline::
luaSpline(const luaSpline &l)
    : luaPolyline()
{
    Spline* spline = dynamic_cast<Spline*>(l.path_);
    path_ = new Spline(*spline);

    for( size_t i = 0; i < l.l_path_.size(); ++i ) {
	l_path_.push_back(l.l_path_[i]->clone());
    }
}

luaSpline*
luaSpline::
clone() const
{
    return new luaSpline(*this);
}


luaSpline::
~luaSpline() {}

const char luaSpline::className[] = "Spline";

Lunar<luaSpline>::RegType luaSpline::member_data[] = {
    member_data(luaSpline, t0),
    member_data(luaSpline, t1),
    {0, 0}
};

Lunar<luaSpline>::RegType luaSpline::methods[] = {
    method(luaSpline, eval),
    method(luaSpline, locate),
    method(luaSpline, dpdt),
    method(luaSpline, length),
    method(luaSpline, partial_length),
    method(luaSpline, point_from_length),
    method(luaSpline, translate),
    method(luaSpline, reverse),
    method(luaSpline, mirror_image),
    {0, 0}
};

Lunar<luaSpline>::MetaType luaSpline::metamethods[] = {
    {0, 0}
};


luaXPath::
luaXPath()
    : luaPath() {}

luaXPath::
~luaXPath() {}

int
luaXPath::
x0(lua_State *L)
{
    int narg = lua_gettop(L);
    // Cast path as an XPath.
    XPath* xpath = xpath_pointer();

    if( narg == 0 ) {
	// This is a getter.
	lua_pushnumber(L, xpath->x0);
	return 1;
    }
    // else
    // Treat as a setter.
    xpath->x0 = luaL_checknumber(L, 1);
    return 0;
}

int
luaXPath::
x1(lua_State *L)
{
    int narg = lua_gettop(L);
    // Cast path as an XPath
    XPath* xpath = xpath_pointer();

    if( narg == 0 ) {
	// This is a getter.
	lua_pushnumber(L, xpath->x1);
	return 1;
    }
    // else
    // Treat as a setter.
    xpath->x1 = luaL_checknumber(L, 1);
    return 0;
}

int
luaXPath::
xeval(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg != 1 ) {
	ostringstream ost;
	ost << "error in call xeval():\n"
	    << "   1 argument expected. " << narg << " arguments received.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
    XPath* xpath = xpath_pointer();
    Vector3 tmp = xpath->xeval(luaL_checknumber(L, 1));
    Lunar<luaVector3>::push(L, new luaVector3(tmp), true);
    return 1; 
}

int
luaXPath::
dydx(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg != 1 ) {
	ostringstream ost;
	ost << "error in call dxdy():\n"
	    << "   1 argument expected. " << narg << " arguments received.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
    XPath* xpath = xpath_pointer();
    lua_pushnumber(L, xpath->dydx(luaL_checknumber(L, 1)));
    return 1;
}

int
luaXPath::
d2ydx2(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg != 1 ) {
	ostringstream ost;
	ost << "error in call d2xdy2():\n"
	    << "   1 argument expected. " << narg << " arguments received.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
    XPath* xpath = xpath_pointer();
    lua_pushnumber(L, xpath->d2ydx2(luaL_checknumber(L, 1)));
    return 1;
}

int
luaXPath::
locate_x(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg != 1 ) {
	ostringstream ost;
	ost << "error in cal locate_x():\n"
	    << "   1 argument expected. " << narg << " arguments received.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
    XPath* xpath = xpath_pointer();
    int result_flag = 0;
    double y = luaL_checknumber(L, 1);
    double x = xpath->locate_x(y, result_flag);
    lua_pushnumber(L, x);
    if( result_flag == SUCCESS )
	lua_pushboolean(L, 1);
    else
	lua_pushboolean(L, 0);
    return 2;
}

string
luaXPoly::
str() const
{
    XPoly *xpoly = dynamic_cast<XPoly*>(path_);
    ostringstream ost;
    ost << "XPoly(" << xpoly->x0 << ", " << xpoly->x1 << ", {";
    for( size_t i =0; i < xpoly->B.size(); ++i ) ost << xpoly->B[i] << ", ";
    ost << "},\n{label=\"" << xpoly->label << "\"})";
    return ost.str();
}

luaXPoly::
luaXPoly(lua_State *L)
    : luaXPath()
{
    // Accept a constructor from Lua of the form:
    // > XPoly(x0, x1, {b1, b2, b3},
    //         {label="label", t0=0.0, t1=1.0})
    //

    double x0, x1;
    handle_x_params_xpath(L, x0, x1);

    vector<double> tmp;
    if( !lua_istable(L, 3) ) {
	cout << "Table of numbers expected.\n";
	exit(1);
    }

    int n = lua_objlen(L, 3);
    for( int i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 3, i);
	tmp.push_back(luaL_checknumber(L, -1));
	lua_pop(L, 1);
    }

    path_ = new XPoly(x0, x1, tmp);

    handle_xpath_defaults(path_, L, 4);

}

luaXPoly::
luaXPoly(const luaXPoly &l)
    : luaXPath()
{
    XPoly* xpoly = dynamic_cast<XPoly*>(l.path_);
    path_ = new XPoly(*xpoly);
}

luaXPoly::
luaXPoly(const XPoly &x)
    : luaXPath()
{
    path_ = new XPoly(x);
}

luaXPoly*
luaXPoly::
clone() const
{
    return new luaXPoly(*this);
}

luaXPoly::
~luaXPoly() {}

int
luaXPoly::
B(lua_State *L)
{
    int narg = lua_gettop(L);

    XPoly *xpoly = dynamic_cast<XPoly*>(path_);

    if( narg == 1 ) {
	// This is a getter.
	int i = luaL_checkint(L, 1);
	lua_pushnumber(L, xpoly->B[i]);
	return 1;
    }
    // else
    // Treat as a setter.
    int i = luaL_checkint(L, 1);
    double val = luaL_checknumber(L, 2);
    xpoly->B[i] = val;
    return 0;
}

int
luaXPoly::
n(lua_State *L)
{
    // Treat as getter only.
    XPoly *xpoly = dynamic_cast<XPoly*>(path_);
    lua_pushinteger(L, xpoly->B.size());
    return 1;
}

const char luaXPoly::className[] = "XPoly";

Lunar<luaXPoly>::RegType luaXPoly::member_data[] = {
    member_data(luaXPoly, t0),
    member_data(luaXPoly, t1),
    member_data(luaXPoly, x0),
    member_data(luaXPoly, x1),
    {0, 0}
};

Lunar<luaXPoly>::RegType luaXPoly::methods[] = {
    method(luaXPoly, B),
    method(luaXPoly, n),
    method(luaXPoly, eval),
    method(luaXPoly, locate),
    method(luaXPoly, dpdt),
    method(luaXPoly, length),
    method(luaXPoly, partial_length),
    method(luaXPoly, point_from_length),
    method(luaXPoly, translate),
    method(luaXPoly, reverse),
    method(luaXPoly, mirror_image),
    // special xpath methods
    method(luaXPoly, xeval),
    method(luaXPoly, dydx),
    method(luaXPoly, d2ydx2),
    method(luaXPoly, locate_x),
    {0, 0}
};

Lunar<luaXPoly>::MetaType luaXPoly::metamethods[] = {
    {0, 0}
};

string
luaXSpline::
str() const
{
    XSpline *xspline = dynamic_cast<XSpline*>(path_);
    ostringstream ost;
    ost << "XSpline(" << xspline->x0 << ", " << xspline->x1 << ", {";
    for( size_t i = 0; i < l_xpoly_.size(); ++i ) ost << endl << l_xpoly_[i]->str() << ", ";
    ost << "},\n{label=\"" << xspline->label << "\"})";
    return ost.str();
}

luaXSpline::
luaXSpline(lua_State *L)
    : luaXPath()
{
    // Accept a constructor from Lua of the form:
    // > XSpline({p1, p2, p3}, {label="label"})

    if( !lua_istable(L, 1) ) {
	cout << "Table of objects expected.\n";
	exit(1);
    }

    int n = lua_objlen(L, 1);
    vector<XPoly> tmp;
    if( n >= 1 ) {
	// Array of XPolys
	for( int i = 1; i <= n; ++i ) {
	    lua_rawgeti(L, 1, i);
	    luaXPoly* lp = Lunar<luaXPoly>::check(L, -1);
	    XPoly xp = XPoly(*(dynamic_cast<XPoly*>(lp->xpath_pointer())));
	    tmp.push_back(xp);
	    l_xpoly_.push_back(lp->clone());
	    lua_pop(L, 1);
	}
    }

    path_ = new XSpline(tmp.front().x0, tmp.back().x1, tmp);
    handle_xpath_defaults(path_, L, 2);
}

luaXSpline::
luaXSpline(const luaXSpline &l)
    : luaXPath()
{
    XSpline* xspline = dynamic_cast<XSpline*>(l.path_);
    path_ = new XSpline(*xspline);
    
    for( size_t i = 0; i < l.l_xpoly_.size(); ++i ) {
	l_xpoly_.push_back(l.l_xpoly_[i]->clone());
    }
}

luaXSpline*
luaXSpline::
clone() const
{
    return new luaXSpline(*this);
}

luaXSpline::
luaXSpline(const XSpline &x)
    : luaXPath()
{
    path_ = new XSpline(x);
    for( size_t i = 0; i < x.X.size(); ++i ) {
	l_xpoly_.push_back(new luaXPoly(x.X[i]));
    }
}

luaXSpline::
~luaXSpline()
{
    for( size_t i = 0; i < l_xpoly_.size(); ++i ) delete l_xpoly_[i];
}

int
luaXSpline::
locate_dydx(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg != 1 ) {
	ostringstream ost;
	ost << "error in call locate_dydx():\n"
	    << "   1 argument expected. " << narg << " arguments received.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
    
    XSpline* xspline = dynamic_cast<XSpline*>(path_);
    int result_flag;
    Vector3 v = xspline->locate_dydx(luaL_checknumber(L, 1), result_flag);
    Lunar<luaVector3>::push(L, new luaVector3(v), true);
    lua_pushboolean(L, result_flag);
    return 2; 
}

const char luaXSpline::className[] = "XSpline";

Lunar<luaXSpline>::RegType luaXSpline::member_data[] = {
    member_data(luaXSpline, t0),
    member_data(luaXSpline, t1),
    member_data(luaXSpline, x0),
    member_data(luaXSpline, x1),
    {0, 0}
};

Lunar<luaXSpline>::RegType luaXSpline::methods[] = {
    method(luaXSpline, eval),
    method(luaXSpline, locate),
    method(luaXSpline, dpdt),
    method(luaXSpline, length),
    method(luaXSpline, partial_length),
    method(luaXSpline, point_from_length),
    method(luaXSpline, translate),
    method(luaXSpline, reverse),
    method(luaXSpline, mirror_image),
    // special xpath methods
    method(luaXSpline, xeval),
    method(luaXSpline, dydx),
    method(luaXSpline, d2ydx2),
    method(luaXSpline, locate_x),
    // special xspline methods
    method(luaXSpline, locate_dydx),
    {0, 0}
};

Lunar<luaXSpline>::MetaType luaXSpline::metamethods[] = {
    {0, 0}
};


string
luaXBezier::
str() const
{
    XBezier *xbez = dynamic_cast<XBezier*>(path_);
    ostringstream ost;
    ost << setprecision(12);
    ost << "XBezier({";
    for( size_t i =0; i < xbez->bez.B.size(); ++i ) ost << xbez->bez.B[i].str(12) << ",\n";
    ost << "},\n{label=\"" << xbez->label << "\"})";
    return ost.str();
}

luaXBezier::
luaXBezier(lua_State *L)
    : luaXPath()
{ 
    // Accept a constructor from Lua of the form:
    // > XBezier({b1, b2, b3},
    //           {label="label"})
    //

    vector<Vector3> tmp;
    if( !lua_istable(L, 1) ) {
	cout << "Table of Vector3 objects expected.\n";
	exit(1);
    }

    int n = lua_objlen(L, 1);
    for( int i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 1, i);
	luaVector3* v = Lunar<luaVector3>::check(L, -1);
	tmp.push_back(v->r_object());
	lua_pop(L, 1);
    }

    path_ = new XBezier(tmp.front().x, tmp.back().x, tmp);

    handle_path_defaults(path_, L, 2);
}

luaXBezier::
luaXBezier(const luaXBezier &l)
    : luaXPath()
{
    XBezier* xbez = dynamic_cast<XBezier*>(l.path_);
    path_ = new XBezier(*xbez);
}

luaXBezier::
luaXBezier(const XBezier &x)
    : luaXPath()
{
    path_ = new XBezier(x);
}

luaXBezier*
luaXBezier::
clone() const
{
    return new luaXBezier(*this);
}

luaXBezier::
~luaXBezier() {}

const char luaXBezier::className[] = "XBezier";

Lunar<luaXBezier>::RegType luaXBezier::member_data[] = {
    member_data(luaXBezier, t0),
    member_data(luaXBezier, t1),
    member_data(luaXBezier, x0),
    member_data(luaXBezier, x1),
    {0, 0}
};

Lunar<luaXBezier>::RegType luaXBezier::methods[] = {
    method(luaXBezier, eval),
    method(luaXBezier, locate),
    method(luaXBezier, dpdt),
    method(luaXBezier, length),
    method(luaXBezier, partial_length),
    method(luaXBezier, point_from_length),
    method(luaXBezier, translate),
    method(luaXBezier, reverse),
    method(luaXBezier, mirror_image),
    // special xpath methods
    method(luaXBezier, xeval),
    method(luaXBezier, dydx),
    method(luaXBezier, d2ydx2),
    method(luaXBezier, locate_x),
    {0, 0}
};

Lunar<luaXBezier>::MetaType luaXBezier::metamethods[] = {
    {0, 0}
};


string
luaXPolyline::
str() const
{
    XPolyline *pline = dynamic_cast<XPolyline*>(path_);
    ostringstream ost;
    ost << "XPolyline({";
    for( size_t i = 0; i < l_xpath_.size(); ++i ) ost << endl << l_xpath_[i]->str() << ", ";
    ost << "},\n{label=\"" << pline->label << "\"})";
    return ost.str();
}

luaXPolyline::
luaXPolyline(lua_State *L)
    : luaXPath()
{
    // Accept a constructor from Lua of the form:
    // > Polyline({path1, path2}, {label="label"})

    if( !lua_istable(L, 1) ) {
	cout << "Table of objects expected.\n";
	exit(1);
    }

    int n = lua_objlen(L, 1);
    vector<XPath*> tmp;
    if( n >= 1 ) {
	// Array of paths
	for( int i = 1; i <= n; ++i ) {
	    lua_rawgeti(L, 1, i);
	    XPath* p = check_xpath(L, -1);
	    luaXPath* lp = new_l_xpath(L, -1);
	    if( p != 0 )
		tmp.push_back(p);
	    if( lp != 0 )
		l_xpath_.push_back(lp);
	    lua_pop(L, 1);
	}
    }

    path_ = new XPolyline(tmp);
    handle_xpath_defaults(path_, L, 2);
}

luaXPolyline::
luaXPolyline(const luaXPolyline &l)
    : luaXPath()
{
    XPolyline* pline = dynamic_cast<XPolyline*>(l.path_);
    path_ = new XPolyline(*pline);
    
    for( size_t i = 0; i < l.l_xpath_.size(); ++i ) {
    	l_xpath_.push_back(l.l_xpath_[i]->clone());
    }
}

luaXPolyline*
luaXPolyline::
clone() const
{
    return new luaXPolyline(*this);
}

luaXPolyline::
~luaXPolyline()
{
    for( size_t i = 0; i < l_xpath_.size(); ++i ) delete l_xpath_[i];
}

const char luaXPolyline::className[] = "XPolyline";

Lunar<luaXPolyline>::RegType luaXPolyline::member_data[] = {
    member_data(luaXPolyline, t0),
    member_data(luaXPolyline, t1),
    member_data(luaXPolyline, x0),
    member_data(luaXPolyline, x1),
    {0, 0}
};

Lunar<luaXPolyline>::RegType luaXPolyline::methods[] = {
    method(luaXPolyline, eval),
    method(luaXPolyline, locate),
    method(luaXPolyline, dpdt),
    method(luaXPolyline, length),
    method(luaXPolyline, partial_length),
    method(luaXPolyline, point_from_length),
    method(luaXPolyline, translate),
    method(luaXPolyline, reverse),
    method(luaXPolyline, mirror_image),
    // special xpath methods
    method(luaXPolyline, xeval),
    method(luaXPolyline, dydx),
    method(luaXPolyline, d2ydx2),
    method(luaXPolyline, locate_x),
    {0, 0}
};

Lunar<luaXPolyline>::MetaType luaXPolyline::metamethods[] = {
    {0, 0}
};


Path* check_path(lua_State *L, int index)
{
    // Just try successive object types
    if( istype(L, index, "Line") ) {
	luaLine* l = Lunar<luaLine>::check(L, index);
	return l->r_pointer();
    }

    if( istype(L, index, "Arc") ) {
	luaArc* a = Lunar<luaArc>::check(L, index);
	return a->r_pointer();
    }

    if( istype(L, index, "Arc3") ) {
	luaArc3* a3 = Lunar<luaArc3>::check(L, index);
	return a3->r_pointer();
    }

    if( istype(L, index, "Bezier") ) {
	luaBezier* b = Lunar<luaBezier>::check(L, index);
	return b->r_pointer();
    }

    if( istype(L, index, "Polyline") ) {
	luaPolyline* p = Lunar<luaPolyline>::check(L, index);
	return p->r_pointer();
    }

    // Nothing found
    return 0;
}

luaPath* new_l_path(lua_State *L, int index)
{
    // Just try successive object types
    if( istype(L, index, "Line") ) {
	luaPath* p = dynamic_cast<luaPath*>(Lunar<luaLine>::check(L, index));
	return p->clone();
    }

    if( istype(L, index, "Arc") ) {
	luaPath* p = dynamic_cast<luaPath*>(Lunar<luaArc>::check(L, index));
	return p->clone();
    }

    if( istype(L, index, "Arc3") ) {
	luaPath* p = dynamic_cast<luaPath*>(Lunar<luaArc3>::check(L, index));
	return p->clone();
    }

    if( istype(L, index, "Bezier") ) {
	luaPath* p = dynamic_cast<luaPath*>(Lunar<luaBezier>::check(L, index));
	return p->clone();
    }

    if( istype(L, index, "Polyline") ) {
	luaPath* p = dynamic_cast<luaPath*>(Lunar<luaPolyline>::check(L, index));
	return p->clone();
    }

    // Nothing found
    return 0;

}

XPath* check_xpath(lua_State *L, int index)
{
    // Just try successive object types
    if( istype(L, index, "XPoly") ) {
	luaXPoly* x = Lunar<luaXPoly>::check(L, index);
	return x->xpath_pointer();
    }

    if( istype(L, index, "XSpline") ) {
	luaXSpline* x = Lunar<luaXSpline>::check(L, index);
	return x->xpath_pointer();
    }

    if( istype(L, index, "XBezier") ) {
	luaXBezier* x = Lunar<luaXBezier>::check(L, index);
	return x->xpath_pointer();
    }

    if( istype(L, index, "XPolyline") ) {
	luaXPolyline* x = Lunar<luaXPolyline>::check(L, index);
	return x->xpath_pointer();
    }
    // Nothing found
    return 0;
}

luaXPath* new_l_xpath(lua_State *L, int index)
{
    // Just try successive object types
    if( istype(L, index, "XPoly") ) {
	luaXPath* p = dynamic_cast<luaXPath*>(Lunar<luaXPoly>::check(L, index));
	return p->clone();
    }

    if( istype(L, index, "XSpline") ) {
	luaXPath* p = dynamic_cast<luaXPath*>(Lunar<luaXSpline>::check(L, index));
	return p->clone();
    }

    if( istype(L, index, "XBezier") ) {
	luaXPath* p = dynamic_cast<luaXPath*>(Lunar<luaXBezier>::check(L, index));
	return p->clone();
    }

    if( istype(L, index, "XPolyline") ) {
	luaXPath* p = dynamic_cast<luaXPath*>(Lunar<luaXPolyline>::check(L, index));
	return p->clone();
    }
    // Nothing found
    return 0;

}

void handle_path_defaults(Path *p, lua_State *L, int index)
{
    if( lua_istable(L, index) ) {
	// Extract optional arguments
	lua_getfield(L, index, "label");
	if( lua_isstring(L, -1) ) {
	    p->label = string(lua_tostring(L, -1));
	}
	else {
	    // No sensible value for lable
	    p->label = "(default label)";
	}
	lua_pop(L, 1);

	lua_getfield(L, index, "t0");
	if( lua_isnumber(L, -1) ) {
	    p->t0 = lua_tonumber(L, -1);
	}
	else {
	    // No sensible value for t0
	    // Assume default 0.0
	    p->t0 = 0.0;
	}
	lua_pop(L, 1);

	lua_getfield(L, index, "t1");
	if( lua_isnumber(L, -1) ) {
	    p->t1 = lua_tonumber(L, -1);
	}
	else {
	    // No sensible value for t1
	    // Assume default 1.0
	    p->t1 = 1.0;
	}
	lua_pop(L, 1);
    }
    else {
	// Not a table passed in.
	// Just set some defaults.
	p->label = "(default label)";
	p->t0 = 0.0;
	p->t1 = 1.0;
    }
}

void handle_x_params_xpath(lua_State *L, double &x0, double &x1)
{
    x0 = luaL_checknumber(L, 1);
    x1 = luaL_checknumber(L, 2);
}

void handle_xpath_defaults(Path *p, lua_State *L, int index)
{

    if( lua_istable(L, index) ) {
	// Extract optional arguments
	lua_getfield(L, index, "label");
	if( lua_isstring(L, -1) ) {
	    string label(string(lua_tostring(L, -1)));
	    p->label = label;
	}
	else {
	    // No sensible value for label
	    p->label = "(default label)";
	}
	lua_pop(L, 1);
    }
    else {
	// Not a table passed in.
	// Just set some defaults.
	p->label = "(default label)";
    }
}

int open_gpath(lua_State *L, int table)
{
    Lunar<luaLine>::Register(L, table);
    Lunar<luaArc>::Register(L, table);
    Lunar<luaArc3>::Register(L, table);
    Lunar<luaBezier>::Register(L, table);
    Lunar<luaNurbs>::Register(L, table);
    Lunar<luaPolyline>::Register(L, table);
    Lunar<luaSpline>::Register(L, table);
    Lunar<luaXPoly>::Register(L, table);
    Lunar<luaXSpline>::Register(L, table);
    Lunar<luaXBezier>::Register(L, table);
    Lunar<luaXPolyline>::Register(L, table);
    return 0;
}


