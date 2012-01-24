%module geometry

%include "std_vector.i"
%include "std_string.i"

%{
#include "geom.hh"
#include "gpath.hh"
#include "gpath_utils.hh"
%}

// Lua does not provide this kind of operator

%ignore Vector3::operator+=;
%ignore Vector3::operator-=;
%ignore Vector3::operator*=;
%ignore Vector3::operator/=;
%ignoer Vector3::operator=;
%ignore operator<<;
%ignore operator+;
%ignore operator-;
%ignore operator*;
%ignore operator/;
%ignore operator=;

%template(Vector3List) std::vector<Vector3>;
%template(vectord) std::vector<double>;
%template(XPathList) std::vector<XPath*>;
%template(XPolyList) std::vector<XPoly>;

%include "geom.hh"
%include "gpath.hh"
%include "gpath_utils.hh"


