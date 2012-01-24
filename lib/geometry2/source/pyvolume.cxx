/// \file pyvolume.cxx
/// \ingroup libgeom2
/// \brief Python-powered C++ parametric-volume class.
/// \author PJ

#include "pyvolume.hh"
#include "surface.hh"
#include "volume.hh"

#include <iostream>
#include <string>
#include <sstream>

PyFunctionVolume::PyFunctionVolume( PyObject* pyfun, string label,
				    double r0, double r1, 
				    double s0, double s1,
				    double t0, double t1 )
    : ParametricVolume(label, r0, r1, s0, s1, t0, t1)
{
    if ( !PyCallable_Check(pyfun) ) {
	PyErr_SetString(PyExc_TypeError, "Need a callable object.");
	pyfunc = NULL;
    } else {
	pyfunc = pyfun;
	Py_INCREF(pyfunc);
    }
    set_surfaces();
    set_and_check_corners();
}
PyFunctionVolume::PyFunctionVolume( const PyFunctionVolume &pv )
    : ParametricVolume(pv.label, pv.r0, pv.r1, pv.s0, pv.s1, pv.t0, pv.t1) 
{
    pyfunc = pv.pyfunc;
    Py_INCREF(pyfunc);
    set_surfaces();
    set_and_check_corners();
}
PyFunctionVolume::~PyFunctionVolume()
{
    Py_XDECREF(pyfunc);
    pyfunc = NULL;
    // Base class deletes the surfaces.
}
PyFunctionVolume* PyFunctionVolume::clone() const
{
    return new PyFunctionVolume(*this);
} 
Vector3 PyFunctionVolume::eval(double r, double s, double t) const
{
    Vector3 p;
    r = r0 + (r1 - r0) * r;
    s = s0 + (s1 - s0) * s;
    t = t0 + (t1 - t0) * t;
    PyObject *arglist = Py_BuildValue("(ddd)", r, s, t);
    PyObject *result = PyEval_CallObject(pyfunc, arglist);
    Py_DECREF(arglist);
    if ( result ) {
	double x, y, z;
	PyArg_ParseTuple(result, "ddd", &x, &y, &z);
	p = Vector3(x, y, z);
    } else {
	p = Vector3(0.0, 0.0, 0.0);
    }
    Py_XDECREF(result);
    return p;
}
string PyFunctionVolume::str() const
{
    ostringstream ost;
    ost << "PyFunctionVolume(" << "<PythonFunction>" << ", ";
    ost << "\"" << label << "\"" << ", ";
    ost << r0 << ", " << r1 << ", " 
	<< s0 << ", " << s1 << ", " 
	<< t0 << ", " << t1 << ")";
    return ost.str();
}
int PyFunctionVolume::set_surfaces()
{
    cout << "PyFunctionVolume::set_surfaces() using MeshPatch objects" << endl;
    int ni = 10;
    int nj = 10;
    int nk = 10;
    int i, j, k;
    double dr = 1.0/(ni-1);
    double ds = 1.0/(nj-1);
    double dt = 1.0/(nk-1);
    double r, s, t;
    vector<Vector3> south_plist;
    vector<Vector3> north_plist;
    for ( k = 0; k < nk; ++k ) {
	t = k * dt;
	for ( i = 0; i < ni; ++i ) {
	    r = i * dr;
	    s = 0.0;
	    south_plist.push_back(eval(r,s,t));
	    s = 1.0;
	    north_plist.push_back(eval(r,s,t));
	}
    }
    south = new MeshPatch(south_plist, ni, nk);
    north = new MeshPatch(north_plist, ni, nk);
    vector<Vector3> west_plist;
    vector<Vector3> east_plist;
    for ( k = 0; k < nk; ++k ) {
	t = k * dt;
	for ( j = 0; j < nj; ++j ) {
	    s = j * ds;
	    r = 0.0;
	    west_plist.push_back(eval(r,s,t));
	    r = 1.0;
	    east_plist.push_back(eval(r,s,t));
	}
    }
    west = new MeshPatch(west_plist, nj, nk);
    east = new MeshPatch(east_plist, nj, nk);
    vector<Vector3> bottom_plist;
    vector<Vector3> top_plist;
    for ( j = 0; j < nj; ++j ) {
	s = j * ds;
	for ( i = 0; i < ni; ++i ) {
	    r = i * dr;
	    t = 0.0;
	    bottom_plist.push_back(eval(r,s,t));
	    t = 1.0;
	    top_plist.push_back(eval(r,s,t));
	}
    }
    bottom = new MeshPatch(bottom_plist, ni, nj);
    top = new MeshPatch(top_plist, ni, nj);
    return 0;
}
