/// \file pysurface.cxx
/// \ingroup libgeom2
/// \brief Python-powered C++ geometric-surface class.
/// \author PJ

#include "pysurface.hh"

#include <iostream>
#include <string>
#include <sstream>

PyFunctionSurface::PyFunctionSurface( PyObject* pyfun, string label,
				      double r0, double r1, double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1)
{
    if ( !PyCallable_Check(pyfun) ) {
	PyErr_SetString(PyExc_TypeError, "Need a callable object.");
	pyfunc = NULL;
    } else {
	pyfunc = pyfun;
	Py_INCREF(pyfunc);
    }
}
PyFunctionSurface::PyFunctionSurface( const PyFunctionSurface &surf )
    : ParametricSurface(surf.label, surf.r0, surf.r1, surf.s0, surf.s1) 
{
    pyfunc = surf.pyfunc;
    Py_INCREF(pyfunc);
}
PyFunctionSurface::~PyFunctionSurface()
{
    Py_XDECREF(pyfunc);
    pyfunc = NULL;
}
PyFunctionSurface* PyFunctionSurface::clone() const
{
    return new PyFunctionSurface(*this);
} 
Vector3 PyFunctionSurface::eval( double r, double s ) const
{
    Vector3 p;
    r = r0 + (r1 - r0) * r;
    s = s0 + (s1 - s0) * s;
    PyObject *arglist = Py_BuildValue("(dd)", r, s);
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
string PyFunctionSurface::str() const
{
    ostringstream ost;
    ost << "PyFunctionSurface(" << "<PythonFunction>" << ", ";
    ost << "\"" << label << "\"" << ", ";
    ost << r0 << ", " << r1 << ", " << s0 << ", " << s1 << ")";
    return ost.str();
}
