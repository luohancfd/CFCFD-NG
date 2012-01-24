/// \file pypath.cxx
/// \ingroup libgeom2
/// \brief Python-powered C++ geometric-path class.
/// \author PJ

#include "pypath.hh"

#include <iostream>
#include <string>
#include <sstream>

PyFunctionPath::PyFunctionPath(PyObject* pyfun, const string label, double t0, double t1)
    : Path(label, t0, t1)
{
    if ( !PyCallable_Check(pyfun) ) {
	PyErr_SetString(PyExc_TypeError, "Need a callable object.");
	pyfunc = NULL;
    } else {
	pyfunc = pyfun;
	Py_INCREF(pyfunc);
    }
}
PyFunctionPath::PyFunctionPath(const PyFunctionPath &pypath)
    : Path(pypath.label, pypath.t0, pypath.t1) 
{
    pyfunc = pypath.pyfunc;
    Py_INCREF(pyfunc);
}
PyFunctionPath::~PyFunctionPath()
{
    Py_XDECREF(pyfunc);
    pyfunc = NULL;
}
PyFunctionPath* PyFunctionPath::clone() const
{
    return new PyFunctionPath(*this);
} 
Vector3 PyFunctionPath::eval(double t) const
{
    Vector3 p;
    t = t0 + (t1 - t0) * t;
    PyObject *arglist = Py_BuildValue("(d)", t); // tuple needed, even if only one arg
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
string PyFunctionPath::str() const
{
    ostringstream ost;
    ost << "PyFunctionPath(" << "<PythonFunction>" << ", ";
    ost << "\"" << label << "\"" << ", ";
    ost << t0 << ", " << t1 << ")";
    return ost.str();
}
