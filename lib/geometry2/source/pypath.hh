/// \file pypath.hh
/// \ingroup libgeom2
/// \brief Declarations for the Python-powered C++ geometric-path class.
/// \author PJ

#ifndef PYPATH_HH
#define PYPATH_HH

#ifndef SWIG
#ifndef PYTHON_H_ALREADY_INCLUDED
#define PYTHON_H_ALREADY_INCLUDED
extern "C" {
#include <Python.h>
}
#endif
#include "gpath.hh"
#endif

/** \brief A geometric path defined by a Python function.
 *
 */
class PyFunctionPath : public Path {
public:
    PyObject *pyfunc;  // Python function that is called by eval().
    PyFunctionPath(PyObject* pyfun, const string label="", double t0=0.0, double t1=1.0);
    /// Construct as a copy of another PyFunctionPath object.
    PyFunctionPath(const PyFunctionPath &pypath);
    virtual ~PyFunctionPath();
    virtual PyFunctionPath* clone() const; 
    virtual Vector3 eval(double t) const;
    virtual string str() const;
};

#endif
