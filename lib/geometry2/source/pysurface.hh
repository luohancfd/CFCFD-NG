/// \file pysurface.hh
/// \ingroup libgeom2
/// \brief Declarations for the Python-powered C++ geometric-surface class.
/// \author PJ

#ifndef PYSURFACE_HH
#define PYSURFACE_HH

#ifndef SWIG
#ifndef PYTHON_H_ALREADY_INCLUDED
#define PYTHON_H_ALREADY_INCLUDED
extern "C" {
#include <Python.h>
}
#endif
#include "surface.hh"
#endif

/** \brief A surface defined by a Python function.
 *
 */
class PyFunctionSurface : public ParametricSurface {
public:
    PyObject *pyfunc;  // Python function that is called by eval().
    PyFunctionSurface( PyObject* pyfunc, string label="",
		     double r0 = 0.0, double r1 = 1.0,
		     double s0 = 0.0, double s1 = 1.0 );
    /// Construct as a copy of another PyFunctionSurface surface.
    PyFunctionSurface( const PyFunctionSurface &surf );
    virtual ~PyFunctionSurface();
    virtual PyFunctionSurface* clone() const; 
    virtual Vector3 eval( double r, double s ) const;
    virtual string str() const;
};

#endif
