/// \file pyvolume.hh
/// \ingroup libgeom2
/// \brief Declarations for the Python-powered C++ parametric class.
/// \author PJ

#ifndef PYVOLUME_HH
#define PYVOLUME_HH

#ifndef SWIG
#ifndef PYTHON_H_ALREADY_INCLUDED
#define PYTHON_H_ALREADY_INCLUDED
extern "C" {
#include <Python.h>
}
#endif
#include "volume.hh"
#endif

/** \brief A volume defined by a Python function.
 *
 */
class PyFunctionVolume : public ParametricVolume {
public:
    PyObject *pyfunc;  // Python function that is called by eval().
    PyFunctionVolume(PyObject* pyfunc, string label="",
		     double r0 = 0.0, double r1 = 1.0,
		     double s0 = 0.0, double s1 = 1.0,
		     double t0 = 0.0, double t1 = 1.0);
    /// Construct as a copy of another PyFunctionSurface surface.
    PyFunctionVolume(const PyFunctionVolume &pv);
    virtual ~PyFunctionVolume();
    virtual PyFunctionVolume* clone() const; 
    virtual Vector3 eval(double r, double s, double t) const;
    virtual string str() const;
protected:
    int set_surfaces();
};

#endif
