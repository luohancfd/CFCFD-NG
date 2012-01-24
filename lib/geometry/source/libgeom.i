/** \file libgeom.i
 * \ingroup geom 
 * \brief SWIG interface header file for geometry library.
 */

%module libgeom 
%{
#include <stdio.h>
#include "../source/geom.h"
#include "../source/gpath.h"
#include "../source/bezier.h"
%}

%typemap(python, in) FILE * {
    if (!PyFile_Check($input)) {
        PyErr_SetString(PyExc_TypeError, "expected PyFile");
        return NULL;
    }
    $1=PyFile_AsFile($input);
}

%include "../source/geom.h"
%include "../source/gpath.h"
%include "../source/bezier.h"
