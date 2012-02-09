/* 
 
 PyLevmar - Python Bindings for Levmar. GPL2.
 Copyright (c) 2006. Alastair Tse <alastair@liquidx.net>
 Copyright (c) 2008. Dustin Lang <dstndstn@gmail.com>
 
 Modified to run numpy ... (c) 2010 Stuart Wilkins <stuwilkins@mac.com>
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 
 $Id: pylevmar.h 42 2010-10-13 23:23:42Z stuwilkins $
 
 */

#ifndef __PYLEVMAR_H
#define __PYLEVMAR_H

#include "lm.h"

#include <Python.h>
#include <numpy/arrayobject.h>
#include "structmember.h"

// Type definitions

typedef struct _pylm_callback_data {
    PyObject *func;
    PyObject *jacf;
} pylm_callback_data;

// Function Prototypes

static PyObject * _pylm_dlevmar_generic(PyObject *mod, PyObject *args, PyObject *kwds,
										char *argstring, char *kwlist[],
										int jacobian, int	bounds);

static PyObject * pylm_dlevmar_der(PyObject *mod, PyObject *args, PyObject *kwds);
static PyObject * pylm_dlevmar_dif(PyObject *mod, PyObject *args, PyObject *kwds);
static PyObject * pylm_dlevmar_bc_der(PyObject *mod, PyObject *args, PyObject *kwds);
static PyObject * pylm_dlevmar_bc_dif(PyObject *mod, PyObject *args, PyObject *kwds);
static PyObject * pylm_dlevmar_chkjac(PyObject *mod, PyObject *args, PyObject *kwds);

void _pylm_callback(PyObject *func, double *p, double *hx, int m, int n, int  jacobian);
void _pylm_func_callback(double *p, double *hx, int m, int n, void *data);
void _pylm_jacf_callback(double *p, double *hx, int m, int n, void *data);


// Python Module Documents

static char *pylm_doc = \
"Python Bindings to Levmar Non-Linear Regression Solver.\n\n" 
"Members:\n\n"
"  levmar.DEFAULT_OPTS = 4 element tuple representing the default\n"
"                        opts into dder an ddif.\n"
"  levmar.INIT_MU      = Initial value for mu\n"
"  levmar.STOP_THRESH  = Stopping threshold\n"
"  levmar.DIFF_DELTA   = Differential delta\n";


static PyMethodDef pylm_functions[] = {
    {
        "dder",
        (PyCFunction)pylm_dlevmar_der,
        METH_VARARGS|METH_KEYWORDS,
        "dlevmar_der(func, jacf, estimate, measurements, itmax,\n"
        "            opts = None, covar = None)\n"
        "-> returns: (result, iterations, run_info)\n"
        "\n"
    },
    {
        "ddif",
        (PyCFunction)pylm_dlevmar_dif,
        METH_VARARGS|METH_KEYWORDS,
        "dlevmar_der(func, estimate, measurements, itmax,\n"
        "            opts = None, covar = None)\n"
        "-> returns: (result, iterations, run_info)\n"
        "\n"
    },
    {
        "dchkjac",
        (PyCFunction)pylm_dlevmar_chkjac,
        METH_VARARGS|METH_KEYWORDS,
        "dlevmar_der(func, jacf, initial, m, data = None) -> tuple\n\n"
        "Check a function and the jacobian and the relative error at point\n\n"
        "initial."
    },
    {
        "ddif_bc",
        (PyCFunction)pylm_dlevmar_bc_dif,
        METH_VARARGS|METH_KEYWORDS,
        "dlevmar_bc_dif(func, estimate, measurements, lower, upper, itmax,\n"
        "               opts = None, covar = None)\n"
        "-> returns: (result, iterations, run_info)\n"
        "\n"
    },
    {
        "dder_bc",
        (PyCFunction)pylm_dlevmar_bc_der,
        METH_VARARGS|METH_KEYWORDS,
        "dlevmar_bc_der(func, estimate, measurements, lower, upper, itmax,\n"
        "               opts = None, covar = None)\n"
        "-> returns: (result, iterations, run_info)\n"
        "\n"
    },
    {NULL, NULL, 0, NULL}
};

#endif
