.. highlight:: python

============================
Fitting of Experimental Data
============================

************
Introduction
************

These routines / classes provide a method for fitting of data using mostly least
squares methods. There are two main methods here. The ``pyspec.fit.fit`` class provides a
class for fitting of data. The ``pyspec.fit.fitdata`` subroutine serves as a wrapper around the ``pyspec.fit.fit`` class. Here the data is taken from the current selected figure. Basically this means that you can fit any set of data which is on the current figure.

*****************
The ``fit`` Class
*****************

Fitting user supplied data
--------------------------

.. autoclass:: pyspec.fit.fit 
   :members:

Fitting data plotted using matplotlib
-------------------------------------

.. autofunction:: pyspec.fit.fitdata

**********************
Standard Fit Functions
**********************

.. automodule:: pyspec.fitfuncs
   :members:

**********************
User Defined Functions
**********************

User defined functions can easily be written. An example of a fit
function to provide a linear (straight line) is shown below::

    def linear(x, p, mode='eval'): 
       if mode == 'eval':
          out = (p[0] * x) + p[1]
       elif mode == 'params':
          out = ['grad','offset']
       elif mode == 'name':
          out = "Linear"
       elif mode == 'guess':
          g = peakguess(x, p)
          out = g[4:6]
       else:
          out = []
       return out

This function is called by the fitting routine each time the function
needs to be evaluated. It is also called other times, to gain
information about the function. The ``mode`` parameter is used to let
the function know what is required. The default value of ``mode``
should be set to ``'eval'``. The possible values of ``mode`` are
discussed below

``eval``
    **Evaluate function**. Here the function should return
    :math:`f(x)` from ``x`` and ``p``, where ``p`` is the parameters
``params``
    **List of parameter names**. Here the
    function should return a list of strings describing the fit
    parameters. This list should have the same length as the number of
    fit parameters.
``name``
    **Name of the function**. Here the function should return a string
    with a text description of the function.
``guess``
    **Parameter guess**. Here the function should return a guess of
    the fit parameters. The function is passed :math:`f(x)` through
    the second parameter ``p``. The function should return it's best
    guess based on the experimental data, and return those
    parameters. The ``peakguess`` function is provided for convenience
    which can guess the *vital statistics* of a peak from the
    data. This option is not essential but if it is not implemented
    then the user must provide a numerical guess to the ``fit`` class.

Finally, the default response should be to return a null list ``[]``. 

