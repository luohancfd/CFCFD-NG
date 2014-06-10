# README.txt for mms

This is a work in progress, starting with Rowan's mms_euler case.
The preparation of the source terms is now more automated so that
we have a lot less to maintain.  The trade is that we now need maxima
installed to run the cases.

Remember to set the case number.
It is a single digit on the only line in the file "case.txt".
1 == supersonic inviscid flow
2 == subsonic viscous flow

When checking the supersonic inviscid case (1), we have been looking at 
the final time and comparing the density field with the analytic solution, 
as plotted in Rowan's thesis.  
The key value is an L2 norm of 0.00476 for the grid with ni=nj=16.

The corresponding value for the viscous case (2) is 0.003648.

----------------------------------------------------------------------------
30-June-2012
Have scaled the perturbation to reduce it away from the centre of the domain.
This is done via the "S" expressions and functions in the Maxima, Python 
and Lua files.
The residual for case 1 is now 
Volume-weighted norms, global:
rho L2 0.00151061874308

----------------------------------------------------------------------------
03-July-2012

A lot of changes added:
 - case.txt now controls a variety of parameters
 - in the code base, 1st order boundary conditions are implemented
 - the scaled perturbation cases are now known as '3' and '4'

With first order boundary conditions in place, the error for
the L2 norm in density on a ni=nj=16 grid has dropped considerably to:
  0.001086
----------------------------------------------------------------------------
10-June-2014

In this change, we now use sympy to generate the source terms. This has introduced
a small difference in the compute L2 norm in density for case 2 (and probably
for the other cases as well). For some time, we have been running case 2 as the
default in the test suite. I have convinced myself that this difference in L2
norm is because of subtle differences in the precision of floating point value representations
when using maxima-generated source terms as compared to sympy-generated source terms.
For the record, I have included some images that compare these differences.
------------------------------------------------------------------------------



