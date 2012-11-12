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
04-Nov-2012

The first crack at 3D.  Let's see what's broken...
