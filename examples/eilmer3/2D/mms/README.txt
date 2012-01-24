# README.txt for mms

This is a work in progress, starting with Rowan's mms_euler case.
The preparation of the source terms is now more automated so that
we have a lot less to maintain.  The trade is that we now need maxima
installed to run the cases.

When checking the supersonic inviscid case, we have been looking at 
the final time and comparing the density field with the analytic solution, 
as plotted in Rowan's thesis.  
The key value is an L2 norm of 0.00476 for the grid with ni=nj=16.

