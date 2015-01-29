Update for 14-Apr-2013 (while at Oxford)
Try out the high-order time stepping to see if we can make things go faster.
Set the gasdynamic_update_scheme to classic_rk3 and cfl to 1.2
Run with the MPI code using 4 procs on the dual-core Toshiba C850.
It now takes 2607 seconds (compared with 7235 seconds) to take 66239 steps
(cf 158935 steps) for 2D/cyl50/ case.

29-Jan-2015
A run on Helmholtz (with 4-cores AMD Phenom II X4 840) now takes 1678 seconds
for 74766 steps.  The total CPU time was 6619 seconds, indicating that all 
cores were in use.  The increase in the number of steps may be related to 
code clean-up in the CFL check some time in the past couple of years.
