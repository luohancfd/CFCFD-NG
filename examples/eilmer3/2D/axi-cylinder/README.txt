Update for 14-Apr-2013
Try out the high-order time stepping to see if we can make things go faster.
Set the gasdynamic_update_scheme to classic_rk3 and cfl to 1.2
Run with the MPI code using 4 procs on the dual-core Toshiba C850.
It now takes 2607 seconds (compared with 7235 seconds) to take 66239 steps
(cf 158935 steps) for 2D/cyl50/ case.
