#! /bin/sh
# sharp_run.sh
# Exercise the Navier-Stokes solver for a sharp 2D body.

# Integrate the solution in time.
time e3shared.exe --job=sharp --run

echo At this point, we should have a final solution in sharp.b0000.t0015

