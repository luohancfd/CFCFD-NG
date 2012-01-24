## \file flux_dict.py
## \ingroup flux_calc
##
## \author R.Gollan and PJ
## \version 31-Jan-2005
##
"""
Dictionary to look up flux-calculator index from name or number.

@var RIEMANN: An exact Riemann-solver-based flux calculator.
    It is slow and only works for ideal gas models at the moment.
@var AUSM: M. S. Liou's AUSM approximate flux calculator.
    Fast but tends to be a bit noisy.
@var EFM: Mike Macrossan's version of Dale Pullin's equilibrium
    flux calculator as coded by Paul Petrie-Repar.
    When you need a dissipative scheme, this is a good one.
@var AUSMDV: A version of Wada and Liou's AUSMDV scheme.
@var ADAPTIVE: A switched AUSMDV/EFM scheme that uses EFM near
    shocks and AUSMDV elsewhere.
    This is a good all-rounder for shock-tunnel work.
@var AUSM_PLUS_UP: An updated version of Liou's AUSM flux
    calculator that includes (1) preconditioning to make it 
    suitable for calculations at all speeds and (2) a 
    pressure-velocity coupling modification to remove odd-even 
    decoupling issues encountered in low speed flows.
    
@undocumented: fluxcalcIndexFromName
"""

RIEMANN  = 0
AUSM     = 1
EFM      = 2
AUSMDV   = 3
ADAPTIVE = 4
AUSM_PLUS_UP = 5

fluxcalcIndexFromName = {
    0: RIEMANN, "0": RIEMANN, "riemann": RIEMANN, "RIEMANN": RIEMANN, "Riemann": RIEMANN,
    1: AUSM, "1": AUSM, "ausm": AUSM, "AUSM": AUSM,
    2: EFM, "2": EFM, "efm": EFM,  "EFM": EFM,
    3: AUSMDV, "3": AUSMDV, "ausmdv": AUSMDV, "AUSMDV": AUSMDV,
    4: ADAPTIVE, "4": ADAPTIVE, "adaptive": ADAPTIVE, "ADAPTIVE": ADAPTIVE,
    5: AUSM_PLUS_UP, "5": AUSM_PLUS_UP, "ausm_plus_up": AUSM_PLUS_UP, "AUSM_PLUS_UP": AUSM_PLUS_UP
}

