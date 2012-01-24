#!/bin/bash

cp marrone.rsi air-7sp-10r.rsi
poshax.x air-Mach_12.3.cfg
gnuplot plot_data.gnu
mv profile_T_rho.eps profile_T_rho_marrone.eps
mv profile_moles.eps profile_moles_marrone.eps

cp gupta_etal.rsi air-7sp-10r.rsi
poshax.x air-Mach_12.3.cfg
gnuplot plot_data.gnu

