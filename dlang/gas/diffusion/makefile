# makefile for the diffusion module

unittest: sutherland_visc sutherland_therm_cond

sutherland_visc : sutherland_visc.d viscosity.d ../gasmodel.d
	dmd -main -unittest sutherland_visc.d viscosity.d ../gasmodel.d -I../

sutherland_therm_cond : sutherland_therm_cond.d therm_cond.d ../gasmodel.d
	dmd -main -unittest sutherland_therm_cond.d therm_cond.d ../gasmodel.d -I../


