# makefile for the diffusion module

unittest: sutherland_viscosity sutherland_therm_cond

sutherland_viscosity : sutherland_viscosity.d viscosity.d ../gasmodel.d
	dmd -main -unittest sutherland_viscosity.d viscosity.d ../gasmodel.d -I../

sutherland_therm_cond : sutherland_therm_cond.d therm_cond.d ../gasmodel.d
	dmd -main -unittest sutherland_therm_cond.d therm_cond.d ../gasmodel.d -I../


