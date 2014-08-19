# makefile for the diffusion module

unittest: sutherland_viscosity

sutherland_viscosity : sutherland_viscosity viscosity.d ../gasmodel.d
	dmd -main -unittest sutherland_viscosity.d viscosity.d ../gasmodel.d -I../


