# makefile for the gas module

DEMO_PROGRAMS := gasmodel_demo idealgas_demo

demo: ${DEMO_PROGRAMS}
	echo "Demo programs built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

gasmodel_demo: gasmodel_demo.d gasmodel.d idealgas.d gasmodelutil.d \
	thermo/thermalEOS.d thermo/perfectgasEOS.d thermo/caloricEOS.d thermo/calperfectgasEOS.d \
	diffusion/viscosity.d diffusion/sutherland_visc.d diffusion/therm_cond.d diffusion/sutherland_therm_cond.d
	dmd -unittest -g gasmodel_demo.d gasmodel.d gasmodelutil.d idealgas.d thermo/thermalEOS.d thermo/perfectgasEOS.d thermo/caloricEOS.d thermo/calperfectgasEOS.d diffusion/viscosity.d diffusion/sutherland_visc.d diffusion/therm_cond.d diffusion/sutherland_therm_cond.d

idealgas_demo: idealgas_demo.d gasmodel.d idealgas.d \
	thermo/thermalEOS.d thermo/perfectgasEOS.d thermo/caloricEOS.d thermo/calperfectgasEOS.d \
	diffusion/viscosity.d diffusion/sutherland_visc.d diffusion/therm_cond.d diffusion/sutherland_therm_cond.d
	dmd -unittest -g idealgas_demo.d gasmodel.d idealgas.d thermo/thermalEOS.d thermo/perfectgasEOS.d thermo/caloricEOS.d thermo/calperfectgasEOS.d diffusion/viscosity.d diffusion/sutherland_visc.d diffusion/therm_cond.d diffusion/sutherland_therm_cond.d
