# makefile for the gas module

unittest: perfectgasEOS

perfectgasEOS : perfectgasEOS.d thermalEOS.d ../gasmodel.d
	dmd -main -unittest perfectgasEOS.d thermalEOS.d ../gasmodel.d -I../


