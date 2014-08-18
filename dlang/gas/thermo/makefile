# makefile for the gas module

unittest: perfectgasEOS calperfectgasEOS

perfectgasEOS : perfectgasEOS.d thermalEOS.d ../gasmodel.d
	dmd -main -unittest perfectgasEOS.d thermalEOS.d ../gasmodel.d -I../

calperfectgasEOS : calperfectgasEOS.d caloricEOS.d ../gasmodel.d
	dmd -main -unittest calperfectgasEOS.d caloricEOS.d ../gasmodel.d -I../



