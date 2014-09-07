# makefile for the gas module

unittest: perfectgasEOS calperfectgasEOS perfectgasmixEOS calperfectgasmixEOS

perfectgasEOS : perfectgasEOS.d thermalEOS.d ../gasmodel.d
	dmd -main -unittest perfectgasEOS.d thermalEOS.d ../gasmodel.d -I../

calperfectgasEOS : calperfectgasEOS.d caloricEOS.d ../gasmodel.d
	dmd -main -unittest calperfectgasEOS.d caloricEOS.d ../gasmodel.d -I../

perfectgasmixEOS : perfectgasmixEOS.d perfectgasEOS.d thermalEOS.d ../gasmodel.d
	dmd -main -unittest perfectgasmixEOS.d perfectgasEOS.d thermalEOS.d ../gasmodel.d -I../

calperfectgasmixEOS : calperfectgasmixEOS.d calperfectgasEOS.d caloricEOS.d ../gasmodel.d
	dmd -main -unittest calperfectgasmixEOS.d calperfectgasEOS.d caloricEOS.d ../gasmodel.d -I../

