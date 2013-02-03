WITH_SPRADIAN ?= 0
WITH_PARADE ?= 0

OPTIONS := --with-photaura --with-equilibrium-air
ifeq ($(WITH_SPRADIAN), 1)
    OPTIONS += --with-spradian
endif
ifeq ($(WITH_PARADE), 1)
    OPTIONS += --with-parade
endif

test: 
	- ./test-spectral-models.py  $(OPTIONS) 

clean: 
	- rm *.txt gas-model.lua rad-model.lua
	- rm -R -f spradian_working_dir_0 parade_working_dir_0
	- rm thermo.* trans.* tmp.*
