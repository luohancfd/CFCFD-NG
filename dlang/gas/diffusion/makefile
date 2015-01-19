# makefile for the diffusion module
UTIL_DIR := ../../util
UTIL_FILES := $(UTIL_DIR)/lua_service.d

LUA := ../../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib
LUAD_DIR := ../../extern/LuaD
LUAD_FILES := $(LUAD_DIR)/luad/*.d \
	$(LUAD_DIR)/luad/c/*.d \
	$(LUAD_DIR)/luad/conversions/*.d

$(LIBLUA):
	cd $(LUA); make linux local

unittest: sutherland_visc sutherland_therm_cond \
	wilke_mixing_visc wilke_mixing_therm_cond \
	cea_viscosity cea_therm_cond
	./sutherland_visc
	./sutherland_therm_cond
	./wilke_mixing_visc
	./wilke_mixing_therm_cond
	./cea_viscosity


sutherland_visc : sutherland_visc.d viscosity.d ../gasmodel.d $(UTIL_FILES) $(LIBLUA)
	dmd -main -unittest sutherland_visc.d viscosity.d ../gasmodel.d -I../ $(UTIL_FILES) ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

sutherland_therm_cond : sutherland_therm_cond.d therm_cond.d ../gasmodel.d $(UTIL_FILES) $(LIBLUA)
	dmd -main -unittest sutherland_therm_cond.d therm_cond.d ../gasmodel.d -I../ $(UTIL_FILES) ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

wilke_mixing_visc : wilke_mixing_visc.d viscosity.d ../gasmodel.d
	dmd -main -unittest wilke_mixing_visc.d sutherland_visc.d viscosity.d ../gasmodel.d -I../

wilke_mixing_therm_cond : wilke_mixing_therm_cond.d therm_cond.d ../gasmodel.d
	dmd -main -unittest wilke_mixing_therm_cond.d sutherland_therm_cond.d therm_cond.d ../gasmodel.d -I../

cea_viscosity : cea_viscosity.d viscosity.d ../gasmodel.d $(UTIL_FILES) $(LIBLUA)
	dmd -main -unittest cea_viscosity.d viscosity.d ../gasmodel.d -I../ $(UTIL_FILES) ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

cea_therm_cond : cea_therm_cond.d therm_cond.d ../gasmodel.d $(UTIL_FILES) $(LIBLUA)
	dmd -main -unittest cea_therm_cond.d therm_cond.d ../gasmodel.d -I../ $(UTIL_FILES) ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

