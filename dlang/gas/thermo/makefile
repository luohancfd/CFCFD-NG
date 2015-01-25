
UTIL_DIR := ../../util
UTIL_FILES := $(UTIL_DIR)/lua_service.d

LUA := ../../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib
LUAD_DIR := ../../extern/LuaD
LUAD_FILES := $(LUAD_DIR)/luad/*.d \
	$(LUAD_DIR)/luad/c/*.d \
	$(LUAD_DIR)/luad/conversions/*.d

GAS_DIR := ..
LIBGAS := ${GAS_DIR}/libgas.a

DLINKFLAGS :=  -w -L-L${LIBLUAPATH} -L-llua -L-ldl

unittest: cea_thermo_curves
	./cea_thermo_curves

${LIBGAS}:
	cd ${GAS_DIR}; make libgas.a

${LIBLUA}:
	cd ${LUA_DIR}; make linux local

cea_thermo_curves : cea_thermo_curves.d $(LIBGAS) $(LIBLUA)
	dmd -main -unittest cea_thermo_curves.d $(LIBGAS) $(LUAD_FILES) -I../../ -I$(UTIL_DIR) $(DLINKFLAGS)
