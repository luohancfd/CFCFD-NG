# makefile for the gas module

DEMO_PROGRAMS := gasmodel_demo idealgas_demo
GAS_DIR := .
GAS_MODEL_FILES := ${GAS_DIR}/gasmodel.d \
	${GAS_DIR}/idealgas.d \
	${GAS_DIR}/gasmodelutil.d \
	${GAS_DIR}/thermo/thermalEOS.d \
	${GAS_DIR}/thermo/perfectgasEOS.d \
	${GAS_DIR}/thermo/caloricEOS.d \
	${GAS_DIR}/thermo/calperfectgasEOS.d \
	${GAS_DIR}/diffusion/viscosity.d \
	${GAS_DIR}/diffusion/sutherland_visc.d \
	${GAS_DIR}/diffusion/therm_cond.d \
	${GAS_DIR}/diffusion/sutherland_therm_cond.d

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/lua_service.d

LUA := ../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib
LUAD_DIR := ../extern/LuaD
LUAD_FILES := $(LUAD_DIR)/luad/*.d \
	$(LUAD_DIR)/luad/c/*.d \
	$(LUAD_DIR)/luad/conversions/*.d

demo: ${DEMO_PROGRAMS}
	echo "Demo programs built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

$(LIBLUA):
	cd $(LUA); make linux local

gasmodel_demo: gasmodel_demo.d ${GAS_MODEL_FILES} $(LIBLUA)
	dmd -unittest -g gasmodel_demo.d ${GAS_MODEL_FILES} ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

idealgas_demo: idealgas_demo.d ${GAS_MODEL_FILES} $(UTIL_FILES) $(LIBLUA)
	dmd -unittest -g idealgas_demo.d ${GAS_MODEL_FILES} $(UTIL_FILES) ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

