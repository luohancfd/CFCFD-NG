# makefile for the gas module

DFLAGS := -w

DEMO_PROGRAMS := gasmodel_demo idealgas_demo

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/lua_service.d

LUA := ../../extern/lua-5.1.4
LUAD_DIR := ../extern/LuaD
LUAD_FILES := $(LUAD_DIR)/luad/*.d \
	$(LUAD_DIR)/luad/c/*.d \
	$(LUAD_DIR)/luad/conversions/*.d

GAS_MODEL_FILES := \
	package.d \
	gas_model.d \
	ideal_gas.d 

THERMO_FILES := \
	thermo/caloric_EOS.d \
	thermo/cal_perf_gas_EOS.d \
	thermo/cal_perf_gas_mix_EOS.d \
	thermo/perf_gas_EOS.d \
	thermo/perf_gas_mix_EOS.d \
	thermo/thermal_EOS.d

DIFFUSION_FILES := \
	diffusion/cea_therm_cond.d \
	diffusion/cea_viscosity.d \
	diffusion/sutherland_therm_cond.d \
	diffusion/sutherland_viscosity.d \
	diffusion/therm_cond.d \
	diffusion/viscosity.d \
	diffusion/wilke_mixing_therm_cond.d \
	diffusion/wilke_mixing_viscosity.d

GAS_PACKAGE_FILES := $(GAS_MODEL_FILES) $(THERMO_FILES) $(DIFFUSION_FILES)

lib: $(GAS_PACKAGE_FILES) $(LUAD_FILES) $(UTIL_FILES)
	dmd -lib $(GAS_PACKAGE_FILES) $(LUAD_FILES) $(UTIL_FILES) -oflibgas $(DFLAGS)

demo: ${DEMO_PROGRAMS}
	echo "Demo programs built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

gasmodel_demo: gasmodel_demo.d ${GAS_MODEL_FILES} $(LIBLUA)
	dmd -unittest -g gasmodel_demo.d ${GAS_MODEL_FILES} ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

idealgas_demo: idealgas_demo.d ${GAS_MODEL_FILES} $(UTIL_FILES) $(LIBLUA)
	dmd -unittest -g idealgas_demo.d ${GAS_MODEL_FILES} $(UTIL_FILES) ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

