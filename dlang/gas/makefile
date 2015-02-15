# makefile for the gas module

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.
DMD ?= dmd

DEMO_PROGRAMS := gas_model_demo ideal_gas_demo
UNITTEST_FILES := therm_perf_gas

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

NM_DIR := ../nm
NM_FILES := $(NM_DIR)/ridder.d

LUA := ../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib
LUAD_DIR := ../extern/LuaD
LUAD_FILES := $(LUAD_DIR)/luad/*.d \
	$(LUAD_DIR)/luad/c/*.d \
	$(LUAD_DIR)/luad/conversions/*.d

GAS_MODEL_FILES := \
	package.d \
	gas_model.d \
	gas_model_util.d \
	ideal_gas.d \
	physical_constants.d \
	therm_perf_gas.d

THERMO_FILES := \
	thermo/caloric_EOS.d \
	thermo/cal_perf_gas_EOS.d \
	thermo/cal_perf_gas_mix_EOS.d \
	thermo/cea_thermo_curves.d \
	thermo/perf_gas_EOS.d \
	thermo/perf_gas_mix_EOS.d \
	thermo/thermal_EOS.d \
	thermo/therm_perf_gas_mix_EOS.d

DIFFUSION_FILES := \
	diffusion/cea_therm_cond.d \
	diffusion/cea_viscosity.d \
	diffusion/sutherland_therm_cond.d \
	diffusion/sutherland_viscosity.d \
	diffusion/therm_cond.d \
	diffusion/viscosity.d \
	diffusion/wilke_mixing_therm_cond.d \
	diffusion/wilke_mixing_viscosity.d

GAS_FILES := $(GAS_MODEL_FILES) $(THERMO_FILES) $(DIFFUSION_FILES)

ifeq ($(DMD), dmd)
    # DFLAGS := -w
    DFLAGS := -O -release -inline -boundscheck=off -w
    DLINKFLAGS := -w -I../ -I../util -L-L$(LIBLUAPATH) -L-llua -L-ldl
    LIBNAME := libgas
endif
ifeq ($(DMD), ldmd2)
    # DFLAGS := -w
    DFLAGS := -O -release -inline -boundscheck=off -w
    DLINKFLAGS := -w -I../ -I../util -L-L$(LIBLUAPATH) -L-llua -L-ldl
    LIBNAME := ./gas
endif

libgas.a : $(GAS_FILES) $(LUAD_FILES) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -lib $(GAS_FILES) $(LUAD_FILES) $(UTIL_FILES) $(NM_FILES) -of$(LIBNAME) $(DFLAGS)

demo: ${DEMO_PROGRAMS}
	echo "Demo programs built."

clean:
	- rm ${DEMO_PROGRAMS} *.o *.a
	- cd $(LUA); make clean

$(LIBLUA):
	cd $(LUA); make linux local

gas_model_demo: gas_model_demo.d $(LUAD_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g gas_model_demo.d $(GAS_FILES) $(LUAD_FILES) $(DLINKFLAGS) $(UTIL_FILES) $(NM_FILES)

ideal_gas_demo: ideal_gas_demo.d $(LUAD_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -g ideal_gas_demo.d $(GAS_FILES) $(LUAD_FILES) $(DLINKFLAGS) $(UTIL_FILES) $(NM_FILES)

therm_perf_gas : therm_perf_gas.d $(LUAD_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES)
	$(DMD) -main -unittest therm_perf_gas.d $(GAS_FILES) $(LUAD_FILES) $(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)
