# makefile for Eilmer4
# Builds main simulation program by default.

# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= dmd

PROGRAMS := e4shared e4prep

DEMO_PROGRAMS := globalconfig_demo flowstate_demo fv_demo block_demo \
	luaflowstate_demo

BC_FILES := bc.d \
	bc_slip_wall.d \
	bc_adiabatic_wall.d \
	bc_fixed_t_wall.d \
	bc_full_face_exchange.d \
	bc_mapped_cell_exchange.d \
	bc_extrapolate_out.d \
	bc_fixed_p_out.d \
	bc_supersonic_in.d \
	bc_menter_correction.d

GAS_DIR := ../gas

GAS_MODEL_FILES := \
	$(GAS_DIR)/package.d \
	$(GAS_DIR)/gas_model.d \
	$(GAS_DIR)/gas_model_util.d \
	$(GAS_DIR)/ideal_gas.d \
	$(GAS_DIR)/physical_constants.d \
	$(GAS_DIR)/therm_perf_gas.d

THERMO_FILES := \
	$(GAS_DIR)/thermo/cea_thermo_curves.d \
	$(GAS_DIR)/thermo/evt_eos.d \
	$(GAS_DIR)/thermo/perf_gas_mix_eos.d \
	$(GAS_DIR)/thermo/pvt_eos.d \
	$(GAS_DIR)/thermo/therm_perf_gas_mix_eos.d

DIFFUSION_FILES := \
	$(GAS_DIR)/diffusion/cea_therm_cond.d \
	$(GAS_DIR)/diffusion/cea_viscosity.d \
	$(GAS_DIR)/diffusion/sutherland_therm_cond.d \
	$(GAS_DIR)/diffusion/sutherland_viscosity.d \
	$(GAS_DIR)/diffusion/therm_cond.d \
	$(GAS_DIR)/diffusion/viscosity.d \
	$(GAS_DIR)/diffusion/wilke_mixing_therm_cond.d \
	$(GAS_DIR)/diffusion/wilke_mixing_viscosity.d

GAS_FILES := $(GAS_MODEL_FILES) $(THERMO_FILES) $(DIFFUSION_FILES)

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

NM_DIR := ../nm
NM_FILES := $(NM_DIR)/ridder.d

GEOM_DIR := ../geom
GEOM_FILES := $(GEOM_DIR)/geom.d \
	$(GEOM_DIR)/univariatefunctions.d \
	$(GEOM_DIR)/gpath.d \
	$(GEOM_DIR)/surface.d \
	$(GEOM_DIR)/luageom.d \
	$(GEOM_DIR)/luagpath.d \
	$(GEOM_DIR)/luasurface.d

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := $(LUA_DIR)/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib
LUAD_DIR := ../extern/LuaD
LUAD_FILES := $(wildcard $(LUAD_DIR)/luad/*.d) \
	$(wildcard $(LUAD_DIR)/luad/c/*.d) \
	$(wildcard $(LUAD_DIR)/luad/conversions/*.d)

INSTALL_DIR ?= $(HOME)/e3bin
REVISION_STRING := $(shell hg identify --id --num --branch --tags)

ifeq ($(DMD), dmd)
    # DFLAGS := -unittest -gc -wi
    # DFLAGS := -profile -vgc -w -O -release -inline -boundscheck=off
    DFLAGS := -gc -w -O -release -inline -boundscheck=off \
	-I.. -I$(LUAD_DIR) -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GZIP_DIR)
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), ldmd2)
    DFLAGS := -w -O -release -inline -boundscheck=off \
	-I.. -I$(LUAD_DIR) -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GZIP_DIR)
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl
endif
ifeq ($(DMD), gdc)
    DFLAGS := -Wall -O2 -frelease -fno-debug \
	-I.. -I$(LUAD_DIR) -I$(NM_DIR) -I$(UTIL_DIR) -I$(GEOM_DIR) -I$(GZIP_DIR)
    OF := -o
    DLINKFLAGS := -L$(LIBLUAPATH) $(LIBLUA) -ldl
endif

default: $(PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Eilmer4 simulation code built."

demo: $(DEMO_PROGRAMS)
	@echo "Source code revision string $(REVISION_STRING)"
	@echo "Demo codes built."

clean:
	- rm *.o
	- rm $(DEMO_PROGRAMS) $(PROGRAMS)
	- rm test-grid.txt.gz test-flow.txt.gz trace.log trace.def
	- rm -r plot/ flow/ grid/
	- rm e4_main_with_rev_string.d e4prep_with_rev_string.d
	- cd $(LUA_DIR); make clean
	- cd $(GAS_DIR); make clean; rm libgas.a

$(LIBLUA):
	cd $(LUA_DIR); make linux local

globalconfig_demo: globalconfig_demo.d globalconfig.d fvcore.d \
	$(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES) $(LUAD_FILES)
	$(DMD) $(DFLAGS) $(OF)globalconfig_demo \
		globalconfig_demo.d globalconfig.d fvcore.d \
		$(GEOM_FILES) $(GAS_FILES) $(LUAD_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)  

flowstate_demo: flowstate_demo.d flowstate.d fvcore.d json_helper.d \
	$(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES) $(LUAD_FILES)
	$(DMD) $(DFLAGS) $(OF)flowstate_demo \
		flowstate_demo.d fvcore.d flowstate.d json_helper.d \
		$(GEOM_FILES) $(GAS_FILES) $(LUAD_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

luaflowstate_demo: luaflowstate_demo.d luaflowstate.d flowstate.d json_helper.d \
	$(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES) $(LUAD_FILES)
	$(DMD) $(DFLAGS) $(OF)luaflowstate_demo \
		luaflowstate_demo.d luaflowstate.d flowstate.d json_helper.d \
		$(GEOM_FILES) $(GAS_FILES) $(LUAD_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)



fv_demo: fv_demo.d fvcore.d fvcell.d fvvertex.d fvinterface.d onedinterp.d fluxcalc.d \
	conservedquantities.d flowstate.d globalconfig.d globaldata.d \
	block.d sblock.d ublock.d json_helper.d \
	$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(UTIL_FILES) $(NM_FILES) $(LUAD_FILES)
	$(DMD) $(DFLAGS) $(OF)fv_demo \
		fv_demo.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		conservedquantities.d flowstate.d onedinterp.d fluxcalc.d \
		globalconfig.d globaldata.d block.d sblock.d ublock.d \
		json_helper.d \
		$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(LUAD_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

block_demo: block_demo.d \
	block.d sblock.d ublock.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d \
	fluxcalc.d onedinterp.d \
	conservedquantities.d flowstate.d globalconfig.d globaldata.d \
	block.d sblock.d ublock.d \
	json_helper.d \
	$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(UTIL_FILES) $(NM_FILES) $(LUAD_FILES)
	$(DMD) $(DFLAGS) $(OF)block_demo \
		block_demo.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		fluxcalc.d onedinterp.d \
		conservedquantities.d flowstate.d \
		globalconfig.d globaldata.d \
		block.d sblock.d ublock.d \
		json_helper.d \
		$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(LUAD_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

e4shared: e4_main.d e4_core.d readconfig.d json_helper.d \
	block.d sblock.d ublock.d \
	fluxcalc.d onedinterp.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d \
	conservedquantities.d flowstate.d globalconfig.d globaldata.d \
	block.d sblock.d ublock.d \
	$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(GZIP_FILES) \
	$(UTIL_FILES) $(NM_FILES) $(LUAD_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' e4_main.d > e4_main_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)e4shared \
		e4_main_with_rev_string.d e4_core.d readconfig.d json_helper.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		fluxcalc.d onedinterp.d \
		conservedquantities.d flowstate.d \
		globalconfig.d globaldata.d \
		block.d sblock.d ublock.d \
		$(BC_FILES) $(GEOM_FILES) $(GAS_FILES) $(LUAD_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)

e4prep: e4prep.d luaflowstate.d flowstate.d json_helper.d \
	$(GEOM_FILES) $(GAS_FILES) $(LIBLUA) $(UTIL_FILES) $(NM_FILES) $(LUAD_FILES) \
	$(GZIP_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' e4prep.d > e4prep_with_rev_string.d
	$(DMD) $(DFLAGS) $(OF)e4prep \
		e4prep_with_rev_string.d luaflowstate.d flowstate.d json_helper.d \
		$(GEOM_FILES) $(GAS_FILES) $(LUAD_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(DLINKFLAGS)
