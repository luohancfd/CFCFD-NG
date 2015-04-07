# makefile for the geometry package
# Builds all of the demonstration programs.

DEMO_PROGRAMS := geom_demo gpath_demo surface_demo univariatefunctions_demo \
	volume_demo \
	luageom_demo luagpath_demo luasurface_demo luaunifunction_demo \
	sgrid_demo luasgrid_demo

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := $(LUA_DIR)/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib
LUAD_DIR := ../extern/LuaD
LUAD_FILES := $(wildcard $(LUAD_DIR)/luad/*.d) \
	$(wildcard $(LUAD_DIR)/luad/c/*.d) \
	$(wildcard $(LUAD_DIR)/luad/conversions/*.d)

DFLAGS := -w -unittest -g -I..
OF := -of
DLINKFLAGS := -L-ldl

DFLAGS_LUA := -w -g -I.. -I$(LUAD_DIR)
DLINKFLAGS_LUA := -L-L$(LIBLUAPATH) -L-llua $(DLINKFLAGS)

demo: $(DEMO_PROGRAMS)
	echo "Demo codes built."

clean:
	- rm *.o
	- rm $(DEMO_PROGRAMS)
	- cd $(LUA_DIR); make clean

$(LIBLUA):
	cd $(LUA_DIR); make linux local

geom_demo: geom_demo.d geom.d
	dmd $(DFLAGS) geom_demo.d $(OF)geom_demo \
		geom.d $(DLINKFLAGS)

gpath_demo: gpath_demo.d gpath.d geom.d
	dmd $(DFLAGS) gpath_demo.d $(OF)gpath_demo \
		gpath.d geom.d $(DLINKFLAGS)

surface_demo: surface_demo.d surface.d gpath.d geom.d
	dmd $(DFLAGS) surface_demo.d $(OF)surface_demo \
		surface.d gpath.d geom.d $(DLINKFLAGS)

volume_demo: volume_demo.d volume.d surface.d gpath.d geom.d
	dmd $(DFLAGS) volume_demo.d $(OF)volume_demo \
		volume.d surface.d gpath.d geom.d $(DLINKFLAGS)

univariatefunctions_demo: univariatefunctions_demo.d univariatefunctions.d
	dmd $(DFLAGS) univariatefunctions_demo.d $(OF)univariatefunctions_demo \
		univariatefunctions.d $(DLINKFLAGS)

sgrid_demo: sgrid_demo.d sgrid.d surface.d gpath.d geom.d univariatefunctions.d \
	$(GZIP_FILES)
	dmd $(DFLAGS) sgrid_demo.d $(OF)sgrid_demo \
		sgrid.d surface.d gpath.d geom.d univariatefunctions.d \
		$(GZIP_FILES) $(DLINKFLAGS)

luageom_demo: luageom_demo.d luageom.d geom.d $(LUAD_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luageom_demo.d $(OF)luageom_demo \
		luageom.d geom.d $(LUAD_FILES) $(DLINKFLAGS_LUA)

luagpath_demo: luagpath_demo.d luagpath.d luageom.d geom.d gpath.d \
	$(UTIL_FILES) $(LUAD_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luagpath_demo.d $(OF)luagpath_demo \
		luagpath.d luageom.d geom.d gpath.d \
		$(UTIL_FILES) $(LUAD_FILES) $(DLINKFLAGS_LUA)

luasurface_demo: luasurface_demo.d luasurface.d luagpath.d luageom.d geom.d gpath.d \
	$(UTIL_FILES) $(LUAD_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luasurface_demo.d $(OF)luasurface_demo \
		luasurface.d luagpath.d luageom.d surface.d geom.d gpath.d \
		$(UTIL_FILES) $(LUAD_FILES) $(DLINKFLAGS_LUA)

luaunifunction_demo: luaunifunction_demo.d luaunifunction.d univariatefunctions.d \
	$(UTIL_FILES) $(LUAD_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luaunifunction_demo.d $(OF)luaunifunction_demo \
		luaunifunction.d univariatefunctions.d \
		$(UTIL_FILES) $(LUAD_FILES) $(DLINKFLAGS_LUA)

luasgrid_demo: luasgrid_demo.d sgrid.d surface.d gpath.d geom.d univariatefunctions.d \
	luasgrid.d luasurface.d luagpath.d luageom.d luaunifunction.d \
	$(GZIP_FILES) $(UTIL_FILES) $(LUAD_FILES)
	dmd $(DFLAGS_LUA) luasgrid_demo.d $(OF)luasgrid_demo \
		sgrid.d surface.d gpath.d geom.d  univariatefunctions.d \
		luasgrid.d luasurface.d luagpath.d luageom.d luaunifunction.d \
		$(GZIP_FILES) $(UTIL_FILES) $(LUAD_FILES) $(DLINKFLAGS_LUA)

