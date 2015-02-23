# makefile for the geometry package
# Builds all of the demonstration programs.

DEMO_PROGRAMS := geom_demo gpath_demo surface_demo univariatefunctions_demo \
	luageom_demo luagpath_demo sgrid_demo

GZIP_DIR := ../extern/gzip
GZIP_FILES := $(GZIP_DIR)/gzip.d

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

DFLAGS_LUA := $(DFLAGS) -I$(LUAD_DIR)
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

luagpath_demo: luagpath_demo.d luagpath.d luageom2.d geom.d gpath.d $(LUAD_FILES) $(LIBLUA)
	dmd $(DFLAGS_LUA) luagpath_demo.d $(OF)luagpath_demo \
		luagpath.d luageom2.d geom.d gpath.d $(LUAD_FILES) $(DLINKFLAGS_LUA)
