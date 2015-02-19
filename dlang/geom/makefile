# makefile for the geometry package
# Builds all of the demonstration programs.

DEMO_PROGRAMS := geom_demo gpath_demo surface_demo

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := $(LUA_DIR)/lib/liblua.a
LIBLUAPATH := $(LUA_DIR)/lib
LUAD_DIR := ../extern/LuaD
LUAD_FILES := $(wildcard $(LUAD_DIR)/luad/*.d) \
	$(wildcard $(LUAD_DIR)/luad/c/*.d) \
	$(wildcard $(LUAD_DIR)/luad/conversions/*.d)

DFLAGS := -w -unittest -g -I.. -I$(LUAD_DIR)
OF := -of
DLINKFLAGS := -L-L$(LIBLUAPATH) -L-llua -L-ldl

demo: $(DEMO_PROGRAMS)
	echo "Demo codes built."

clean:
	- rm *.o
	- rm $(DEMO_PROGRAMS)
	- cd $(LUA_DIR); make clean

$(LIBLUA):
	cd $(LUA_DIR); make linux local

geom_demo: geom_demo.d geom.d $(LUAD_FILES) $(LIBLUA)
	dmd $(DFLAGS) geom_demo.d $(OF)geom_demo \
		geom.d $(LUAD_FILES) $(DLINKFLAGS)

gpath_demo: gpath_demo.d gpath.d geom.d $(LUAD_FILES) $(LIBLUA)
	dmd $(DFLAGS) gpath_demo.d $(OF)gpath_demo \
		gpath.d geom.d $(LUAD_FILES) $(DLINKFLAGS)

surface_demo: surface_demo.d surface.d gpath.d geom.d $(LUAD_FILES) $(LIBLUA)
	dmd $(DFLAGS) surface_demo.d $(OF)surface_demo \
		surface.d gpath.d geom.d $(LUAD_FILES) $(DLINKFLAGS)
