# makefile for kinetics module
# Can be used to perform unit tests
# and build stand-alone programs.

UNITTEST_FILES := rate_constant

GAS_DIR := ../gas
LIBGAS := ${GAS_DIR}/libgas.a

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := ${LUA_DIR}/lib/liblua.a
LIBLUAPATH := ${LUA_DIR}/lib
LUAD_DIR := ../extern/LuaD
LUAD_FILES := ${LUAD_DIR}/luad/*.d \
	${LUAD_DIR}/luad/c/*.d \
	${LUAD_DIR}/luad/conversions/*.d

DFLAGS := -main -unittest -gc -wi -I../ 
DLINKFLAGS :=  -w -L-L${LIBLUAPATH} -L-llua -L-ldl

unittests: $(UNITTEST_FILES)

${LIBGAS}:
	cd ${GAS_DIR}; make libgas.a

${LIBLUA}:
	cd ${LUA_DIR}; make linux local

rate_constant: rate_constant.d $(LIBGAS) $(LIBLUA) $(LUAD_FILES)
	dmd $(DFLAGS) rate_constant.d ${LIBGAS} ${LIBLUA} $(LUAD_FILES) $(DLINKFLAGS)
