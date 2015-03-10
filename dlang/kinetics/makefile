# makefile for kinetics module
# Can be used to perform unit tests
# and build stand-alone programs.

INSTALL_DIR ?= $(HOME)/e3bin

UNITTEST_FILES := rate_constant_test \
	reaction_test

UTIL_DIR := ../util
UTIL_FILES := $(UTIL_DIR)/lua_service.d \
	$(UTIL_DIR)/msg_service.d

NM_DIR := ../nm
NM_FILES := $(NM_DIR)/ridder.d

GAS_DIR := ../gas
LIBGAS := ${GAS_DIR}/libgas.a

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := ${LUA_DIR}/lib/liblua.a
LIBLUAPATH := ${LUA_DIR}/lib
LUAD_DIR := ../extern/LuaD
LUAD_FILES := ${LUAD_DIR}/luad/*.d \
	${LUAD_DIR}/luad/c/*.d \
	${LUAD_DIR}/luad/conversions/*.d

KINETICS_FILES := \
	chemistry_update.d \
	rate_constant.d \
	reaction.d \
	reaction_mechanism.d


DFLAGS := -main -unittest -gc -wi -I../ 
DLINKFLAGS :=  -w -L-L${LIBLUAPATH} -L-llua -L-ldl

install-prep-chem: prep_chem.lua reaction.lua lex_elems.lua
	cp reaction.lua lex_elems.lua $(INSTALL_DIR)
	cp prep_chem.lua $(INSTALL_DIR)/prep-chem; chmod +x $(INSTALL_DIR)/prep-chem

unittests: $(UNITTEST_FILES)

${LIBGAS}:
	cd ${GAS_DIR}; make libgas.a

${LIBLUA}:
	cd ${LUA_DIR}; make linux local

libkinetics.a : $(KINETICS_FILES) $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LUAD_FILES)
	dmd -lib -oflibkinetics.a $(KINETICS_FILES) $(GAS_FILES) $(UTIL_FILES) $(NM_FILES) $(LUAD_FILES) $(DFLAGS)

chemistry_update_test: chemistry_update.d libkinetics.a $(LIBGAS) $(LIBLUA) $(LUAD_FILES) $(NM_FILES)
	dmd $(DFLAGS) -ofchemistry_update_test chemistry_update.d libkinetics.a ${LIBGAS} ${LIBLUA} $(LUAD_FILES) $(NM_FILES) $(DLINKFLAGS)

rate_constant_test: rate_constant.d $(LIBGAS) $(LIBLUA) $(LUAD_FILES)
	dmd $(DFLAGS) -ofrate_constant_test rate_constant.d ${LIBGAS} ${LIBLUA} $(LUAD_FILES) $(DLINKFLAGS)

reaction_test: reaction.d $(LIBGAS) $(LIBLUA) $(LUAD_FILES)
	dmd $(DFLAGS) -ofreaction_test reaction.d rate_constant.d ${LIBGAS} ${LIBLUA} $(LUAD_FILES) $(DLINKFLAGS)

reaction_mechanism_test: reaction_mechanism.d libkinetics.a $(LIBGAS) $(LIBLUA) $(LUAD_FILES)
	dmd $(DFLAGS) -ofreaction_mechanism_test reaction.d rate_constant.d ${LIBGAS} ${LIBLUA} $(LUAD_FILES) $(DLINKFLAGS)


