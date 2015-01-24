# makefile for Eilmer3
# Builds all of the demonstration programs.

DEMO_PROGRAMS := globalconfig_demo flowstate_demo fv_demo block_demo

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
LIBGAS := ${GAS_DIR}/libgas.a

UTIL_DIR := ../util

GEOM_DIR := ../geom
GEOM_FILES := ${GEOM_DIR}/geom.d

GZIP_DIR := ../extern/gzip
GZIP_FILES := ${GZIP_DIR}/gzip.d

LUA_DIR := ../../extern/lua-5.1.4
LIBLUA := ${LUA_DIR}/lib/liblua.a
LIBLUAPATH := ${LUA_DIR}/lib
LUAD_DIR := ../extern/LuaD
LUAD_FILES := ${LUAD_DIR}/luad/*.d \
	${LUAD_DIR}/luad/c/*.d \
	${LUAD_DIR}/luad/conversions/*.d

DFLAGS := -unittest -gc -wi
DLINKFLAGS :=  -w -L-L${LIBLUAPATH} -L-llua -L-ldl

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}
	- rm test-grid.txt.gz test-flow.txt.gz

${LIBGAS}:
	cd ${GAS_DIR}; make libgas.a

${LIBLUA}:
	cd ${LUA_DIR}; make linux local

globalconfig_demo: globalconfig_demo.d globalconfig.d fvcore.d ${GEOM_FILES} ${LIBGAS} ${LIBLUA}
	dmd ${DLINKFLAGS} -I${GEOM_DIR} -I${UTIL_DIR} -I../ \
		globalconfig_demo.d globalconfig.d fvcore.d \
		${GEOM_FILES} ${LIBGAS} ${LUAD_FILES}  

flowstate_demo: flowstate_demo.d flowstate.d fvcore.d \
	${GEOM_FILES} ${LIBGAS} ${LIBLUA}
	dmd ${DLINKFLAGS} -g -I${GEOM_DIR} -I${UTIL_DIR} -I../ \
		flowstate_demo.d fvcore.d flowstate.d \
		${GEOM_FILES} ${LIBGAS} ${LUAD_FILES} 

fv_demo: fv_demo.d fvcore.d fvcell.d fvvertex.d fvinterface.d onedinterp.d fluxcalc.d \
	conservedquantities.d flowstate.d globalconfig.d globaldata.d \
	block.d sblock.d ublock.d \
	${BC_FILES} ${GEOM_FILES} ${LIBGAS} ${LIBLUA} ${GZIP_FILES}
	dmd ${DLINKFLAGS} -g -I${UTIL_DIR} -I${GEOM_DIR} -I${GZIP_DIR} -I../ \
		fv_demo.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		conservedquantities.d flowstate.d onedinterp.d fluxcalc.d \
		globalconfig.d globaldata.d block.d sblock.d ublock.d \
		${BC_FILES} ${GEOM_FILES} ${LIBGAS} ${LUAD_FILES} ${GZIP_FILES} 

block_demo: block_demo.d \
	block.d sblock.d ublock.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d \
	conservedquantities.d flowstate.d globalconfig.d globaldata.d \
	block.d sblock.d ublock.d \
	${BC_FILES} ${GEOM_FILES} ${LIBGAS} ${LIBLUA} ${GZIP_FILES}
	dmd ${DLINKFLAGS} -I${GEOM_DIR} -I${UTIL_DIR} -I${GZIP_DIR} -I../ \
		block_demo.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		conservedquantities.d flowstate.d \
		globalconfig.d globaldata.d \
		block.d sblock.d ublock.d \
		${BC_FILES} ${GEOM_FILES} ${LIBGAS} ${LUAD_FILES} ${GZIP_FILES}
