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
# The following is copied from $GAS_DIR/makefile
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

GEOM_DIR := ../geom
GEOM_FILES := ${GEOM_DIR}/geom.d

GZIP_DIR := ../extern/gzip
GZIP_FILES := ${GZIP_DIR}/gzip.d

DFLAGS := -unittest -gc -wi

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}
	- rm test-grid.txt.gz test-flow.txt.gz

globalconfig_demo: globalconfig_demo.d globalconfig.d fvcore.d ${GEOM_FILES}
	dmd ${DFLAGS} -I${GEOM_DIR} -I${GAS_DIR} \
		globalconfig_demo.d globalconfig.d fvcore.d \
		${GEOM_FILES} ${GAS_MODEL_FILES} 

flowstate_demo: flowstate_demo.d flowstate.d fvcore.d \
	${GEOM_FILES} ${GAS_MODEL_FILES}
	dmd -unittest -g -I${GAS_DIR} -I${GEOM_DIR} \
		flowstate_demo.d fvcore.d flowstate.d \
		${GEOM_FILES} ${GAS_MODEL_FILES} 

fv_demo: fv_demo.d fvcore.d fvcell.d fvvertex.d fvinterface.d onedinterp.d fluxcalc.d \
	conservedquantities.d flowstate.d globalconfig.d globaldata.d \
	block.d sblock.d ublock.d \
	${BC_FILES} ${GEOM_FILES} ${GAS_MODEL_FILES} ${GZIP_FILES}
	dmd ${DFLAGS} -I${GAS_DIR} -I${GEOM_DIR} -I${GZIP_DIR} \
		fv_demo.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		conservedquantities.d flowstate.d onedinterp.d fluxcalc.d \
		globalconfig.d globaldata.d block.d sblock.d ublock.d \
		${BC_FILES} ${GEOM_FILES} ${GAS_MODEL_FILES} ${GZIP_FILES} 

block_demo: block_demo.d \
	block.d sblock.d ublock.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d \
	conservedquantities.d flowstate.d globalconfig.d globaldata.d \
	block.d sblock.d ublock.d \
	${BC_FILES} ${GEOM_FILES} ${GAS_MODEL_FILES} ${GZIP_FILES}
	dmd ${DFLAGS} -I${GEOM_DIR} -I${GAS_DIR} -I${GZIP_DIR} \
		block_demo.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		conservedquantities.d flowstate.d \
		globalconfig.d globaldata.d \
		block.d sblock.d ublock.d \
		${BC_FILES} ${GEOM_FILES} ${GAS_MODEL_FILES} ${GZIP_FILES}
