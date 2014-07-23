# makefile for Eilmer3
# Builds all of the demonstration programs.

DEMO_PROGRAMS := globalconfig_demo flowstate_demo fv_demo block_demo

GAS_DIR := ../gas
GAS_MODEL_FILES := ${GAS_DIR}/gasmodel.d ${GAS_DIR}/gasmodelutil.d ${GAS_DIR}/idealgas.d
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

globalconfig_demo: globalconfig_demo.d globalconfig.d ${GEOM_FILES}
	dmd ${DFLAGS} -I${GEOM_DIR} -I${GAS_DIR} \
		globalconfig_demo.d globalconfig.d \
		${GEOM_FILES} ${GAS_MODEL_FILES} 

flowstate_demo: flowstate_demo.d flowstate.d \
	${GEOM_FILES} ${GAS_MODEL_FILES}
	dmd -unittest -g -I${GAS_DIR} -I${GEOM_DIR} \
		flowstate_demo.d flowstate.d \
		${GEOM_FILES} ${GAS_MODEL_FILES} 

fv_demo: fv_demo.d fvcore.d fvcell.d fvvertex.d fvinterface.d \
	conservedquantities.d flowstate.d globalconfig.d globaldata.d block.d \
	${GEOM_FILES} ${GAS_MODEL_FILES}
	dmd ${DFLAGS} -I${GAS_DIR} -I${GEOM_DIR} \
		fv_demo.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		conservedquantities.d flowstate.d \
		globalconfig.d globaldata.d block.d \
		${GEOM_FILES} ${GAS_MODEL_FILES} 

block_demo: block_demo.d \
	block.d sblock.d ublock.d \
	bc.d \
	fvcore.d fvcell.d fvvertex.d fvinterface.d \
	conservedquantities.d flowstate.d globalconfig.d globaldata.d block.d \
	${GEOM_FILES} ${GAS_MODEL_FILES} ${GZIP_FILES}
	dmd ${DFLAGS} -I${GEOM_DIR} -I${GAS_DIR} -I${GZIP_DIR} \
		block_demo.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		conservedquantities.d flowstate.d \
		globalconfig.d globaldata.d \
		block.d sblock.d ublock.d \
		bc.d \
		${GEOM_FILES} ${GAS_MODEL_FILES} ${GZIP_FILES}
