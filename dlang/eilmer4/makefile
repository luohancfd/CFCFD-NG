# makefile for Eilmer3
# Builds all of the demonstration programs.

DEMO_PROGRAMS := globalconfig_demo flowstate_demo fv_demo block_demo

GAS_MODEL_FILES := ../gas/gasmodel.d ../gas/gasmodelutil.d ../gas/idealgas.d

GEOM_FILES := ../geom/geom.d

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

globalconfig_demo: globalconfig_demo.d globalconfig.d ../geom/geom.d
	dmd -unittest -g -I../geom/ -I../gas/ \
		globalconfig_demo.d globalconfig.d \
		${GEOM_FILES} ${GAS_MODEL_FILES} 

flowstate_demo: flowstate_demo.d flowstate.d \
	../gas/gasmodel.d ../gas/gasmodelutil.d ../gas/idealgas.d
	dmd -unittest -g -I../gas/ -I../geom/ \
		flowstate_demo.d flowstate.d \
		${GEOM_FILES} ${GAS_MODEL_FILES} 

fv_demo: fv_demo.d fvcore.d fvcell.d fvvertex.d fvinterface.d \
	conservedquantities.d flowstate.d globalconfig.d globaldata.d block.d \
	../gas/gasmodel.d ../gas/gasmodelutil.d ../gas/idealgas.d ../geom/geom.d
	dmd -unittest -g -I../gas/ -I../geom/ \
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
	../gas/gasmodel.d ../gas/gasmodelutil.d ../gas/idealgas.d ../geom/geom.d
	dmd -unittest -g -I../gas/ -I../geom/ \
		block_demo.d \
		fvcore.d fvcell.d fvvertex.d fvinterface.d \
		conservedquantities.d flowstate.d \
		globalconfig.d globaldata.d \
		block.d sblock.d ublock.d \
		bc.d \
		${GEOM_FILES} ${GAS_MODEL_FILES} 
