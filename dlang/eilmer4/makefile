# makefile for Eilmer3
# Builds all of the demonstration programs.

DEMO_PROGRAMS := globalconfig_demo flowstate_demo

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

globalconfig_demo: globalconfig_demo.d globalconfig.d
	dmd -unittest -g globalconfig_demo.d globalconfig.d

flowstate_demo: flowstate_demo.d flowstate.d \
	../gas/gasmodel.d ../gas/gasmodelutil.d ../gas/idealgas.d
	dmd -unittest -g -I../gas/ -I../geom/ \
		flowstate_demo.d flowstate.d \
		../gas/gasmodel.d ../gas/gasmodelutil.d ../gas/idealgas.d \
		../geom/geom.d
