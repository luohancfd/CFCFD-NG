# makefile for Eilmer3
# Builds all of the demonstration programs.

DEMO_PROGRAMS := globaldata_demo cflowcondition_demo

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

globaldata_demo: globaldata_demo.d globaldata.d
	dmd -unittest -g globaldata_demo.d globaldata.d

cflowcondition_demo: cflowcondition_demo.d cflowcondition.d \
	../gas/gas_model.d ../gas/gas_model_util.d ../gas/ideal_gas.d
	dmd -unittest -g -I../gas/ -I../geom/ \
		cflowcondition_demo.d cflowcondition.d \
		../gas/gas_model.d ../gas/gas_model_util.d ../gas/ideal_gas.d \
		../geom/geom.d
