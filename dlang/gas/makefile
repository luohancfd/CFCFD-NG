# makefile for the gas module

DEMO_PROGRAMS := gas_model_demo ideal_gas_demo

demo: ${DEMO_PROGRAMS}
	echo "Demo programs built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

gas_model_demo: gas_model_demo.d gas_model.d ideal_gas.d
	dmd -unittest -g gas_model_demo.d gas_model.d ideal_gas.d

ideal_gas_demo: ideal_gas_demo.d gas_model.d ideal_gas.d
	dmd -unittest -g ideal_gas_demo.d gas_model.d ideal_gas.d
