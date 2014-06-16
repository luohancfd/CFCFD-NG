# makefile for the geometry package
# Builds all of the demonstration programs.

DEMO_PROGRAMS := geom_demo

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

geom_demo: geom_demo.d geom.d
	dmd -unittest -g geom_demo.d geom.d
