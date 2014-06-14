# makefile for the Numerical Methods package
# Builds all of the demonstration programs.

DEMO_PROGRAMS := ridder_demo linesearch_demo nelmin_demo

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

ridder_demo: ridder_demo.d ridder.d
	dmd -unittest -g ridder_demo.d

# Don't know why we have to include linesearch.d on the command line.
linesearch_demo: linesearch_demo.d linesearch.d
	dmd -unittest -g linesearch_demo.d linesearch.d

nelmin_demo: nelmin_demo.d nelmin.d
	dmd -unittest -g nelmin_demo.d nelmin.d
