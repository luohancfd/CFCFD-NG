# makefile for the Numerical Methods package
# Builds all of the demonstration programs.

demo: ridder_demo
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ./ridder_demo

ridder_demo: ridder_demo.d ridder.d
	dmd -unittest ridder_demo.d
