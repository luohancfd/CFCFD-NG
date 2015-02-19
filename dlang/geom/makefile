# makefile for the geometry package
# Builds all of the demonstration programs.

DEMO_PROGRAMS := geom_demo gpath_demo surface_demo

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}

geom_demo: geom_demo.d geom.d
	dmd -unittest -g geom_demo.d geom.d

gpath_demo: gpath_demo.d gpath.d geom.d
	dmd -unittest -g gpath_demo.d gpath.d geom.d

surface_demo: surface_demo.d surface.d gpath.d geom.d
	dmd -unittest -g surface_demo.d surface.d gpath.d geom.d
