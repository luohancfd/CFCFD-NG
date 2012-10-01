// l_io.hh

#ifndef L_IO_HH
#define L_IO_HH

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "l_kernel.hh"
#include "l_diaph.hh"
#include "l_piston.hh"
#include "l_slug.hh"

int print_simulation_status(FILE *strm, const char* efname, int step, SimulationData& SD,
			    std::vector<GasSlug>& A, std::vector<DiaphragmData>& Diaph,
			    std::vector<PistonData>& Pist, double cfl_max, 
			    double cfl_tiny, double time_tiny );
int log_event( const char* efname, const char* event_message );
int L_write_cell_history(GasSlug& A, FILE * hisfile);
int L_write_x_history(double xloc, std::vector<GasSlug>& A, FILE* hisfile);

#endif
