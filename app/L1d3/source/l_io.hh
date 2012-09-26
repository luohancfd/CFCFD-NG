// l_io.hh

#ifndef L_IO_HH
#define L_IO_HH

#include "l_kernel.hh"
#include "l_diaph.hh"
#include "l_piston.hh"

int print_simulation_status( FILE *strm, const char* efname, int step, SimulationData& SD,
			     vector<slug_data>& A, vector<DiaphragmData>& Diaph,
			     vector<PistonData>& Pist, double cfl_max, 
			     double cfl_tiny, double time_tiny );
int log_event( const char* efname, const char* event_message );
int L_set_slug_parameters(slug_data* A, int indx, SimulationData& SD, 
			  ConfigParser& dict, int echo_input);
std::string write_iface_values_to_string(struct L_cell& cell);
int scan_iface_values_from_string(char* bufptr, struct L_cell& cell);
std::string write_cell_values_to_string(struct L_cell& cell);
int scan_cell_values_from_string(char* bufptr, struct L_cell& cell);
int L_read_solution(struct slug_data*A, FILE* infile);
int L_write_solution(struct slug_data* A, FILE* outfile);
int L_write_cell_history(struct slug_data* A, FILE * hisfile);
int L_write_x_history(double xloc, std::vector<slug_data> &A,
                      int nslug, FILE* hisfile);

#endif
