// l_io.hh

#ifndef L_IO_HH
#define L_IO_HH

#include "l_diaph.hh"

int print_simulation_status( FILE *strm, const char* efname, int step, simulation_data *SD,
			     vector<slug_data> &A, vector<DiaphragmData> &Diaph,
			     vector<piston_data> &Pist, double cfl_max, 
			     double cfl_tiny, double time_tiny );
int log_event( const char* efname, const char* event_message );
int L_set_case_parameters(simulation_data *SD, ConfigParser& dict, int echo_input);
int set_piston_parameters(struct piston_data* B, int indx, 
			  ConfigParser& dict, 
			  double dt_init, int echo_input);
int L_set_slug_parameters(slug_data* A, int indx, simulation_data* SD, 
			  ConfigParser& dict, int echo_input);
int read_piston_solution(struct piston_data* B, FILE* infile);
int write_piston_solution(struct piston_data* B, FILE* outfile);
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
