// l_io.hh

#ifndef L_IO_HH
#define L_IO_HH

int print_simulation_status( FILE *strm, char *efname, int step, simulation_data *SD,
			     vector<slug_data> &A, vector<diaphragm_data> &Diaph,
			     vector<piston_data> &Pist, double cfl_max, 
			     double cfl_tiny, double time_tiny );
int log_event( char *efname, char *event_message );
int L_set_case_parameters(simulation_data *SD, tube_data *T, 
			  ConfigParser &dict, int echo_input);
int set_piston_parameters(struct piston_data *B, int indx, ConfigParser &dict, 
			  double dt_init, int echo_input);
int set_diaphragm_parameters(diaphragm_data *D, int indx,
                             ConfigParser &dict, int echo_input);
int L_set_slug_parameters(slug_data *A, int indx, simulation_data *SD, 
			  ConfigParser &dict, int echo_input);
int L_read_area(struct tube_data *tube, FILE * gf);
int L_write_area(struct tube_data *tube, FILE * gf);
int read_piston_solution(struct piston_data *B, FILE * infile);
int write_piston_solution(struct piston_data *B, FILE * outfile);
int read_diaphragm_solution(struct diaphragm_data *D, FILE * infile);
int write_diaphragm_solution(struct diaphragm_data *D, FILE * outfile);
int L_read_solution(struct slug_data *A, FILE * infile);
int L_write_solution(struct slug_data *A, FILE * outfile, int nsp);
int L_write_cell_history(struct slug_data *A, FILE * hisfile);
int L_write_x_history(double xloc, std::vector<slug_data> &A,
                      int nslug, FILE * hisfile);

#endif
