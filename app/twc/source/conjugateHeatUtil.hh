#ifndef CHUTIL_HH
#define CHUTIL_HH


void wall_printer(vector<double> &wall, list_of_inputs &inputs);
bool file_exists(const char *filename);
int write_row(list_of_vars &vars, list_of_inputs &inputs, int row);
int write_column(list_of_vars &vars, list_of_inputs &inputs, int column);
int write_soln(list_of_vars &vars, list_of_inputs &inputs, double time_elapsed, int print_number);
int write_hist(list_of_vars &vars, list_of_inputs &inputs, double time_elapsed);
int old_file_delete_catcher(list_of_inputs &inputs) ;
double linear_interp(vector<double> in, vector<double> out, double searchval);
#endif
