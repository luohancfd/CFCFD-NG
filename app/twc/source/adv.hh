#ifndef ADV_SOLN_HH
#define ADV_SOLN_HH

using namespace std;

// Functions
int flux_calc( list_of_inputs &inputs, list_of_vars &vars, vector<double> q_wall);
int update_coeff_matrix(list_of_vars &vars, list_of_inputs &inputs);
int increment_in_time(list_of_vars &vars, list_of_inputs &inputs);
int read_lookup_flux(list_of_inputs &inputs, list_of_vars &vars); 

#endif

