// main.hh
// Function prototype for those functions near the main function.

#ifndef MAIN_HH
#define MAIN_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

void ensure_directory_is_present(string pathname);
void do_system_cmd(string commandstring);
int prepare_to_integrate(size_t start_tindx);
int call_udf(double t, size_t step, std::string udf_fn_name);
int add_udf_source_vector_for_cell(FV_Cell *cell, size_t gtl, double t);
int integrate_blocks_in_sequence(void);
int write_solution_data(std::string tindxstring);
int integrate_in_time(double target_time);
int finalize_simulation(void);
int gasdynamic_explicit_increment_with_fixed_grid(double dt);
int gasdynamic_increment_with_moving_grid(double dt, bool &finished_time_stepping, bool &do_cfl_check_now);
int gasdynamic_separate_explicit_viscous_increment();
int do_bad_cell_count(size_t gtl);
int write_finishing_data(global_data *G, std::string filename);
int check_radiation_scaling(void);
int radiation_calculation(void);
void perform_radiation_transport(void);

#endif
