// main.hh
// Function prototype for those functions near the main function.

#ifndef MAIN_HH
#define MAIN_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#define CHECK_RADIATION_SCALING 0

void ensure_directory_is_present( string pathname );
void do_system_cmd( string commandstring );
int prepare_to_integrate( int start_tindx );
int call_udf( double t, int step, std::string udf_fn_name );
int udf_source_vector_for_cell( FV_Cell *cell, double t );
int integrate_blocks_in_sequence( void );
int integrate_in_time( double target_time );
int finalize( void );
int gasdynamic_inviscid_increment( void );
int gasdynamic_viscous_increment( void );
int do_bad_cell_count( void );
int write_finishing_data( global_data *G, std::string filename );
int check_radiation_scaling( void );

#endif
