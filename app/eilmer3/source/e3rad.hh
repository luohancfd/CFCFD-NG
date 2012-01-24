// e3rad.hh
// Function prototypes for e3rad.

#ifndef E3RAD_HH
#define E3RAD_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

void ensure_directory_is_present( string pathname );
void do_system_cmd( string commandstring );
int prepare_for_radiation_calculation( int start_tindx );
int finalize_e3rad( void );
int write_finishing_data( struct global_data *G, std::string filename );

#endif
