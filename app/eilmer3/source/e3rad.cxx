/** \file e3rad.cxx
 * \ingroup eilmer3
 * \brief OpenMP radiation calculator for eilmer3 solutions
 *
 * \author Daniel Potter
 *         See main.cxx for eilmer3 contributors
 *
 * \version 25-Jan-2010 -- ported from main.cxx
 */

/*-----------------------------------------------------------------*/
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif
extern "C" {
#include <popt.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include "../../../lib/util/source/time_to_go.h"
}
#include "../../../lib/util/source/lua_service.hh"
#include "cell.hh"
#include "block.hh"
#include "kernel.hh"
#include "init.hh"
#include "bc.hh"
#include "bc_user_defined.hh"
#include "exch2d.hh"
#include "piston.hh"
#include "radiation_transport.hh"
#include "e3rad.hh"
#include "main.hh"

//-----------------------------------------------------------------
// Global data
//
int history_just_written, output_just_written, av_output_just_written;
int program_return_flag = 0;
int output_counter = 0; // counts the number of flow-solutions written
int zip_files = 1; // flag to indicate if flow and grid files are to be gzipped
int master;
int max_wall_clock = 0;
time_t start, now; // wall-clock timer

lua_State *L; // for the uder-defined procedures

//-----------------------------------------------------------------

void ensure_directory_is_present( string pathname )
{
    string commandstring;
    int cmd_return_status;
    if ( master ) {
	if ( access(pathname.c_str(), F_OK) != 0 ) {
	    commandstring = "mkdir " + pathname;
	    cmd_return_status = system(commandstring.c_str());
	    if ( cmd_return_status != 0 ) {
		cerr << "Problem with system cmd: " << commandstring << endl;
		cerr << "Quitting simulation." << endl;
		exit(FILE_ERROR);
	    }
	}
    }
#   ifdef _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#   endif
    return;
} // end ensure_directory_is_present()


void do_system_cmd( string commandstring )
{
    int cmd_return_status;
    cmd_return_status = system(commandstring.c_str());
    if ( cmd_return_status != 0 ) {
	cerr << "Problem with system cmd: " << commandstring << endl;
	cerr << "Quitting simulation." << endl;
	exit(FILE_ERROR);
    }
    return;
} // end do_system_cmd()

//-----------------------------------------------------------------

void usage(poptContext optCon, int exitcode, char *error, char *addl) {
    // Print program usage message, copied from popt example.
    poptPrintUsage(optCon, stderr, 0);
    if (error) fprintf(stderr, "%s: %s0", error, addl);
    exit(exitcode);
}

//-----------------------------------------------------------------
// Begin here...
//
int main(int argc, char **argv)
{
    global_data &G = *get_global_data_ptr();
    int do_run_simulation = 0;
    int start_tindx = 0; 
    char c, job_name[132], text_buf[132];
    char log_file_name[132];

    poptContext optCon;   /* context for parsing command-line options */
    struct poptOption optionsTable[] = {
	{ "job", 'f', POPT_ARG_STRING, NULL, 'f',
	  "job_name is typically the same as root_file_name", 
	  "<job_name>" },
	{ "run", 'r', POPT_ARG_NONE, NULL, 'r',
	  "run the main simulation time-stepper", 
	  NULL },
	{ "tindx", 't', POPT_ARG_STRING, NULL, 't',
	  "start with this set of flow data", 
	  "<int>" },
	{ "zip-files", 'z', POPT_ARG_NONE, NULL, 'z',
	  "use gzipped flow and grid files", 
	  NULL },
	{ "no-zip-files", 'a', POPT_ARG_NONE, NULL, 'a',
	  "use ASCII (not gzipped) flow and grid files", 
	  NULL },
	{ "no-complain", 'n', POPT_ARG_NONE, NULL, 'n',
	  "suppress complaints about bad data in cells", 
	  NULL },
	{ "verbose", 'v', POPT_ARG_NONE, NULL, 'v',
	  "verbose messages at startup", 
	  NULL },
	{ "max-wall-clock", 'w', POPT_ARG_STRING, NULL, 'w',
	  "maximum wall-clock time in seconds", 
	  "<seconds>" },
	POPT_AUTOHELP
	POPT_TABLEEND
    };

    // ----------
    // INITIALIZE
    // ----------
    program_return_flag = SUCCESS;

    master = 1;
    sprintf(log_file_name, "e3rad.log");
    printf("e3rad: C++,shared-memory version.\n");
#   ifdef _OPENMP
    printf("OpenMP version using %d thread(s).\n", omp_get_max_threads());
#   endif
    if ((G.logfile = fopen(log_file_name, "w")) == NULL) {
        printf( "\nCould not open %s; BAILING OUT\n", log_file_name );
        exit(FILE_ERROR);
    }
    // Configuration.
    //
    printf("e3rad: process command-line options...\n");
    optCon = poptGetContext(NULL, argc, (const char**)argv, optionsTable, 0);
    if (argc < 2) {
	poptPrintUsage(optCon, stderr, 0);
	goto Quit;
    }
    /* Now do options processing as per the popt example. */
    while ((c = poptGetNextOpt(optCon)) >= 0) {
	/* printf("received option %c\n", c); */
	switch (c) {
	case 'f':
	    strcpy(job_name, poptGetOptArg(optCon));
	    break;
	case 'r':
	    do_run_simulation = 1;
	    break;
	case 'z':
	    zip_files = 1;
	    break;
	case 'a':
	    zip_files = 0;
	    break;
	case 't':
	    start_tindx = atoi(poptGetOptArg(optCon));
	    break;
	case 'n':
	    set_bad_cell_complain_flag(0);
	    break;
	case 'v':
	    set_verbose_flag(1);
	    break;
	case 'w':
	    strcpy( text_buf, poptGetOptArg(optCon) );
	    sscanf( text_buf, "%d", &max_wall_clock );
	    break;
	}
    }
    if ( master ) {
	printf( "job_name=%s\n", job_name );
	if ( max_wall_clock > 0 ) {
	    printf( "max_wall_clock=%d seconds\n", max_wall_clock );
	} else {
	    printf( "Run will not be limited by wall clock.\n" );
	}
    }

    G.base_file_name = string(job_name);
    // Read the static configuration parameters from the INI file
    read_config_parameters(G.base_file_name+".config", master ); 
    // Read the time-step control parameters from a separate file.
    // These will also be read at the start of each time step.
    read_control_parameters(G.base_file_name+".control", master, 1);

    // At this point, we know the number of blocks in the calculation.
    // All blocks in same process.
    G.parallel = 0;
    G.rank = 0;
    set_block_range(0, G.nblock - 1);

    // Real work is delegated to functions that know what to do...
    //
    if ( do_run_simulation == 1 ) {
	printf("Run simulation...\n");
	/* The simulation proper. */
	if (prepare_for_radiation_calculation(start_tindx) != SUCCESS) goto Quit;

	RadiationTransportModel * rtm = get_radiation_transport_model_ptr();
	global_data &G = *get_global_data_ptr();
	int jb;
	Block * bdp;
	
	// Determine if a scaled or complete radiation call is required
	if ( ( (G.step / get_radiation_update_frequency()) * 
		get_radiation_update_frequency() == G.step) ) {
	    // recompute
	    rtm->compute_Q_rad_for_flowfield();
	    // store the radiation scaling parameters for each cell
#	    ifdef _OPENMP
#	    pragma omp parallel for private(jb) schedule(runtime)
#	    endif
	    for ( jb = first_block(); jb <= last_block(); ++jb ) {
		bdp = get_block_data_ptr( jb );
		if ( bdp->active != 1 ) continue;
		bdp->apply( &FV_Cell::store_rad_scaling_params, "store-rad-scaling-params" );
	    }
	}
	else {
	    // rescale
#	    ifdef _OPENMP
#	    pragma omp parallel for private(jb) schedule(runtime)
#	    endif
	    for ( jb = first_block(); jb <= last_block(); ++jb ) {
		bdp = get_block_data_ptr(jb);
		if ( bdp->active != 1 ) continue;
		bdp->apply( &FV_Cell::rescale_Q_rE_rad, "rescale-Q_rE_rad" );
	    }
	}
	finalize_e3rad();
    } else {
	printf( "NOTHING DONE -- because you didn't ask...\n" );
    }

    // Finalization.
    //
    Quit: /* nop */;
    fclose(G.logfile);
    eilmer_finalize();
    printf("e3rad: done.\n");
    return program_return_flag;
}   /* end main */

//------------------------------------------------------------------------------

int prepare_for_radiation_calculation( int start_tindx )
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    Block *bdp;
    int jb;
    string filename, commandstring, jbstring, tindxstring;
    char jbcstr[10], tindxcstr[10];
    FILE *fp;

    if ( G.npiston > 0 ) {
	filename = G.base_file_name;
	filename += ".piston0";
	if ((fp = fopen(filename.c_str(), "r")) == NULL) {
	    cerr << "Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	double temporary_time;
	for ( int jp = 0; jp < G.npiston; ++jp ) {
	    temporary_time = G.pistons[jp]->read_state( fp );
	    cout << "temporary_time=" << temporary_time 
		 << " Brendan, FIX-ME please. What are we doing with the piston code?" << endl;
	}
	fclose( fp );
    }

    // Allocate and initialise memory in parallel.
    for ( jb = first_block(); jb <= last_block(); ++jb ) {
	bdp = get_block_data_ptr(jb);
        if ( bdp->array_alloc(G.dimensions) != 0 ) exit( MEMORY_ERROR );
	bdp->bind_interfaces_to_cells(G.dimensions);
    }

    // Read block grid and flow data; write history-file headers.
    // Note that the global simulation time is set by the last flow data read.
    sprintf( tindxcstr, "t%04d", start_tindx);
    tindxstring = tindxcstr;
    for ( jb = first_block(); jb <= last_block(); ++jb ) {
        printf( "----------------------------------\n" );
	bdp = get_block_data_ptr(jb);
	sprintf( jbcstr, ".b%04d", jb );
	jbstring = jbcstr;
	// Read grid from the tindx=0 files, always.
	filename = "grid/t0000/"+G.base_file_name+".grid"+jbstring+".t0000";
        bdp->read_grid(filename, G.dimensions, zip_files);
	// Read flow data from the specified tindx files.
	filename = "flow/"+tindxstring+"/"+G.base_file_name+".flow"+jbstring+"."+tindxstring;
        bdp->read_solution(filename, &(G.sim_time), G.dimensions, zip_files);
	// History file header is only written for a fresh start.
	ensure_directory_is_present("hist");
	filename = "hist/" + G.base_file_name + ".hist"+jbstring;
	if ( access(filename.c_str(), F_OK) != 0 ) {
	    // History file does not yet exist; write header.
	    bdp->write_history( filename, G.sim_time, 1 );
	}
    }
    output_counter = start_tindx;
    filename = G.base_file_name; filename += ".times";
    if ((G.timestampfile = fopen(filename.c_str(), "a")) == NULL) {
        cerr << "Could not open " << filename << "; BAILING OUT" << endl;
        exit(FILE_ERROR);
    }
    // The zeroth entry in the timestampfile has been written already by e3prep.py.
    // fprintf( G.timestampfile, "# count sim_time dt_global\n" );
    // fprintf( G.timestampfile, "%04d %e %e\n", output_counter, G.sim_time, G.dt_global );
    // fflush( G.timestampfile );

    // Prepare data within the primary finite-volume cells.
    // This includes both geometric data and flow state data.
    for (jb = first_block(); jb <= last_block(); ++jb) {
	bdp = get_block_data_ptr(jb);
	bdp->compute_primary_cell_geometric_data(G.dimensions);
	bdp->compute_distance_to_nearest_wall_for_all_cells(G.dimensions);
	bdp->compute_secondary_cell_geometric_data(G.dimensions);
	bdp->set_base_qdot( G );  // this need be done only once per block
	bdp->identify_reaction_zones( G );
	bdp->identify_turbulent_zones( G );
        bdp->apply( &FV_Cell::encode_conserved, bdp->omegaz, "encode_conserved" );
        // Even though the following call appears redundant at this point,
        // fills in some gas properties such as Prandtl number that is
        // needed for both the cfd_check and the BLomax turbulence model.
        bdp->apply( &FV_Cell::decode_conserved, bdp->omegaz, "decode_conserved" );
    }

    // Exchange boundary cell geometry information so that we can
    // next calculate secondary-cell geometries.
    // NOTE: now copying all data as required for radiation calculation
    for (jb = first_block(); jb <= last_block(); ++jb) {
        exchange_shared_boundary_data(jb, COPY_ALL_CELL_DATA);
    }

    // Start up the Lua interpreter and load the external file
    // containing the user-defined functions -- if appropriate.
    if ( G.udf_file.length() > 0 ) {
	L = luaL_newstate();
	luaL_openlibs(L); // load the standard libraries

	// Set up some global variables that might be handy in 
	// the Lua environment.
	lua_pushinteger(L, G.nblock);
	lua_setglobal(L, "nblock");
	int nsp = get_gas_model_ptr()->get_number_of_species();
	int nmodes = get_gas_model_ptr()->get_number_of_modes();
	lua_pushinteger(L, nsp);
	lua_setglobal(L, "nsp");
	lua_pushinteger(L, nmodes);
	lua_setglobal(L, "nmodes");

	// Register functions so that they are accessible 
	// from the Lua environment.
	lua_pushcfunction(L, luafn_sample_flow);
	lua_setglobal(L, "sample_flow");
	lua_pushcfunction(L, luafn_locate_cell);
	lua_setglobal(L, "locate_cell");

	// Presume that the user-defined functions are in the specified file.
	if ( luaL_loadfile(L, G.udf_file.c_str()) || lua_pcall(L, 0, 0, 0) ) {
	    handle_lua_error(L, "Could not run user file: %s", lua_tostring(L, -1));
	}
	lua_settop(L, 0); // clear the stack
    }
    
    // Initialise radiation transport if appropriate
    if ( get_radiation_flag() ) {
    	if ( get_radiation_transport_model_ptr()->initialise() ) {
	    cerr << "Problem with initialisation of radiation transport data\n";
	    cerr << "Exiting program." << endl;
	    exit(FAILURE);
	}
    }
    else {
    	cerr << "No radiation models present.\n";
	cerr << "Quitting simulation." << endl;
	exit(FAILURE);
    }
    
    // Apply viscous THEN inviscid boundary conditions to match environment for
    // radiation calculation in gasdynamic_inviscid_increment()
    for ( jb = first_block(); jb <= last_block(); ++jb ) {
	bdp = get_block_data_ptr(jb);
	if ( bdp->active != 1 ) continue;
	if ( get_viscous_flag() ) {
	    apply_viscous_bc( *bdp, G.sim_time, G.dimensions );
	}
	apply_inviscid_bc( *bdp, G.sim_time, G.dimensions );
    }
    
    start = time(NULL); // start of wallclock timing
    
    return SUCCESS;
} // end prepare_for_radiation_calculation()

//---------------------------------------------------------------------------

int finalize_e3rad( void )
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    Block *bdp;
    string filename, commandstring, foldername, jbstring, jsstring;
    char jbcstr[10],jscstr[10];
    int jb, js;
    int final_s;

    // Write out the final solution even if it has just been written
    // as part of the main time-stepping loop.
    output_counter = 9999;
    fprintf( G.timestampfile, "%04d %e %e\n", output_counter, G.sim_time, G.dt_global );
    fflush( G.timestampfile );
    foldername = "flow/t9999";
    ensure_directory_is_present(foldername);
    for ( jb = first_block(); jb <= last_block(); ++jb ) {
	bdp = get_block_data_ptr(jb);
	sprintf( jbcstr, ".b%04d", jb ); jbstring = jbcstr;
	filename = foldername+"/"+G.base_file_name+".flow"+jbstring+".t9999";
	bdp->write_solution(filename, G.sim_time, G.dimensions, zip_files);
    }
    // Compute, store and write heat-flux data, if viscous simulation
    if ( get_viscous_flag() ) {
	foldername = "heat/t9999";
	ensure_directory_is_present(foldername);
	// Loop over blocks
	for ( jb = first_block(); jb <= last_block(); ++jb ) {
	    bdp = get_block_data_ptr(jb);
	    sprintf( jbcstr, ".b%04d", jb ); jbstring = jbcstr;
	    final_s = ((G.dimensions == 3)? BOTTOM : WEST);
	    // Loop over boundaries/surfaces
	    for ( js = NORTH; js <= final_s; ++js ) {
		sprintf( jscstr, ".s%04d", js ); jsstring = jscstr;
		filename = foldername+"/"+ G.base_file_name+".heat" \
				+jbstring+jsstring+".t9999";
		bdp->bcp[js]->compute_surface_heat_flux();
		bdp->bcp[js]->write_surface_heat_flux(filename,G.sim_time);
		if ( zip_files == 1 ) do_system_cmd("gzip -f "+filename);
	    }
	}
    }
    output_just_written = 1;
    // For the history files, we don't want to double-up on solution data.
    if ( !history_just_written ) {
        for ( jb = first_block(); jb <= last_block(); ++jb ) {
	    bdp = get_block_data_ptr(jb);
	    sprintf( jbcstr, ".b%04d", jb ); jbstring = jbcstr;
	    filename = "hist/"+G.base_file_name+".hist"+jbstring;
            bdp->write_history( filename, G.sim_time );
	}
	for ( int jp = 0; jp < G.npiston; ++jp ) {
#           if 0
	    G.pistons[jp]->write_state( piston_out_file, G.sim_time );
#           else
	    printf( "Brendan, please fix piston stuff.\n" );
	    exit( NOT_IMPLEMENTED_ERROR );
#           endif
	}
        history_just_written = 1;
    }
    printf( "\nTotal number of steps = %d\n", G.step );

    filename = G.base_file_name; filename += ".finish";
    write_finishing_data( &G, filename );
    fclose( G.timestampfile );

    for ( jb = first_block(); jb <= last_block(); ++jb ) {
	bdp = get_block_data_ptr(jb);
	bdp->array_cleanup(G.dimensions);
    }
    return SUCCESS;
} // end finalize_e3rad()

/// \brief Write out simulation data at end of simulation, such as: no. steps, final time.
int write_finishing_data( struct global_data *G, std::string filename )
{
    FILE *fp;
    if ((fp = fopen(filename.c_str(), "a")) == NULL) {
	cerr << "Could not open " << filename << "; BAILING OUT" << endl;
	exit( FILE_ERROR );
    }
    fprintf(fp, "[simulation_end]\n");
    fprintf(fp, "final_time = %.12e \n", G->sim_time);
    fprintf(fp, "dt = %.12e \n", G->dt_allow);
    fprintf(fp, "no_steps = %d\n", G->step);
    fclose( fp );
    return SUCCESS;

}

/*============= End of e3rad.cxx ===============*/
