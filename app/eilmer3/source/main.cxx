/** \file main.cxx
 * \ingroup eilmer3
 * \brief Multiple-Block Complete Navier-Stokes -- main program.
 *
 * \author PA Jacobs for the core stepper, geometry description and
 *         the boring house-keeping bits such as macros in the viscous derivatives.
 *         Various postgrad students and post docs for the interesting bits.
 *         Andrew McGhee and Richard Goozee -- MPI parallel computing.
 *         Chris Craddock and Ian Johnston -- first crack at thermochemistry, 3D
 *         Paul Petrie and Ian Johnston -- flux calculators, grid management
 *         Rowan Gollan  -- thermochemistry, radiation
 *         Andrew Denman -- 3D viscous effects and LES turbulence
 *         Joseph Tang -- hierarchical grids
 *         Brendan O'Flaherty -- piston, thermochemistry, real gas models
 *         Jan-Pieter Nap -- help with the k-omega turbulence
 *         Wilson Chan -- more work on the turbulence models
 *         Dan Potter -- radiation, thermochemistry
 *
 *         ...and quite a few people have spent lots of effort "kicking the tyres",
 *         suffering through buggy versions of new features and generally helping
 *         by making new and interesting applications.
 *         Mike Wendt -- lots of enthusiasm in the early days when the code was
 *                       in its most buggy state
 *         Rainer Kirchhartz -- patience with the turbulent scramjet models
 *         Katsu Tanimizu -- the 3D scramjet work
 *
 * \version 20-Mar-2006 -- New-generation using some C++.
 * \version    Jan-2008 -- split out material to cns_kernel.cxx, much like Elmer2.
 * \version 25-Nov-2008 -- ported last of the big pieces (3D viscous terms)
 */

/*-----------------------------------------------------------------*/
// Intel MPI requires mpi.h included BEFORE stdio.h
#ifdef _MPI
#   include <mpi.h>
#endif
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
#include "bc_extrapolate_out.hh"
#include "bc_user_defined.hh"
#include "bc_menter_correction.hh"
#include "exch2d.hh"
#include "visc.hh"
#include "visc3D.hh"
#ifdef _MPI
#   include "exch_mpi.hh"
#endif
#include "piston.hh"
#include "implicit.hh"
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
    int run_status = SUCCESS;
    char c, job_name[132], text_buf[132];
    char log_file_name[132], mpimap_file_name[132];

#   ifdef _MPI
    MPI_Status status;
    int node_name_len;
    char node_name[MPI_MAX_PROCESSOR_NAME];
#   endif

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
	{ "mpimap", 'm', POPT_ARG_STRING, NULL, 'm',
	  "use this specific MPI map of blocks to rank", 
	  "<mpimap_file>" },
	POPT_AUTOHELP
	POPT_TABLEEND
    };

    // Initialization
    //
    program_return_flag = SUCCESS; // We'll start out optimistically :)
#   ifdef _MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &(G.num_mpi_proc));
    MPI_Comm_rank(MPI_COMM_WORLD, &(G.my_mpi_rank));
    UNUSED_VARIABLE(status);
    G.mpi_parallel = 1;
#   else
    G.mpi_parallel = 0;
    G.num_mpi_proc = 0;
    G.my_mpi_rank = 0;
#   endif

    // Each process writes to a separate log file,
    // but all processes may write to stdout as well.
#   ifdef _MPI
    master = (G.my_mpi_rank == 0);
    sprintf(log_file_name, "e3mpi.%04d.log", G.my_mpi_rank);
    if (master) {
	printf("e3main: C++,MPI version.\n");
    }
    MPI_Get_processor_name(node_name, &node_name_len);
    printf("e3main: process %d of %d active on node %s\n", 
	   G.my_mpi_rank, G.num_mpi_proc, node_name );
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   else
    master = 1;
    sprintf(log_file_name, "e3shared.log");
    printf("e3main: C++,shared-memory version.\n");
#   ifdef _OPENMP
    printf("OpenMP version using %d thread(s).\n", omp_get_max_threads());
#   error "OpenMP version not functional..."    
#   endif
#   endif
    if ((G.logfile = fopen(log_file_name, "w")) == NULL) {
        printf( "\nCould not open %s; BAILING OUT\n", log_file_name );
        exit(FILE_ERROR);
    }
    // Configuration.
    //
    // printf("e3main: process command-line options...\n");
    optCon = poptGetContext(NULL, argc, (const char**)argv, optionsTable, 0);
    if (argc < 2) {
	poptPrintUsage(optCon, stderr, 0);
	goto Quit;
    }
    strcpy(job_name, "");
    strcpy(mpimap_file_name, "");
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
	case 'm':
	    strcpy(mpimap_file_name, poptGetOptArg(optCon));
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
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    if ( SUCCESS != read_config_parameters(G.base_file_name+".config", master) ) goto Quit; 
    // Read the time-step control parameters from a separate file.
    // These will also be read at the start of each time step.
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    if ( SUCCESS != read_control_parameters(G.base_file_name+".control", master, 1) ) goto Quit;
    // At this point, we know the number of blocks in the calculation.
    // Depending on whether we are running all blocks in the one process
    // or we are running a subset of blocks in this process, talking to
    // the other processes via MPI, we need to decide what blocks belong
    // to the current process.
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    if ( SUCCESS != assign_blocks_to_mpi_rank(mpimap_file_name, master) ) goto Quit;

    // Real work is delegated to functions that know what to do...
    //
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    if ( do_run_simulation == 1 ) {
	if ( master ) printf("Run simulation...\n");
	/* The simulation proper. */
	run_status = prepare_to_integrate(start_tindx);
	if (run_status != SUCCESS) goto Quit;
	if ( G.sequence_blocks ) {
	    run_status = integrate_blocks_in_sequence();
	    if (run_status != SUCCESS) goto Quit;
	} else {
	    run_status = integrate_in_time(-1.0);
	    if (run_status != SUCCESS) goto Quit;
	}
	finalize_simulation();
    } else {
	printf( "NOTHING DONE -- because you didn't ask...\n" );
    }

    // Finalization.
    //
    Quit: /* nop */;
    fclose(G.logfile);
    eilmer_finalize();
#   ifdef _MPI
    printf("e3main: end of process %d.\n", G.my_mpi_rank);
    if (run_status != SUCCESS) 
	MPI_Abort(MPI_COMM_WORLD, program_return_flag);
    else
	MPI_Finalize();
#   else
    printf("e3main: done.\n");
#   endif
    return program_return_flag;
} // end main

//------------------------------------------------------------------------------

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
    // We really do need this barrier.
    // Even if we are not the master node, we have to wait for the master
    // to check and (possibly) make the directory.
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


int prepare_to_integrate( int start_tindx )
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    Block *bdp;
    int iface;
    string filename, commandstring, jbstring, tindxstring, ifstring;
    char jbcstr[10], tindxcstr[10], ifcstr[10];
    FILE *fp;

    if ( G.npiston > 0 ) {
	filename = G.base_file_name;
	filename += ".piston0";
	if ((fp = fopen(filename.c_str(), "r")) == NULL) {
	    cerr << "Could not open " << filename << "; BAILING OUT" << endl;
	    return FILE_ERROR;
	}
	double temporary_time;
	for ( int jp = 0; jp < G.npiston; ++jp ) {
	    temporary_time = G.pistons[jp]->read_state( fp );
	    cout << "temporary_time=" << temporary_time 
		 << " Brendan, FIX-ME please. What are we doing with the piston code?" << endl;
	}
	fclose( fp );
    }

    // Allocate and initialise memory.
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
        if ( bdp->array_alloc(G.dimensions) != 0 ) exit( MEMORY_ERROR );
	bdp->bind_interfaces_to_cells(G.dimensions);
    }
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
    if ( allocate_send_and_receive_buffers() != 0 ) exit( MEMORY_ERROR );
#   endif

    // Read block grid and flow data; write history-file headers.
    // Note that the global simulation time is set by the last flow data read.
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    sprintf( tindxcstr, "t%04d", start_tindx);
    tindxstring = tindxcstr;
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
        if ( get_verbose_flag() ) printf( "----------------------------------\n" );
	sprintf( jbcstr, ".b%04d", bdp->id );
	jbstring = jbcstr;
	// Read grid from the tindx=0 files, always.
	filename = "grid/t0000/"+G.base_file_name+".grid"+jbstring+".t0000";
        if (bdp->read_grid(filename, G.dimensions, zip_files) != SUCCESS) {
	    return FAILURE;
	}
	// Read flow data from the specified tindx files.
	filename = "flow/"+tindxstring+"/"+G.base_file_name+".flow"+jbstring+"."+tindxstring;
        if (bdp->read_solution(filename, &(G.sim_time), G.dimensions, zip_files) != SUCCESS) {
	    return FAILURE;
	}
    }

    // History file header is only written for a fresh start.
    ensure_directory_is_present("hist"); // includes Barrier

    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
        if ( get_verbose_flag() ) printf( "----------------------------------\n" );
	sprintf( jbcstr, ".b%04d", bdp->id );
	jbstring = jbcstr;
	filename = "hist/" + G.base_file_name + ".hist"+jbstring;
	if ( access(filename.c_str(), F_OK) != 0 ) {
	    // History file does not yet exist; write header.
	    bdp->write_history( filename, G.sim_time, 1 );
	}
	// Read in heat-flux vectors if present
	for ( iface = NORTH; iface <= ((G.dimensions == 3)? BOTTOM : WEST); ++iface ) {
	    sprintf( ifcstr, ".s%04d", iface );
	    ifstring = ifcstr;
	    filename = "heat/"+tindxstring+"/"+G.base_file_name+".heat"+jbstring+ifstring+"."+tindxstring;
	    bdp->bcp[iface]->read_surface_heat_flux( filename, G.dimensions, zip_files );
	}
    }
    output_counter = start_tindx;
    filename = G.base_file_name; filename += ".times";
    if ( master ) {
	if ((G.timestampfile = fopen(filename.c_str(), "a")) == NULL) {
	    cerr << "Could not open " << filename << "; BAILING OUT" << endl;
	    return FILE_ERROR;
	}
	// The zeroth entry in the timestampfile has been written already by e3prep.py.
	// fprintf( G.timestampfile, "# count sim_time dt_global\n" );
	// fprintf( G.timestampfile, "%04d %e %e\n", output_counter, G.sim_time, G.dt_global );
	// fflush( G.timestampfile );
    }

#   ifdef _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#   endif

    // Prepare data within the primary finite-volume cells.
    // This includes both geometric data and flow state data.
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
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
    // FIX-ME generalize for handling several blocks per MPI process
    // FIX-ME Be careful that we don't break the usage for integrate_blocks_in_sequence
#   ifdef _MPI
    mpi_exchange_boundary_data(COPY_CELL_LENGTHS);
#   else
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
        exchange_shared_boundary_data(jb, COPY_CELL_LENGTHS);
    }
#   endif

#   if 0
    // Piston stuff not presently functional.
    int status;
    block_data** bda = (block_data**) calloc(G.nblock, sizeof(block_data*));
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bda[jb] = G.my_blocks[jb];
    }
    for (int jp = 0; jp < G.npiston; ++jp) {
	cout << "Attempting to set piston[" << jp << "] face values from block data... ";
	status = G.pistons[jp]->initialise_piston_faces(G, bda);

	if (status != SUCCESS) {
	    cout << "Failure!" << endl;
	    return FAILURE;
	}
	cout << "Success!" << endl;
    }
#   endif

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

    start = time(NULL); // start of wallclock timing
    return SUCCESS;
} // end prepare_to_integrate()

//------------------------------------------------------------------------

int call_udf( double t, int step, std::string udf_fn_name )
{
    lua_getglobal(L, udf_fn_name.c_str());  // function to be called
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, t); lua_setfield(L, -2, "t");
    lua_pushinteger(L, step); lua_setfield(L, -2, "step");
    int number_args = 1; // table of {t, step}
    int number_results = 0; // no results returned on the stack.
    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L, "error running user-defined function: %s\n",
			 lua_tostring(L, -1));
    }
    lua_settop(L, 0); // clear the stack
    return SUCCESS;
} // end call_udf()


/// \brief Add to the components of the source vector, Q, via a Lua udf.
/// This is done (occasionally) just after the for inviscid source vector calculation.
int udf_source_vector_for_cell( FV_Cell *cell, double t )
{
    // Call the user-defined function which returns a table 
    // of source term values.
    // These are added to the inviscid source terms 
    // that were computed earlier in the time step.

    int nsp = get_gas_model_ptr()->get_number_of_species();
    int nmodes = get_gas_model_ptr()->get_number_of_modes();

    lua_getglobal(L, "source_vector");  // Lua function to be called
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, t); lua_setfield(L, -2, "t");

    // Pack the interesting cell data in a table with named fields.
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, cell->pos.x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, cell->pos.y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, cell->pos.z); lua_setfield(L, -2, "z");
    lua_pushnumber(L, cell->volume); lua_setfield(L, -2, "vol");
    lua_pushnumber(L, cell->fs->gas->p); lua_setfield(L, -2, "p");
    lua_pushnumber(L, cell->fs->gas->rho); lua_setfield(L, -2, "rho"); 
    lua_pushnumber(L, cell->fs->vel.x); lua_setfield(L, -2, "u"); 
    lua_pushnumber(L, cell->fs->vel.y); lua_setfield(L, -2, "v");
    lua_pushnumber(L, cell->fs->vel.z); lua_setfield(L, -2, "w");
    lua_pushnumber(L, cell->fs->gas->a); lua_setfield(L, -2, "a");
    lua_newtable(L); // A table for the temperatures
    for ( int i = 0; i < nmodes; ++i ) {
	lua_pushinteger(L, i);
	lua_pushnumber(L, cell->fs->gas->T[i]);
	lua_settable(L, -3);
    }
    // At this point, the table of temperatures should be TOS.
    lua_setfield(L, -2, "T");
    lua_newtable(L); // Another table for the mass fractions
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_pushnumber(L, cell->fs->gas->massf[isp]);
	lua_settable(L, -3);
    }
    // At this point, the table of mass fractions should be TOS.
    lua_setfield(L, -2, "massf");
    
    // After all of this we should have ended up with the cell-data table at TOS 
    // with another table (containing t...} and function-name 
    // at pseudo-index locations -2 and -3 respectively.

    int number_args = 2; // table of {t...}, cell-data-table
    int number_results = 1; // one table of results returned on the stack.
    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L, "error running user source_terms function: %s\n",
			 lua_tostring(L, -1));
    }
    // Assume that there is a single table at the TOS
    // cout << "Get mass, momentum and energy source terms." << endl;
    lua_getfield(L, -1, "mass"); cell->Q->mass += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_x"); cell->Q->momentum.x += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_y"); cell->Q->momentum.y += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_z"); cell->Q->momentum.z += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "total_energy"); cell->Q->total_energy += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "romega"); cell->Q->omega += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "rtke"); cell->Q->tke += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "radiation"); cell->Q_rE_rad += lua_tonumber(L, -1); lua_pop(L, 1);
    cell->Q->total_energy += cell->Q_rE_rad;

    // cout << "Table of sources for mass fractions" << endl;
    lua_getfield(L, -1, "species"); // put mass-fraction table at TOS 
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_gettable(L, -2);
	cell->Q->massf[isp] += lua_tonumber(L, -1);
	lua_pop(L, 1); // remove the number to leave the table at TOS
    }
    lua_pop(L, 1); // remove species table from top-of-stack

    // cout << "Table of sources for individual energies" << endl;
    lua_getfield(L, -1, "energies"); // put energy table at TOS 
    for ( int i = 1; i < nmodes; ++i ) {
	lua_pushinteger(L, i);
	lua_gettable(L, -2);
	cell->Q->energies[i] += lua_tonumber(L, -1);
	lua_pop(L, 1); // remove the number to leave the table at TOS
    }
    lua_pop(L, 1); // remove energies table from top-of-stack

    lua_settop(L, 0); // clear the stack
    // cout << "End of udf_source_vector_for_cell()" << endl;
    return SUCCESS;
} // end udf_source_vector_for_cell()


//---------------------------------------------------------------------------

int integrate_blocks_in_sequence( void )
// This procedure integrates the blocks two-at-a-time in sequence.
//
// The idea is to approximate the space-marching approach of sm_3d
// and, hopefully, achieve a significant speed-up over integration
// of the full array of blocks.
//
// It is assumed that we have the restricted case 
// of all the blocks making a single line (west to east) with
// supersonic flow at the west face of block 0 and 
// extrapolate_out on the east face of the final block.
//
// The above assumptions still stands but the code has been modified
// such that the calcluation now does two active blocks at a time.
// The calculation initially does blocks 0 and 1, and then moves over
// one block for the next step calculating blocks 1 and 2 next.
// This process is then repeated for the full calculation

{
    global_data &G = *get_global_data_ptr();
    Block *bdp;
    double time_slice = G.max_time / (G.nblock - 1);
    BoundaryCondition *bcp_save;
    int status_flag = SUCCESS;

    // This procedure does not work with MPI jobs.
    if ( G.mpi_parallel ) {
	cerr << "Error, we have not implemented block-sequence integration with MPI." << endl;
	exit(NOT_IMPLEMENTED_ERROR);
    }

    // Initially deactivate all blocks
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	G.bd[jb].active = 0;
    }

    cout << "Integrate Block 0 and Block 1" << endl;

    // Start by setting up block 0
    bdp = &(G.bd[0]);
    bdp->active = 1;
    // Apply the assumed SupINBC to the west face and propogate across the block
    bdp->bcp[WEST]->apply_inviscid(0.0);
    bdp->propagate_data_west_to_east( G.dimensions );
    bdp->apply( &FV_Cell::encode_conserved, bdp->omegaz, "encode_conserved" );
    // Even though the following call appears redundant at this point,
    // fills in some gas properties such as Prandtl number that is
    // needed for both the cfd_check and the BLomax turbulence model.
    bdp->apply( &FV_Cell::decode_conserved, bdp->omegaz, "decode_conserved" );

    // Now set up block 1
    bdp = &(G.bd[1]);
    bdp->active = 1;
    // Save the original east boundary condition and apply the temporary
    // ExtrapolateOutBC for the calculation
    bcp_save = bdp->bcp[EAST];
    bdp->bcp[EAST] = new ExtrapolateOutBC(*bdp, EAST, 0);
    // Read in data from block 0 and propogate across the block
    exchange_shared_boundary_data( 1, COPY_FLOW_STATE);
    bdp->propagate_data_west_to_east( G.dimensions );
    bdp->apply( &FV_Cell::encode_conserved, bdp->omegaz, "encode_conserved" );
    bdp->apply( &FV_Cell::decode_conserved, bdp->omegaz, "decode_conserved" );

    // Integrate just the first two blocks in time, hopefully to steady state.
    G.bd[0].active = 1;
    G.bd[1].active = 1;
    integrate_in_time( time_slice );

    // The rest of the blocks.
    for ( int jb = 2; jb < (G.nblock); ++jb ) {
	// jb-2, jb-1, jb
	// jb-2 is the block to be deactivated, jb-1 has been iterated and now
	// becomes the left most block and jb is the new block to be iterated
	cout << "Integrate Block " << jb << endl;
	// Make the block jb-2 inactive.
	bdp = &(G.bd[jb-2]);
	bdp->active = 0;

	// block jb-1 - reinstate the previous boundary condition on east face
	// but leave the block active
	bdp = &(G.bd[jb-1]);
	delete bdp->bcp[EAST];
	bdp->bcp[EAST] = bcp_save;

	// Set up new block jb to be integrated
	bdp = &(G.bd[jb]);
	bdp->active = 1;
	if ( jb < G.nblock-1 ) {
	    // Cut off the east boundary of the current block 
	    // from the downstream blocks if there are any.
	    bcp_save = bdp->bcp[EAST];
	    bdp->bcp[EAST] = new ExtrapolateOutBC(*bdp, EAST, 0);
	}
	// Now copy the starting data into the WEST ghost cells
	// and propagate it across the current block.
	exchange_shared_boundary_data( jb, COPY_FLOW_STATE );
	bdp->propagate_data_west_to_east( G.dimensions );
	bdp->apply( &FV_Cell::encode_conserved, bdp->omegaz, "encode_conserved" );
	bdp->apply( &FV_Cell::decode_conserved, bdp->omegaz, "decode_conserved" );
	// Integrate just the two currently active blocks in time,
	// hopefully to steady state.
	G.bd[jb-1].active = 1;
	G.bd[jb].active = 1;
	integrate_in_time( jb*time_slice );
    }
    // Before leaving, we want all blocks active for output.
    for ( int jb = 0; jb < G.nblock; ++jb ) {
	G.bd[jb].active = 1;
    }
    return status_flag;
} // end integrate_blocks_in_sequence()

//---------------------------------------------------------------------------

int integrate_in_time( double target_time )
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    Block *bdp;
    int active_blocks;
    string jbstring, jsstring, tindxstring;
    char jbcstr[10], jscstr[10], tindxcstr[10];
    string filename, commandstring, foldername;
    std::vector<double> dt_record;
    double stopping_time;
    int finished_time_stepping;
    int viscous_terms_are_on;
    int cfl_result, do_cfl_check_now;
    int js, final_s;
#   ifdef _MPI
    int cfl_result_2;
#   endif
    int status_flag = SUCCESS;
    dt_record.resize(G.my_blocks.size()); // Just the block local to this process.

    if ( master ) {
	printf( "Integrate in time\n" );
	fflush(stdout);
    }
    if ( target_time <= 0.0 ) {
	stopping_time = G.max_time;
    } else {
	stopping_time = target_time;
    }
    // The next time for output...
    G.t_plot = G.sim_time + G.dt_plot;
    G.t_his = G.sim_time + G.dt_his;

    // Flags to indicate that the saved output is fresh.
    // On startup or restart, it is assumed to be so.
    output_just_written = 1;
    history_just_written = 1;
    av_output_just_written = 1;

    G.step = 0; // Global Iteration Count
    do_cfl_check_now = 0;
    if ( G.heat_time_stop == 0.0 ) {
	// We don't want heating at all.
	set_heat_factor( 0.0 );
    } else if ( G.sim_time >= G.heat_time_start && G.sim_time < G.heat_time_stop ) {
	// Apply full heat-source effects because both time limits are set
	// and we are within them.
	set_heat_factor( 1.0 );
    } else {
	// Looks like we want heating at some period in time but it is not now.
	set_heat_factor( 0.0 );
    }
    if ( get_viscous_flag() ) {
	// We have requested viscous effects but they may be delayed.
	if ( G.viscous_time_delay > 0.0 && G.sim_time < G.viscous_time_delay ) {
	    // We will initially turn-down viscous effects and
	    // only turn them up when the delay time is exceeded.
	    set_viscous_factor( 0.0 );
	    viscous_terms_are_on = 0;
	} else {
	    // No delay in applying full viscous effects.
	    set_viscous_factor( 1.0 );
	    viscous_terms_are_on = 1;
	}
    } else {
	// We haven't requested viscous effects at all.
	set_viscous_factor( 0.0 );
	viscous_terms_are_on = 0;
    }

    // Spatial filter may be applied occasionally.
    if ( G.do_filter ) {
	if ( G.sim_time > G.filter_tstart ) {
	    G.filter_next_time = G.sim_time + G.filter_dt;
	    if ( G.filter_next_time > G.filter_tend ) G.do_filter = 0;
	} else {
	    G.filter_next_time = G.filter_tstart;
	}
    }

    // Normally, we can terminate upon either reaching 
    // a maximum time or upon reaching a maximum iteration count.
    finished_time_stepping = (G.sim_time >= stopping_time || G.step >= G.max_step);

    //----------------------------------------------------------------
    //                 Top of main time-stepping loop
    //----------------------------------------------------------------
    while ( !finished_time_stepping ) {
	if ( (G.step/G.control_count)*G.control_count == G.step ) {
	    // Reparse the time-step control parameters as frequently as specified.
	    read_control_parameters(G.base_file_name+".control", master, 0);
	}
	// One of the control parameters is G.max_time, so we need to
	// make sure that the stopping_time is updated accordingly.
	if ( target_time <= 0.0 ) {
	    stopping_time = G.max_time;
	} else {
	    stopping_time = target_time;
	}

        // Continue taking steps so long as there is at least one active block.
#       ifdef _MPI
        // FIX-ME -- at the moment we assume all blocks are active
	// in the distributed-memory version.
	active_blocks = G.nblock;
#       else
        active_blocks = 0;
	for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
            if ( bdp->active == 1 ) ++active_blocks;
        }
#       endif
        if ( active_blocks == 0 ) {
            printf( "There are no active blocks at step %d\n", G.step );
	    status_flag = FAILURE;
            break;
        }

        // 0. Alter configuration setting if necessary.
	if ( get_viscous_flag() == 1 && viscous_terms_are_on == 0 && 
	     G.sim_time >= G.viscous_time_delay ) {
	    // We want to turn on the viscous effects only once (if requested)
	    // and, when doing so in the middle of a simulation, 
	    // reduce the time step to ensure that these new terms
	    // do not upset the stability of the calculation.
	    printf( "Turning on viscous effects.\n" );
	    viscous_terms_are_on = 1;
	    if ( G.step > 0 ) {
		printf( "Start softly with viscous effects.\n" );
		set_viscous_factor( 0.0 );  
		G.dt_global *= 0.2;
	    } else {
		printf( "Start with full viscous effects.\n" );
		set_viscous_factor( 1.0 );
	    }
	}
	if ( viscous_terms_are_on == 1 && get_viscous_factor() < 1.0 ) {
	    incr_viscous_factor( get_viscous_factor_increment() );
	    printf( "Increment viscous_factor to %f\n", get_viscous_factor() );
	    do_cfl_check_now = 1;
	}
	if ( G.heat_time_stop > 0.0 ) {
	    // We want heating at some time.
	    if ( G.sim_time >= G.heat_time_start && G.sim_time < G.heat_time_stop ) {
		// We are within the period of heating but
		// we may be part way through turning up the heat gradually.
		if ( get_heat_factor() < 1.0 ) {
		    incr_heat_factor( get_heat_factor_increment() );
		    printf( "Increment heat_factor to %f\n", get_heat_factor() );
		    do_cfl_check_now = 1;
		} 
	    } else {
		// We are outside the period of heating.
		set_heat_factor( 0.0 );
	    }
	}
	// 0.a. call out to user-defined function
	if ( G.udf_file.length() > 0 ) {
	    call_udf( G.sim_time, G.step, "at_timestep_start" );
	}
	// 0.b. We need to remember some flow information for later computing residuals.
	for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
	    bdp->init_residuals( G.dimensions );
	} // end for jb loop 


	// 1. Set the size of the time step.
	if ( G.step == 0 ) {
	    // When starting a new calculation,
	    // set the global time step to the initial value.
	    do_cfl_check_now = 0;
		// if we are using sequence_blocks we don't want to reset the dt_global for each block
		if ( G.sequence_blocks && G.dt_global != 0 ) {
		    /* do nothing i.e. keep dt_global from previous block */ ;
		} else { 
		    G.dt_global = G.dt_init;
		}
	} else if ( !G.fixed_time_step && 
		    (G.step/G.cfl_count)*G.cfl_count == G.step ) {
	    // Check occasionally 
	    do_cfl_check_now = 1;
	} // end if (G.step == 0 ...

	if ( do_cfl_check_now == 1 ) {
	    // Adjust the time step to be the minimum allowed
	    // for any active block. 
	    G.dt_allow = 1.0e6; /* outrageously large so that it gets replace immediately */
	    G.cfl_max = 0.0;
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active == 1 ) {
		    cfl_result = bdp->determine_time_step_size( G.cfl_target, G.dimensions );
		    if ( cfl_result != 0 ) {
			program_return_flag = DT_SEARCH_FAILED;
			status_flag = FAILURE;
			goto conclusion;
		    }
		}
	    } // end for jb loop
	    // If we arrive here, cfl_result will be zero, indicating that all local blocks 
	    // have successfully determined a suitable time step.
#           ifdef _MPI
	    MPI_Allreduce( &cfl_result, &cfl_result_2, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	    if ( cfl_result_2 != 0 ) {
		program_return_flag = DT_SEARCH_FAILED;
		status_flag = FAILURE;
		goto conclusion;
	    }
#           endif
	    // Get an overview of the allowable timestep.
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		dt_record[jb] = 0.0;
	    }
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active != 1 ) continue;
		if ( bdp->dt_allow < G.dt_allow ) G.dt_allow = bdp->dt_allow;
		dt_record[jb] = bdp->dt_allow;
		if ( bdp->cfl_max > G.cfl_max ) G.cfl_max = bdp->cfl_max;
		if ( bdp->cfl_min < G.cfl_min ) G.cfl_min = bdp->cfl_min;
	    }
#           ifdef _MPI
	    // Finding the minimum allowable time step is very important for stability.
	    MPI_Allreduce(MPI_IN_PLACE, &(G.dt_allow), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	    MPI_Allreduce(MPI_IN_PLACE, &(G.cfl_max), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	    MPI_Allreduce(MPI_IN_PLACE, &(G.cfl_min), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#           endif
	    // Change the actual time step, as needed.
	    if ( G.dt_allow <= G.dt_global ) {
		// If we need to reduce the time step, do it immediately.
		G.dt_global = G.dt_allow;
	    } else {
		// Make the transitions to larger time steps gentle.
		G.dt_global = MINIMUM(G.dt_global*1.5, G.dt_allow);
	    }
	    do_cfl_check_now = 0;  // we have done our check for now
	} // end if do_cfl_check_now 

        if ( (G.step > 5) && (G.cfl_max < G.cfl_tiny) ) {
            // Record the worst CFL
            G.cfl_tiny = G.cfl_max;
            G.time_tiny = G.sim_time;
        }
        if ( G.cfl_max > 100.0 ) {
            // If the CFL has gone crazy, bail out.
            printf( "\nCFL = %e: breaking main loop\n", G.cfl_max );
	    status_flag = FAILURE;
            break;
        }


        // 2. Attempt a time step.

	// 2a.
	// explicit or implicit update of the inviscid terms.
	int break_loop2 = 0;
	if ( get_implicit_flag() == 0 ) {
	    // explicit update of inviscid terms
	    break_loop2 = gasdynamic_inviscid_increment();
	} else if ( get_implicit_flag() == 1 ) {
	    // point implicit update of inviscid terms
	    break_loop2 = gasdynamic_point_implicit_inviscid_increment();
	} else if ( get_implicit_flag() == 2 ) {
	    // fully implicit update of inviscid terms
	    break_loop2 = gasdynamic_fully_implicit_inviscid_increment();
	}
	if ( break_loop2 ) {
	    printf("Breaking main loop:\n");
	    printf("    time step failed at inviscid increment.\n");
	    status_flag = FAILURE;
	    break;
	}
        
	// 2b. Piston step.
	//     Using the chosen timestep and calculated flow state of the blocks, 
	//     calculate and apply the effect of the piston(s) motion
#       if 0
	int status;
	block_data** bda = (block_data**) calloc(G.nblock, sizeof(block_data*));
	for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bda[jp] = G.my_blocks[jb];
	}

	// Set all cells active
	activate_cells(G, bda);
	for( int jp = 0; jp < G.npiston; ++jp ) {
	    // Ignore this projectile if it's inactive.
	    if( G.pistons[jp]->is_inactive() ) 
		continue;
	    
	    status = G.pistons[jp]->configure_cells(G, bda);
	    if (status == END_OF_BLOCK) {
		cout << "As indicated, one of the pistons has reached the\n"
		     << "end of a block with a free end.  The only sensible thing\n"
		     << "to do is to terminate the program and write out the solution.\n";
		status_flag = FAILURE;
		goto conclusion;
	    }
	    else if (status != SUCCESS) {
		cout << "Error trying to configure piston boundary cells.\n";
		cout << "Bailing out!" << endl;
		status_flag = FAILURE;
		goto conclusion;
	    }
	}
		
	for (int jp = 0; jp < G.npiston; ++jp) {
	    status = G.pistons[jp]->change_of_state_due_to_motion(G, bda, G.dt_global);
	    if( status == PISTON_FAILURE ) {
		cout << "Error trying to update piston motion.\n";
		cout << "Bailing out!" << endl;
		status_flag = FAILURE;
		goto conclusion;
	    }
	}

	// Set all cells active
	activate_cells(G, bda);
	for( int jp = 0; jp < G.npiston; ++jp ) {
	    // Ignore this projectile if it's inactive.
	    if( G.pistons[jp]->is_inactive() ) 
		continue;
	    
	    status = G.pistons[jp]->configure_cells(G, bda);
	    if (status == END_OF_BLOCK) {
		cout << "As indicated, one of the pistons has reached the\n"
		     << "end of a block with a free end.  The only sensible thing\n"
		     << "to do is to terminate the program and write out the solution.\n";
		status_flag = FAILURE;
		goto conclusion;
	    }
	    else if (status != SUCCESS) {
		cout << "Error trying to configure piston boundary cells.\n";
		cout << "Bailing out!" << endl;
		status_flag = FAILURE;
		goto conclusion;
	    }
	}
#       endif

	// 2c.
	if ( get_viscous_flag() == 1 ) {
	    // We now have the option of explicit or point implicit update
	    // of the viscous terms, thanks to Ojas.
	    int break_loop = 0;
	    if ( get_implicit_flag() == 0 ) {
	    	// explicit update of viscous terms
		break_loop = gasdynamic_viscous_increment();
	    }
	    else if ( get_implicit_flag() == 1 ) {
	    	// point implicit update of viscous terms
		break_loop = gasdynamic_point_implicit_viscous_increment();
	    }
	    else if ( get_implicit_flag() == 2 ) {
	    	// fully implicit update of viscous terms
	    	break_loop = gasdynamic_fully_implicit_viscous_increment();
	    }
	    if ( break_loop ) {
		printf("Breaking main loop:\n");
		printf("    time step failed at viscous increment.\n");
		status_flag = FAILURE;
		break;
	    }
#           if SEPARATE_UPDATE_FOR_K_OMEGA_SOURCE == 1
	    if ( get_k_omega_flag() == 1 ) {
		for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		    bdp = G.my_blocks[jb];
		    if ( bdp->active != 1 ) continue;
		    bdp->apply( &FV_Cell::update_k_omega_properties, G.dt_global, "k_omega_increment" );
		}
	    }
#           endif
	}

        // 2d. Chemistry step. 
	//     Allow finite-rate evolution of species due
        //     to chemical reactions
        if ( get_reacting_flag() == 1 && G.sim_time >= G.reaction_time_start ) {
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active != 1 ) continue;
		bdp->apply( &FV_Cell::chemical_increment, G.dt_global, "chemical_increment" );
	    }
	}

	// 2e. Thermal step.
	//     Allow finite-rate evolution of thermal energy
	//     due to transfer between thermal energy modes.
	if ( get_energy_exchange_flag() == 1 && G.sim_time >= G.reaction_time_start  ) {
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active != 1 ) continue;
		bdp->apply( &FV_Cell::thermal_increment, G.dt_global, "thermal_increment" );
	    }
	}


        // 3. Update the time record and (occasionally) print status.
        ++G.step;
        output_just_written = 0;
        history_just_written = 0;
	    av_output_just_written = 0;
        G.sim_time += G.dt_global;

        if ( ((G.step / G.print_count) * G.print_count == G.step) && master ) {
            // Print the current time-stepping status.
            now = time(NULL);
            printf("Step=%7d t=%10.3e dt=%10.3e %s\n",
		   G.step, G.sim_time, G.dt_global,
		   time_to_go(start, now, G.step, G.max_step, G.dt_global, G.sim_time, G.max_time) );
            fprintf(G.logfile, "Step=%7d t=%10.3e dt=%10.3e %s\n",
		   G.step, G.sim_time, G.dt_global,
		   time_to_go(start, now, G.step, G.max_step, G.dt_global, G.sim_time, G.max_time) );
            fprintf(G.logfile, "CFL_min = %e, CFL_max = %e, dt_allow = %e\n",
		   G.cfl_min, G.cfl_max, G.dt_allow );
            fprintf(G.logfile, "Smallest CFL_max so far = %e at t = %e\n",
		   G.cfl_tiny, G.time_tiny );
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
                if ( bdp->active != 1 ) continue;
                fprintf(G.logfile, " dt[%d]=%e", bdp->id, dt_record[jb] );
            }
	    for ( int jp = 0; jp < G.npiston; ++jp ) {
		fprintf(G.logfile, "%s\n", G.pistons[jp]->string_repr().c_str() );
	    }
	    if ( active_blocks == 1 ) {
		fprintf(G.logfile, "\nThere is %d active block.\n", active_blocks );
	    } else {
		fprintf(G.logfile, "\nThere are %d active blocks.\n", active_blocks );
	    }
	    fflush( stdout );
        } // end if


        // 4. (Occasionally) Write out an intermediate solution
        if ( (G.sim_time >= G.t_plot) && !output_just_written ) {
	    ++output_counter;
	    if ( master ) {
		fprintf( G.timestampfile, "%04d %e %e\n", output_counter, G.sim_time, G.dt_global );
		fflush( G.timestampfile );
	    }
	    sprintf( tindxcstr, "t%04d", output_counter ); tindxstring = tindxcstr;
	    foldername = "flow/"+tindxstring;
	    ensure_directory_is_present(foldername); // includes Barrier
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		sprintf( jbcstr, ".b%04d", bdp->id ); jbstring = jbcstr; 
		filename = foldername+"/"+ G.base_file_name+".flow"+jbstring+"."+tindxstring;
		bdp->write_solution(filename, G.sim_time, G.dimensions, zip_files);
	    }
	    // Compute, store and write heat-flux data, if viscous simulation
	    if ( get_viscous_flag() ) {
		foldername = "heat/"+tindxstring;
		ensure_directory_is_present(foldername); // includes Barrier
		// Loop over blocks
		for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		    bdp = G.my_blocks[jb];
		    sprintf( jbcstr, ".b%04d", bdp->id ); jbstring = jbcstr;
		    final_s = ((G.dimensions == 3)? BOTTOM : WEST);
		    // Loop over boundaries/surfaces
		    for ( js = NORTH; js <= final_s; ++js ) {
			sprintf( jscstr, ".s%04d", js ); jsstring = jscstr;
			filename = foldername+"/"+ G.base_file_name+".heat" \
			                +jbstring+jsstring+"."+tindxstring;
			bdp->bcp[js]->compute_surface_heat_flux();
			bdp->bcp[js]->write_surface_heat_flux(filename,G.sim_time);
			if ( zip_files == 1 ) do_system_cmd("gzip -f "+filename);
		    }
		}
	    }
            output_just_written = 1;
            G.t_plot += G.dt_plot;
        }

        if ( (G.sim_time >= G.t_his) && !history_just_written ) {
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		sprintf(jbcstr, ".b%04d", bdp->id); jbstring = jbcstr;
		filename = "hist/"+G.base_file_name+".hist"+jbstring;
                bdp->write_history( filename, G.sim_time );
		bdp->print_forces( G.logfile, G.sim_time, G.dimensions );
	    }
	    for ( int jp = 0; jp < G.npiston; ++jp ) {
		filename = G.base_file_name;
		filename += ".piston";
#               if 0
		G.pistons[jp]->write_state( fp, G.sim_time );
#               else
		printf( "Brendan, please fix piston stuff...\n" );
		exit( NOT_IMPLEMENTED_ERROR );
#               endif
	    }
            history_just_written = 1;
            G.t_his += G.dt_his;
        }

	// Velocity profile recording (Andrew Denman  19-June-2003)
	//    Limited to taking velocity profiles along a constant xline
	//    defined by x_record parameter in block block_record.
	//    Calls function write_profile()
	if (G.do_record == 1) {
	    if ( fmod(double(G.step), double(G.step_record)) == 0 && 
		 G.sim_time > G.t_record) {
		printf("\nWriting Profile\n\n");
		bdp = get_block_data_ptr(G.block_record);
#               if 0
		write_profile( bdp, G.x_record, &G.n_record, fp, &G.t_start, 0 );
#               else
                printf("*** We have tried to use a non-existant function write_profile()\n");
		exit( NOT_IMPLEMENTED_ERROR );
#               endif
		printf( "main G.n_record = %d\n", G.n_record );
		/* increment counter used for file names */
		G.n_record++;
	    }
	}
	
        // 5. For steady-state approach, check the residuals for mass and energy.
        if ( (G.step / G.print_count) * G.print_count == G.step ) {
            G.mass_residual = 0.0;
            G.energy_residual = 0.0;
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
                if ( bdp->active == 1 ) {
		    bdp->compute_residuals( G.dimensions );
		    fprintf( G.logfile, "RESIDUAL mass block %d max: %e at (%g,%g,%g)\n",
			     bdp->id, bdp->mass_residual, bdp->mass_residual_loc.x,
			     bdp->mass_residual_loc.y, bdp->mass_residual_loc.z );
		    fprintf( G.logfile, "RESIDUAL energy block %d max: %e at (%g,%g,%g)\n",
			     bdp->id, bdp->energy_residual, bdp->energy_residual_loc.x,
			     bdp->energy_residual_loc.y, bdp->energy_residual_loc.z );
		}
            }
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
                if ( bdp->active == 1 ) {
		    if ( bdp->mass_residual > G.mass_residual ) {
			G.mass_residual = bdp->mass_residual;
		    }
		    if ( bdp->energy_residual > G.energy_residual ) {
			G.energy_residual = bdp->energy_residual;
		    }
		}
            }
#           ifdef _MPI
	    MPI_Allreduce(MPI_IN_PLACE, &(G.mass_residual), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	    MPI_Allreduce(MPI_IN_PLACE, &(G.energy_residual), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#           endif
            fprintf( G.logfile, "RESIDUAL mass global max: %e step %d time %g\n",
		     G.mass_residual, G.step, G.sim_time );
            fprintf( G.logfile, "RESIDUAL energy global max: %e step %d time %g\n",
		     G.energy_residual, G.step, G.sim_time );
	    fflush( G.logfile );
        }

	// 6. Spatial filter may be applied occasionally.
	if ( G.do_filter && G.sim_time > G.filter_next_time ) {
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active != 1 ) continue;
		bdp->apply_spatial_filter( G.filter_alpha, G.filter_npass, G.dimensions );
		bdp->apply( &FV_Cell::encode_conserved, bdp->omegaz, "encode_conserved" );
	    }
	    G.filter_next_time = G.sim_time + G.filter_dt;
	    if ( G.filter_next_time > G.filter_tend ) G.do_filter = 0;
	}

	// 7. Call out to user-defined function.
	if ( G.udf_file.length() > 0 ) {
	    call_udf( G.sim_time, G.step, "at_timestep_end" );
	}


        // 8. Loop termination criteria:
        //    (1) reaching a maximum simulation time or target time
        //    (2) reaching a maximum number of steps
        //    (3) finding that the "halt_now" parameter has been set 
	//        in the control-parameter file.
        //        This provides a semi-interactive way to terminate the 
        //        simulation and save the data.
	//    (4) Exceeding a maximum number of wall-clock seconds.
	//    (5) Exceeding an allowable delta(f_rad) / f_rad_org factor
	//    Note that the max_time and max_step control parameters can also
	//    be found in the control-parameter file (which may be edited
	//    while the code is running).
        if ( G.sim_time >= stopping_time ) {
            finished_time_stepping = 1;
            if ( master ) printf( "Integration stopped: reached maximum simulation time.\n" );
        }
        if ( G.step >= G.max_step ) {
            finished_time_stepping = 1;
            if ( master ) printf( "Integration stopped: reached maximum number of steps.\n" );
        }
        if ( G.halt_now == 1 ) {
            finished_time_stepping = 1;
            if ( master ) printf( "Integration stopped: Halt set in control file.\n" );
        }
	now = time(NULL);
	if ( max_wall_clock > 0 && ( (int)(now - start) > max_wall_clock ) ) {
            finished_time_stepping = 1;
            if ( master ) printf( "Integration stopped: reached maximum wall-clock time.\n" );
	}
#       if CHECK_RADIATION_SCALING
	if ( get_radiation_flag() ) {
	    if ( check_radiation_scaling() ) {
	    	finished_time_stepping = 1;
	    	if ( master ) printf( "Integration stopped: radiation source term needs updating.\n" );
	    }
	}
#       endif
	    	
    } // end while
    //----------------------------------------------------------------
    //                Bottom of main time-stepping loop
    //----------------------------------------------------------------

 conclusion:
    dt_record.resize(0);
    return status_flag;
} // end integrate_in_time()


int finalize_simulation( void )
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    Block *bdp;
    string filename, commandstring, foldername, jbstring, jsstring;
    char jbcstr[10],jscstr[10];
    int js, final_s;

    // Write out the final solution even if it has just been written
    // as part of the main time-stepping loop.
    output_counter = 9999;
    if ( master ) {
	fprintf( G.timestampfile, "%04d %e %e\n", output_counter, G.sim_time, G.dt_global );
	fflush( G.timestampfile );
    }
    foldername = "flow/t9999";
    ensure_directory_is_present(foldername); // includes Barrier
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	sprintf( jbcstr, ".b%04d", bdp->id ); jbstring = jbcstr;
	filename = foldername+"/"+G.base_file_name+".flow"+jbstring+".t9999";
	bdp->write_solution(filename, G.sim_time, G.dimensions, zip_files);
    }
    // Compute, store and write heat-flux data, if viscous simulation
    if ( get_viscous_flag() ) {
	foldername = "heat/t9999";
	ensure_directory_is_present(foldername); // includes Barrier
	// Loop over blocks
	for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
	    sprintf( jbcstr, ".b%04d", bdp->id ); jbstring = jbcstr;
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
	for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
	    sprintf( jbcstr, ".b%04d", bdp->id ); jbstring = jbcstr;
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
    if ( master ) printf( "\nTotal number of steps = %d\n", G.step );

    filename = G.base_file_name; filename += ".finish";
    if ( master ) {
	write_finishing_data( &G, filename );
	fclose( G.timestampfile );
    }
#   ifdef _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#   endif

    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	bdp->array_cleanup(G.dimensions);
    }
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
    if ( delete_send_and_receive_buffers() != 0 ) exit( MEMORY_ERROR );
#   endif
    return SUCCESS;
} // end finalize_simulation()

//------------------------------------------------------------------------

int gasdynamic_inviscid_increment( void )
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    Block *bdp;
    int most_bad_cells;
    int attempt_number, step_failed;

    // Record the current values of the conserved variables
    // in preparation for applying the predictor and corrector
    // stages of the time step.
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
        if ( bdp->active != 1 ) continue;
	bdp->apply( &FV_Cell::record_conserved, "record_conserved" );
    }

    attempt_number = 0;
    do {
#       ifdef _MPI
        // Before we try to exchange data, everyone's data should be up-to-date.
	MPI_Barrier( MPI_COMM_WORLD );
#       endif
	++attempt_number;
	step_failed = 0;

	//  Predictor Stage for gas-dynamics
#       ifdef _MPI
	mpi_exchange_boundary_data(COPY_FLOW_STATE);
#       else
	for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
	    if ( bdp->active != 1 ) continue;
	    exchange_shared_boundary_data( jb, COPY_FLOW_STATE );
	}
#       endif
	for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
	    if ( bdp->active != 1 ) continue;
	    apply_inviscid_bc( *bdp, G.sim_time, G.dimensions );
	}

	if ( get_flux_calculator() == FLUX_ADAPTIVE ) {
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active != 1 ) continue;
		bdp->detect_shock_points( G.dimensions );
	    }
	}
	
	// Non-local radiation transport needs to be performed a-priori for parallelization.
	// Note that Q_rad is not re-evaluated for corrector step.
	if ( get_radiation_flag() ) {
	    RadiationTransportModel * rtm = get_radiation_transport_model_ptr();
	    global_data &G = *get_global_data_ptr();
	    Block * bdp;
	    
	    // Determine if a scaled or complete radiation call is required
	    int ruf = get_radiation_update_frequency();
	    if ( ( ruf == 0 ) || ( ( G.step / ruf ) * ruf != G.step ) ) {
		// rescale
#		ifdef _OPENMP
#		pragma omp parallel for private(jb) schedule(runtime)
#		endif
		for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		    bdp = G.my_blocks[jb];
		    if ( bdp->active != 1 ) continue;
		    bdp->apply( &FV_Cell::rescale_Q_rE_rad, "rescale-Q_rE_rad" );
		}
	    }
	    else {
		// recompute
		rtm->compute_Q_rad_for_flowfield();
		// store the radiation scaling parameters for each cell
#		ifdef _OPENMP
#		pragma omp parallel for private(jb) schedule(runtime)
#		endif
		for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		    bdp = G.my_blocks[jb];
		    if ( bdp->active != 1 ) continue;
		    bdp->apply( &FV_Cell::store_rad_scaling_params, "store-rad-scaling-params" );
		}
	    }
	}

	for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
	    if ( bdp->active != 1 ) continue;
	    bdp->inviscid_flux( G.dimensions );
	    bdp->apply( &FV_Cell::inviscid_source_vector, bdp->omegaz, "inviscid-source-vector" );
	    if ( G.udf_source_vector_flag == 1 ) {
		bdp->apply( udf_source_vector_for_cell, G.sim_time, "udf-source-vector" );
	    }
	    bdp->apply( &FV_Cell::time_derivatives, 0, G.dimensions, "time-derivatives-level-0" );
	    bdp->apply( &FV_Cell::predictor_update, G.dt_global, "predictor-step" );
	    bdp->apply( &FV_Cell::decode_conserved, bdp->omegaz, "decode-conserved" );
#           define WILSON_OMEGA_FILTER 0
#           if WILSON_OMEGA_FILTER == 1
            if ( get_k_omega_flag() ) apply_wilson_omega_correction( *bdp );
#           endif

	} // end of for jb...

	// Corrector Stage
	if ( get_Torder_flag() == 2 ) {
#           ifdef _MPI
	    MPI_Barrier( MPI_COMM_WORLD );
	    mpi_exchange_boundary_data(COPY_FLOW_STATE);
#           else
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active != 1 ) continue;
		exchange_shared_boundary_data( jb, COPY_FLOW_STATE );
	    }
#           endif
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active != 1 ) continue;
		apply_inviscid_bc( *bdp, G.sim_time, G.dimensions );
	    }
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active != 1 ) continue;
		bdp->inviscid_flux( G.dimensions );
		bdp->apply( &FV_Cell::inviscid_source_vector, bdp->omegaz, "inviscid-source-vector" );
		if ( G.udf_source_vector_flag == 1 ) {
		    bdp->apply( udf_source_vector_for_cell, G.sim_time, "udf-source-vector" );
		}
		bdp->apply( &FV_Cell::time_derivatives, 1, G.dimensions, "time_derivatives-level-1" );
		bdp->apply( &FV_Cell::corrector_update, G.dt_global, "corrector-step" );
		bdp->apply( &FV_Cell::decode_conserved, bdp->omegaz, "decode-conserved" );
#               if WILSON_OMEGA_FILTER == 1
                if ( get_k_omega_flag() ) apply_wilson_omega_correction( *bdp );
#               endif
	    } // end for jb loop
	} // end if (corrector stage)

	// 2d. Check the record of bad cells and if any cells are bad, 
	//     fail this attempt at taking a step,
	//     set everything back to the initial state and
	//     reduce the time step for the next attempt
	most_bad_cells = do_bad_cell_count();
	if ( ADJUST_INVALID_CELL_DATA == 0 && most_bad_cells > 0 ) {
	    step_failed = 1;
	}
#       ifdef _MPI
	MPI_Allreduce(MPI_IN_PLACE, &(step_failed), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#       endif
	if ( step_failed ) {
	    G.dt_global = G.dt_reduction_factor * G.dt_global;
	    printf("Attempt %d failed: reducing dt to %e.\n", attempt_number, G.dt_global);
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active != 1 ) continue;
		bdp->apply( &FV_Cell::restore_conserved, "restore_conserved" );
		bdp->apply( &FV_Cell::decode_conserved, bdp->omegaz, "decode_conserved" );
#               if WILSON_OMEGA_FILTER == 1
                if ( get_k_omega_flag() ) apply_wilson_omega_correction( *bdp );
#               endif
	    }
	}

    } while (attempt_number < 3 && step_failed == 1);

    return step_failed;
} // end gasdynamic_inviscid_increment()


int gasdynamic_viscous_increment( void )
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    Block *bdp;
    // Record the current values of the conserved variables
    // in preparation for applying the predictor and corrector
    // stages of the time step.
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
        if ( bdp->active != 1 ) continue;
	bdp->apply( &FV_Cell::record_conserved, "record-conserved" );
	// cout << "After record_conserved in viscous_increment. jb=" << jb << endl;
	// cout << "mu_t=" << bdp->get_ifi(bdp->imin,bdp->jmin)->mu_t << endl;
	// cout << "k_t=" << bdp->get_ifi(bdp->imin,bdp->jmin)->k_t << endl;
	// print_data_for_cell( bdp->get_cell(bdp->imin,bdp->jmin), 1 );
	// print_data_for_interface( bdp->get_ifi(bdp->imin,bdp->jmin), 1 );
    }
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	if ( bdp->active != 1 ) continue;
	bdp->clear_fluxes_of_conserved_quantities( G.dimensions );
	apply_viscous_bc( *bdp, G.sim_time, G.dimensions );
	if ( get_k_omega_flag() ) apply_menter_boundary_correction( *bdp );
#       define WILSON_OMEGA_FILTER 0
#       if WILSON_OMEGA_FILTER == 1
        if ( get_k_omega_flag() ) apply_wilson_omega_correction( *bdp );
#       endif
	if ( G.dimensions == 2 ) {
	    viscous_derivatives_2D( bdp );
	} else {
	    viscous_derivatives_3D( bdp );
	}
#       define DEBUG_KOMEGA_FLAG 0
#       if DEBUG_KOMEGA_FLAG == 1
        printf("Revised debug mode .... W.C. and P.J.  12-Feb-2009\n");
        if ( jb == 0 ) {
            double dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
            double dtkedx, dtkedy, dtkedz;
            double domegadx, domegady, domegadz;
            int i = 10;
            int j = 5;
            int k;
            FILE * debug_file;
            // for 3D couette debug exercise,
	    // k == 2 is right on the bottom boundary
	    // k == 11 is on top boundary
            if ( G.dimensions == 2 ) k = 0; else k = 11;
            //
            FV_Cell *cell = bdp->get_cell(i,j,k);
            bdp = get_block_data_ptr(jb);
            //
            if ( G.dimensions == 2 ) {
                dudx = 0.25 * (cell->vtx[0]->dudx + cell->vtx[1]->dudx +
                               cell->vtx[2]->dudx + cell->vtx[3]->dudx);
                dudy = 0.25 * (cell->vtx[0]->dudy + cell->vtx[1]->dudy +
                               cell->vtx[2]->dudy + cell->vtx[3]->dudy);
                dudz = 0.25 * (cell->vtx[0]->dudz + cell->vtx[1]->dudz +
                               cell->vtx[2]->dudz + cell->vtx[3]->dudz);
                dvdx = 0.25 * (cell->vtx[0]->dvdx + cell->vtx[1]->dvdx +
                               cell->vtx[2]->dvdx + cell->vtx[3]->dvdx);
                dvdy = 0.25 * (cell->vtx[0]->dvdy + cell->vtx[1]->dvdy +
                               cell->vtx[2]->dvdy + cell->vtx[3]->dvdy);
                dvdz = 0.25 * (cell->vtx[0]->dvdz + cell->vtx[1]->dvdz +
                               cell->vtx[2]->dvdz + cell->vtx[3]->dvdz);
                dwdx = 0.25 * (cell->vtx[0]->dwdx + cell->vtx[1]->dwdx +
                               cell->vtx[2]->dwdx + cell->vtx[3]->dwdx);
                dwdy = 0.25 * (cell->vtx[0]->dwdy + cell->vtx[1]->dwdy +
                               cell->vtx[2]->dwdy + cell->vtx[3]->dwdy);
                dwdz = 0.25 * (cell->vtx[0]->dwdz + cell->vtx[1]->dwdz +
                               cell->vtx[2]->dwdz + cell->vtx[3]->dwdz);
                dtkedx = 0.25 * (cell->vtx[0]->dtkedx + cell->vtx[1]->dtkedx +
                                 cell->vtx[2]->dtkedx + cell->vtx[3]->dtkedx);
                dtkedy = 0.25 * (cell->vtx[0]->dtkedy + cell->vtx[1]->dtkedy +
                                 cell->vtx[2]->dtkedy + cell->vtx[3]->dtkedy);
                dtkedz = 0.25 * (cell->vtx[0]->dtkedz + cell->vtx[1]->dtkedz +
                                 cell->vtx[2]->dtkedz + cell->vtx[3]->dtkedz);
                domegadx = 0.25 * (cell->vtx[0]->domegadx + cell->vtx[1]->domegadx +
                                   cell->vtx[2]->domegadx + cell->vtx[3]->domegadx);
                domegady = 0.25 * (cell->vtx[0]->domegady + cell->vtx[1]->domegady +
                                   cell->vtx[2]->domegady + cell->vtx[3]->domegady);
                domegadz = 0.25 * (cell->vtx[0]->domegadz + cell->vtx[1]->domegadz +
                                   cell->vtx[2]->domegadz + cell->vtx[3]->domegadz);
            } else {
                dudx = 0.125 * (cell->vtx[0]->dudx + cell->vtx[1]->dudx +
                                cell->vtx[2]->dudx + cell->vtx[3]->dudx +
                                cell->vtx[4]->dudx + cell->vtx[5]->dudx +
                                cell->vtx[6]->dudx + cell->vtx[7]->dudx);
                dudy = 0.125 * (cell->vtx[0]->dudy + cell->vtx[1]->dudy +
                                cell->vtx[2]->dudy + cell->vtx[3]->dudy +
                                cell->vtx[4]->dudy + cell->vtx[5]->dudy +
                                cell->vtx[6]->dudy + cell->vtx[7]->dudy);
                dudz = 0.125 * (cell->vtx[0]->dudz + cell->vtx[1]->dudz +
                                cell->vtx[2]->dudz + cell->vtx[3]->dudz +
                                cell->vtx[4]->dudz + cell->vtx[5]->dudz +
                                cell->vtx[6]->dudz + cell->vtx[7]->dudz);
                dvdx = 0.125 * (cell->vtx[0]->dvdx + cell->vtx[1]->dvdx +
                                cell->vtx[2]->dvdx + cell->vtx[3]->dvdx +
                                cell->vtx[4]->dvdx + cell->vtx[5]->dvdx +
                                cell->vtx[6]->dvdx + cell->vtx[7]->dvdx);
                dvdy = 0.125 * (cell->vtx[0]->dvdy + cell->vtx[1]->dvdy +
                                cell->vtx[2]->dvdy + cell->vtx[3]->dvdy +
                                cell->vtx[4]->dvdy + cell->vtx[5]->dvdy +
                                cell->vtx[6]->dvdy + cell->vtx[7]->dvdy);
                dvdz = 0.125 * (cell->vtx[0]->dvdz + cell->vtx[1]->dvdz + 
                                cell->vtx[2]->dvdz + cell->vtx[3]->dvdz +
                                cell->vtx[4]->dvdz + cell->vtx[5]->dvdz +
                                cell->vtx[6]->dvdz + cell->vtx[7]->dvdz);
                dwdx = 0.125 * (cell->vtx[0]->dwdx + cell->vtx[1]->dwdx +
                                cell->vtx[2]->dwdx + cell->vtx[3]->dwdx +
                                cell->vtx[4]->dwdx + cell->vtx[5]->dwdx +
                                cell->vtx[6]->dwdx + cell->vtx[7]->dwdx);
                dwdy = 0.125 * (cell->vtx[0]->dwdy + cell->vtx[1]->dwdy +
                                cell->vtx[2]->dwdy + cell->vtx[3]->dwdy +
                                cell->vtx[4]->dwdy + cell->vtx[5]->dwdy +
                                cell->vtx[6]->dwdy + cell->vtx[7]->dwdy);
                dwdz = 0.125 * (cell->vtx[0]->dwdz + cell->vtx[1]->dwdz +
                                cell->vtx[2]->dwdz + cell->vtx[3]->dwdz +
                                cell->vtx[4]->dwdz + cell->vtx[5]->dwdz +
                                cell->vtx[6]->dwdz + cell->vtx[7]->dwdz);
                dtkedx = 0.125 * (cell->vtx[0]->dtkedx + cell->vtx[1]->dtkedx +
                                  cell->vtx[2]->dtkedx + cell->vtx[3]->dtkedx +
                                  cell->vtx[4]->dtkedx + cell->vtx[5]->dtkedx +
                                  cell->vtx[6]->dtkedx + cell->vtx[7]->dtkedx);
                dtkedy = 0.125 * (cell->vtx[0]->dtkedy + cell->vtx[1]->dtkedy +
                                  cell->vtx[2]->dtkedy + cell->vtx[3]->dtkedy +
                                  cell->vtx[4]->dtkedy + cell->vtx[5]->dtkedy +
                                  cell->vtx[6]->dtkedy + cell->vtx[7]->dtkedy);
                dtkedz = 0.125 * (cell->vtx[0]->dtkedz + cell->vtx[1]->dtkedz +
                                  cell->vtx[2]->dtkedz + cell->vtx[3]->dtkedz +
                                  cell->vtx[4]->dtkedz + cell->vtx[5]->dtkedz +
                                  cell->vtx[6]->dtkedz + cell->vtx[7]->dtkedz);
                domegadx = 0.125 * (cell->vtx[0]->domegadx + cell->vtx[1]->domegadx +
                                    cell->vtx[2]->domegadx + cell->vtx[3]->domegadx +
                                    cell->vtx[4]->domegadx + cell->vtx[5]->domegadx +
                                    cell->vtx[6]->domegadx + cell->vtx[7]->domegadx);
                domegady = 0.125 * (cell->vtx[0]->domegady + cell->vtx[1]->domegady +
                                    cell->vtx[2]->domegady + cell->vtx[3]->domegady +
                                    cell->vtx[4]->domegady + cell->vtx[5]->domegady +
                                    cell->vtx[6]->domegady + cell->vtx[7]->domegady);
                domegadz = 0.125 * (cell->vtx[0]->domegadz + cell->vtx[1]->domegadz +
                                    cell->vtx[2]->domegadz + cell->vtx[3]->domegadz +
                                    cell->vtx[4]->domegadz + cell->vtx[5]->domegadz +
                                    cell->vtx[6]->domegadz + cell->vtx[7]->domegadz);
            }    
            // bdp->get_cell(i,j,k)
            // Open debug.dat for writing data
            debug_file = fopen("debug.dat","a");
            if ( G.step == 0 ) {
                fprintf( debug_file, "t = %12.5e\n\n", G.sim_time);
		fprintf( debug_file, "i=%d, j=%d, k=%d\n", i, j, k);
                fprintf( debug_file, "pos.x    : %12.5e\n", cell->pos.x );
                fprintf( debug_file, "pos.y    : %12.5e\n", cell->pos.y );
                fprintf( debug_file, "pos.z    : %12.5e\n\n", cell->pos.z );
                fprintf( debug_file, "dudx     : %12.5e\n", dudx );
                fprintf( debug_file, "dudy     : %12.5e\n", dudy );
                fprintf( debug_file, "dudz     : %12.5e\n\n", dudz );
                fprintf( debug_file, "dvdx     : %12.5e\n", dvdx );
                fprintf( debug_file, "dvdy     : %12.5e\n", dvdy );
                fprintf( debug_file, "dvdz     : %12.5e\n\n", dvdz );
                fprintf( debug_file, "dwdx     : %12.5e\n", dwdx );
                fprintf( debug_file, "dwdy     : %12.5e\n", dwdy );
                fprintf( debug_file, "dwdz     : %12.5e\n\n", dwdz );
                fprintf( debug_file, "dtkedx   : %12.5e\n", dtkedx );
                fprintf( debug_file, "dtkedy   : %12.5e\n", dtkedy );
                fprintf( debug_file, "dtkedz   : %12.5e\n\n", dtkedz );
                fprintf( debug_file, "domegadx : %12.5e\n", domegadx );
                fprintf( debug_file, "domegady : %12.5e\n", domegady );
                fprintf( debug_file, "domegadz : %12.5e\n\n", domegadz );
		for ( int ii=0; ii<8; ++ii ) {
		    FV_Vertex* vtx = cell->vtx[ii];
		    fprintf( debug_file, "vtx[%d].pos.x   : %12.5e   pos.y   : %12.5e   pos.z   : %12.5e   area: %12.5e  volume: %12.5e\n", 
			     ii, vtx->pos.x, vtx->pos.y, vtx->pos.z, vtx->area, vtx->volume);
		    fprintf( debug_file, "vtx[%d].dudx    : %12.5e   dudy    : %12.5e   dudz    : %12.5e\n", 
			     ii, vtx->dudx, vtx->dudy, vtx->dudz);
		    fprintf( debug_file, "vtx[%d].dvdx    : %12.5e   dvdy    : %12.5e   dvdz    : %12.5e\n", 
			     ii, vtx->dvdx, vtx->dvdy, vtx->dvdz);
		    fprintf( debug_file, "vtx[%d].dwdx    : %12.5e   dwdy    : %12.5e   dwdz    : %12.5e\n", 
			     ii, vtx->dwdx, vtx->dwdy, vtx->dwdz);
		    fprintf( debug_file, "vtx[%d].dtkedx  : %12.5e   dtkedy  : %12.5e   dtkedz  : %12.5e\n", 
			     ii, vtx->dtkedx, vtx->dtkedy, vtx->dtkedz);
		    fprintf( debug_file, "vtx[%d].domegadx: %12.5e   domegady: %12.5e   domegadz: %12.5e\n\n", 
			     ii, vtx->domegadx, vtx->domegady, vtx->domegadz);
		}
		fprintf( debug_file, "--------------------------------\n");
		FV_Interface *iface = bdp->get_sifi(i,j,k);
		fprintf( debug_file, "sifi.pos=(%g,%g,%g) area=%g n=(%g,%g,%g)\n", 
			 iface->pos.x, iface->pos.y, iface->pos.z,
			 iface->area, iface->n.x, iface->n.y, iface->n.z);
		iface = bdp->get_sifj(i,j,k);
		fprintf( debug_file, "sifj.pos=(%g,%g,%g) area=%g n=(%g,%g,%g)\n", 
			 iface->pos.x, iface->pos.y, iface->pos.z,
			 iface->area, iface->n.x, iface->n.y, iface->n.z);
		iface = bdp->get_sifk(i,j,k);
		fprintf( debug_file, "sifk.pos=(%g,%g,%g) area=%g n=(%g,%g,%g)\n", 
			 iface->pos.x, iface->pos.y, iface->pos.z,
			 iface->area, iface->n.x, iface->n.y, iface->n.z);
		fprintf( debug_file, "--------------------------------\n");
            }
            fclose( debug_file );
        }
#       endif
	estimate_turbulence_viscosity( &G, bdp );
	if ( G.dimensions == 2 ) {
	    viscous_flux_2D( bdp );
	} else {
	    viscous_flux_3D( bdp );
	}
	bdp->apply( &FV_Cell::viscous_source_vector, "viscous-source-vector" );
	bdp->apply( &FV_Cell::time_derivatives, 0, G.dimensions, "time-derivatives-level-0,viscous" );
	bdp->apply( &FV_Cell::predictor_update, G.dt_global, "predictor-step,viscous" );
	bdp->apply( &FV_Cell::decode_conserved, bdp->omegaz, "decode-conserved,viscous" );
#       if WILSON_OMEGA_FILTER == 1
        if ( get_k_omega_flag() ) apply_wilson_omega_correction( *bdp );
#       endif
    } // end of for jb...
    return 0;
} // int gasdynamic_viscous_increment()


int do_bad_cell_count( void )
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    Block *bdp;
    int bad_cell_count, most_bad_cells;

    most_bad_cells = 0;
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	bad_cell_count = bdp->count_invalid_cells( G.dimensions );
	if ( bad_cell_count > most_bad_cells ) {
	    most_bad_cells = bad_cell_count;
	} 
	if ( bad_cell_count > G.max_invalid_cells ) {
	    printf( "   Too many bad cells (i.e. %d > %d) in block[%d].\n", 
		    bad_cell_count, G.max_invalid_cells, jb );
	}
    } // end for jb loop
#   ifdef _MPI
    MPI_Allreduce(MPI_IN_PLACE, &most_bad_cells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#   endif
    return most_bad_cells;
} // end do_bad_cell_count()


/// \brief Write out simulation data at end of simulation, such as: no. steps, final time.
int write_finishing_data( global_data *G, std::string filename )
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

/// \brief Check that the maximum delta(f_rad) / f_rad_org value is not exceedingly large
int check_radiation_scaling( void )
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    Block *bdp;
    double f_max, block_f_max, global_f_max = 0.0;
    FV_Cell *cellp;
    
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
    	block_f_max = 0.0;
	for ( int k = bdp->kmin; k <= bdp->kmax; ++k ) {
	    for ( int j = bdp->jmin; j <= bdp->jmax; ++j ) {
		for ( int i = bdp->imin; i <= bdp->imax; ++i ) {
		    cellp = bdp->get_cell(i,j,k);
		    f_max = cellp->rad_scaling_ratio();
		    if ( f_max > block_f_max) block_f_max = f_max;
		}
	    }
	}
	if ( block_f_max > global_f_max ) global_f_max = block_f_max;
    }
#   ifdef _MPI
    MPI_Allreduce(MPI_IN_PLACE, &global_f_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#   endif

    // printf("Global radiation f_max = %f\n", global_f_max);
    if ( global_f_max > 1000000.0 )
    	return FAILURE;
    else
    	return SUCCESS;
}
