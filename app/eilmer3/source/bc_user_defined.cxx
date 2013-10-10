// bc_user_defined.cxx

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <valarray>

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "bc_user_defined.hh"
#include "bc_menter_correction.hh"
#include "lua_service_for_e3.hh"

//----------------------------------------------------------------------------

UserDefinedBC::UserDefinedBC(Block *bdp, int which_boundary,
			     const std::string filename_,
			     bool is_wall, 
			     bool sets_conv_flux, bool sets_visc_flux)
    : BoundaryCondition(bdp, which_boundary, USER_DEFINED),
      filename(filename_)
{
    is_wall_flag = is_wall;
    sets_conv_flux_flag = sets_conv_flux;
    if ( sets_conv_flux ) ghost_cell_data_available = false;
    sets_visc_flux_flag = sets_visc_flux;
    start_interpreter();
} // end constructor

UserDefinedBC::UserDefinedBC(const UserDefinedBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      filename(bc.filename) 
{
    is_wall_flag = bc.is_wall_flag;
    sets_conv_flux_flag = bc.sets_conv_flux_flag;
    ghost_cell_data_available = bc.ghost_cell_data_available;
    sets_visc_flux_flag = bc.sets_visc_flux_flag;
    start_interpreter();
}

UserDefinedBC::UserDefinedBC()
    : BoundaryCondition(0, 0, USER_DEFINED),
      filename("")
{
    is_wall_flag = false;
    sets_conv_flux_flag = false;
    sets_visc_flux_flag = false;
    // Cannot do much useful here because we don't have a filename.
}

UserDefinedBC & UserDefinedBC::operator=(const UserDefinedBC &bc)
{
    // Don't start a new interpreter for assignment to self.
    if ( this != &bc ) {
	BoundaryCondition::operator=(bc);
	is_wall_flag = bc.is_wall_flag;
	sets_conv_flux_flag = bc.sets_conv_flux_flag;
	if ( sets_conv_flux_flag ) ghost_cell_data_available = false;
	sets_visc_flux_flag = bc.sets_visc_flux_flag;
	filename = bc.filename;
	start_interpreter();
    }
    return *this;
}

UserDefinedBC::~UserDefinedBC() 
{
    lua_close(L);
}

void UserDefinedBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "filename= " << filename << endl;
    cout << lead_in << "sets_conv_flux= " << sets_conv_flux_flag << endl;
    cout << lead_in << "sets_visc_flux= " << sets_visc_flux_flag << endl;
    return;
}

int UserDefinedBC::apply_convective(double t)
{
    size_t i, j, k;
    FV_Cell *cell, *dest_cell;
    FV_Interface *IFace;
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[NORTH];
		if ( sets_conv_flux_flag ) {
		    eval_conv_flux_udf(t, i, j, k, IFace);
		} else {
		    // Ghost cells
		    eval_ghost_cell_udf(t, i, j, k, IFace);
		    dest_cell = bd.get_cell(i,j+1,k);
		    if ( fdata1 ) {
			dest_cell->copy_values_from(*fdata1);
			// Clean up memory allocated fdata1 and fdata2. 
			// This allocation was done inside unpack_flow_table()
			// which was called by eval_flux_udf().
			delete fdata1; fdata1 = 0;
		    }
		    dest_cell = bd.get_cell(i,j+2,k);
		    if ( fdata2 ) {
			dest_cell->copy_values_from(*fdata2);
			delete fdata2; fdata2 = 0;
		    }
		}
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[EAST];
		if ( sets_conv_flux_flag ) {
		    eval_conv_flux_udf(t, i, j, k, IFace);
		} else {
		    eval_ghost_cell_udf(t, i, j, k, IFace);
		    dest_cell = bd.get_cell(i+1,j,k);
		    if ( fdata1 ) {
			dest_cell->copy_values_from(*fdata1);
			delete fdata1; fdata1 = 0;
		    }
		    dest_cell = bd.get_cell(i+2,j,k);
		    if ( fdata2 ) {
			dest_cell->copy_values_from(*fdata2);
			delete fdata2; fdata2 = 0;
		    }
		}
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[SOUTH];
		if ( sets_conv_flux_flag ) {
		    eval_conv_flux_udf( t, i, j, k, IFace );
		} else {
		    eval_ghost_cell_udf(t, i, j, k, IFace);
		    dest_cell = bd.get_cell(i,j-1,k);
		    if ( fdata1 ) {
			dest_cell->copy_values_from(*fdata1);
			delete fdata1; fdata1 = 0;
		    }
		    dest_cell = bd.get_cell(i,j-2,k);
		    if ( fdata2 ) {
			dest_cell->copy_values_from(*fdata2);
			delete fdata2; fdata2 = 0;
		    }
		}
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[WEST];
		if ( sets_conv_flux_flag ) {
		    eval_conv_flux_udf(t, i, j, k, IFace);
		} else {
		    eval_ghost_cell_udf(t, i, j, k, IFace);
		    dest_cell = bd.get_cell(i-1,j,k);
		    if ( fdata1 ) {
			dest_cell->copy_values_from(*fdata1);
			delete fdata1; fdata1 = 0;
		    }
		    dest_cell = bd.get_cell(i-2,j,k);
		    if ( fdata2 ) {
			dest_cell->copy_values_from(*fdata2);
			delete fdata2; fdata2 = 0;
		    }
		}
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[TOP];
		if ( sets_conv_flux_flag ) {
		    eval_conv_flux_udf(t, i, j, k, IFace);
		} else {
		    eval_ghost_cell_udf(t, i, j, k, IFace);
		    dest_cell = bd.get_cell(i,j,k+1);
		    if ( fdata1 ) {
			dest_cell->copy_values_from(*fdata1);
			delete fdata1; fdata1 = 0;
		    }
		    dest_cell = bd.get_cell(i,j,k+2);
		    if ( fdata2 ) {
			dest_cell->copy_values_from(*fdata2);
			delete fdata2; fdata2 = 0;
		    }
		}
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[BOTTOM];
		if ( sets_conv_flux_flag ) {
		    eval_conv_flux_udf(t, i, j, k, IFace);
		} else {
		    eval_ghost_cell_udf(t, i, j, k, IFace);
		    dest_cell = bd.get_cell(i,j,k-1);
		    if ( fdata1 ) {
			dest_cell->copy_values_from(*fdata1);
			delete fdata1; fdata1 = 0;
		    }
		    dest_cell = bd.get_cell(i,j,k-2);
		    if ( fdata2 ) {
			dest_cell->copy_values_from(*fdata2);
			delete fdata2; fdata2 = 0;
		    }
		}
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_inviscid not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    } // end switch

    return SUCCESS;
} // end apply_convective()

int UserDefinedBC::apply_viscous( double t )
{
    size_t i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[NORTH];
		eval_iface_udf(t, i, j, k, IFace, cell);
		if ( sets_visc_flux_flag ) {
		    eval_visc_flux_udf(t, i, j, k, IFace);
		}
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[EAST];
		eval_iface_udf(t, i, j, k, IFace, cell);
		if ( sets_visc_flux_flag ) {
		    eval_visc_flux_udf(t, i, j, k, IFace);
		}
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[SOUTH];
		eval_iface_udf(t, i, j, k, IFace, cell);
		if ( sets_visc_flux_flag ) {
		    eval_visc_flux_udf(t, i, j, k, IFace);
		}
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[WEST];
		eval_iface_udf(t, i, j, k, IFace, cell);
		if ( sets_visc_flux_flag ) {
		    eval_visc_flux_udf(t, i, j, k, IFace);
		}
	    } // end j loop
	} // for k
 	break;
    case BOTTOM:
	k = bd.kmin;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[BOTTOM];
		eval_iface_udf(t, i, j, k, IFace, cell);
		if ( sets_visc_flux_flag ) {
		    eval_visc_flux_udf(t, i, j, k, IFace);
		}
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[TOP];
		eval_iface_udf(t, i, j, k, IFace, cell);
		if ( sets_visc_flux_flag ) {
		    eval_visc_flux_udf(t, i, j, k, IFace);
		}
	    } // end j loop
	} // for k
 	break;
    default:
	printf( "Error: apply_viscous not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    }
    return SUCCESS;
} // end UserDefinedBC::apply_viscous()

int UserDefinedBC::start_interpreter()
{
    // Set up an instance of the Lua interpreter and run the specified script file
    // to define the flow() and wall() functions.
    gmodel = get_gas_model_ptr();
    nsp = gmodel->get_number_of_species();
    nmodes = gmodel->get_number_of_modes();

#   if ECHO_ALL
    cout << "UserDefinedBC(): start a Lua interpreter on file " << filename << endl;
#   endif
    L = luaL_newstate();
    luaL_openlibs(L); // load the standard libraries

    // Set up some global variables that might be handy in 
    // the Lua environment.
    lua_pushinteger(L, bdp->id);  // we may need the id for the current block
    lua_setglobal(L, "block_id");
    lua_pushinteger(L, nsp);
    lua_setglobal(L, "nsp");
    lua_pushinteger(L, nmodes);
    lua_setglobal(L, "nmodes");
    lua_pushinteger(L, bdp->nni);
    lua_setglobal(L, "nni");
    lua_pushinteger(L, bdp->nnj);
    lua_setglobal(L, "nnj");
    lua_pushinteger(L, bdp->nnk);
    lua_setglobal(L, "nnk");
    lua_pushinteger(L, NORTH);
    lua_setglobal(L, "NORTH");
    lua_pushinteger(L, EAST);
    lua_setglobal(L, "EAST");
    lua_pushinteger(L, SOUTH);
    lua_setglobal(L, "SOUTH");
    lua_pushinteger(L, WEST);
    lua_setglobal(L, "WEST");
    lua_pushinteger(L, TOP);
    lua_setglobal(L, "TOP");
    lua_pushinteger(L, BOTTOM);
    lua_setglobal(L, "BOTTOM");

    // Register functions so that they are accessible 
    // from the Lua environment.
    register_luafns(L);

    // Presume that the user-defined functions are in the specified file.
    if ( luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0) ) {
	handle_lua_error(L, "Could not run user file: %s", lua_tostring(L, -1));
    }
    lua_settop(L, 0); // clear the stack
    return SUCCESS;
} // end start_interpreter()

int UserDefinedBC::eval_conv_flux_udf(double t, size_t i, size_t j, size_t k, FV_Interface *IFace)
{
    // Call the user-defined function which returns a table 
    // of fluxes of conserved quantities.
    // This give the user complete control over what happens at the boundary.
    // cout << "UserDefinedBC::eval_conv_flux_udf() Begin" << endl;
    double x = IFace->pos.x; 
    double y = IFace->pos.y;
    double z = IFace->pos.z;
    double csX = IFace->n.x; 
    double csY = IFace->n.y;
    double csZ = IFace->n.z;

    lua_getglobal(L, "convective_flux");  // function to be called
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, t); lua_setfield(L, -2, "t");
    lua_pushnumber(L, x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, z); lua_setfield(L, -2, "z");
    lua_pushnumber(L, csX); lua_setfield(L, -2, "csX");
    lua_pushnumber(L, csY); lua_setfield(L, -2, "csY");
    lua_pushnumber(L, csZ); lua_setfield(L, -2, "csZ");
    lua_pushinteger(L, i); lua_setfield(L, -2, "i");
    lua_pushinteger(L, j); lua_setfield(L, -2, "j");
    lua_pushinteger(L, k); lua_setfield(L, -2, "k");
    lua_pushinteger(L, which_boundary); lua_setfield(L, -2, "which_boundary");
    int number_args = 1; // table of {t x y z csX csY csZ i j k which_boundary}
    int number_results = 1; // one table of results returned on the stack.
    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L, "error running user flow function: %s\n",
			 lua_tostring(L, -1));
    }
    ConservedQuantities &F = *(IFace->F);
    // Assume that there is a single table at the TOS
    lua_getfield(L, -1, "species"); // put mass-fraction table at TOS 
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_gettable(L, -2);
	F.massf[isp] = lua_tonumber(L, -1);
	lua_pop(L, 1); // remove the number to leave the table at TOS
    }
    lua_pop(L, 1); // remove F_species table from top-of-stack
    lua_getfield(L, -1, "renergies"); // put individual-energies table at TOS 
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	lua_pushinteger(L, imode);
	lua_gettable(L, -2);
	F.energies[imode] = lua_tonumber(L, -1);
	lua_pop(L, 1); // remove the number to leave the table at TOS
    }
    lua_pop(L, 1); // remove F_species table from top-of-stack
    // cout << "Now get mass, momentum and energy fluxes." << endl;
    lua_getfield(L, -1, "mass"); F.mass = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_x"); F.momentum.x = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_y"); F.momentum.y = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_z"); F.momentum.z = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "total_energy"); F.total_energy = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "romega"); F.omega = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "rtke"); F.tke = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_settop(L, 0); // clear the stack
    return SUCCESS;
} // end eval_conv_flux_udf()

int UserDefinedBC::eval_ghost_cell_udf(double t, size_t i, size_t j, size_t k, FV_Interface *IFace)
{
    // Call the user-defined function which returns two tables 
    // of flow conditions that can be copied into the ghost cells.

    double x = IFace->pos.x; 
    double y = IFace->pos.y;
    double z = IFace->pos.z;
    double csX = IFace->n.x; 
    double csY = IFace->n.y;
    double csZ = IFace->n.z;
    // cout << "UserDefinedBC::eval_inviscid_udf() Begin" << endl;

    lua_getglobal(L, "ghost_cell");  // function to be called
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, t); lua_setfield(L, -2, "t");
    lua_pushnumber(L, x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, z); lua_setfield(L, -2, "z");
    lua_pushnumber(L, csX); lua_setfield(L, -2, "csX");
    lua_pushnumber(L, csY); lua_setfield(L, -2, "csY");
    lua_pushnumber(L, csZ); lua_setfield(L, -2, "csZ");
    lua_pushinteger(L, i); lua_setfield(L, -2, "i");
    lua_pushinteger(L, j); lua_setfield(L, -2, "j");
    lua_pushinteger(L, k); lua_setfield(L, -2, "k");
    lua_pushinteger(L, which_boundary); lua_setfield(L, -2, "which_boundary");
    int number_args = 1; // table of {t x y z csX csY csZ i j k which_boundary}
    int number_results = 2;
    // We are expecting two tables of results returned on the stack,
    // one for each ghost cell.
    // Instead of either table, the user-defined function may have
    // returned nil to indicate that the original ghost-cell data
    // should be left alone.  It may have come from an ealrier exchange.
    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L, "error running user flow function: %s\n",
			 lua_tostring(L, -1));
    }
    // get function results off the stack (in reverse-order)
    if ( lua_isnil(L, -1) ) {
	// The user-defined function has returned nil so
	// we don't want to interfere with the data that
	// may have bee set up elsewhere.
	fdata2 = 0;
    } else {
	fdata2 = unpack_flow_table(); // TOS has the second flow table
    }
    lua_pop(L, 1); // remove the second flow table to leave the first at TOS
    if ( lua_isnil(L, -1) ) {
	// The user-defined function has returned nil so
	// we don't want to interfere with the data that
	// may have bee set up elsewhere.
	fdata1 = 0;
    } else {
	fdata1 = unpack_flow_table();
    }
    lua_settop(L, 0); // clear the stack
    // cout << "UserDefinedBC::eval_inviscid_udf() End" << endl;
    return SUCCESS;
} // end eval_ghost_cell_udf()

CFlowCondition * UserDefinedBC::unpack_flow_table(void)
{
    // Unpack the table at TOS
    // This table (with named fields) represents the flow condition 
    // for a ghost cell {p u v w T massf tke omega mu_t k_t S}
    // where massf is a table of mass-fractions indexed from 0
    // cout << "UserDefinedBC::unpack_flow_table() Begin" << endl;
    std::vector<double> T;
    // cout << "Table of individual-energy-mode temperatures" << endl;
    lua_getfield(L, -1, "T"); 
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	lua_pushinteger(L, imode);
	lua_gettable(L, -2);
	T.push_back( lua_tonumber(L, -1) );
	lua_pop(L, 1); // remove the number to leave the table at TOS
    }
    lua_pop(L, 1); // remove T table from top-of-stack
    std::vector<double> massf;
    lua_getfield(L, -1, "massf"); // put mass-fraction table at TOS 
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_gettable(L, -2);
	massf.push_back( lua_tonumber(L, -1) );
	lua_pop(L, 1); // remove the number to leave the table at TOS
    }
    lua_pop(L, 1); // remove mf table from top-of-stack
    lua_getfield(L, -1, "p"); double p = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "u"); double u = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "v"); double v = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "w"); double w = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "tke"); double tke = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "omega"); double omega = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "mu_t"); double mu_t = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "k_t"); double k_t = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "S"); int S = lua_tointeger(L, -1); lua_pop(L, 1);
    if ( omega == 0.0 ) omega = 1.0;
    CFlowCondition *cfc = new CFlowCondition( gmodel, p, u, v, w, T, massf, "", 
					      tke, omega, mu_t, k_t, S);
    return cfc;
} // end unpack_flow_table()

int UserDefinedBC::eval_iface_udf(double t, size_t i, size_t j, size_t k, 
				  FV_Interface *IFace, const FV_Cell *cell)
{
    double x = IFace->pos.x; 
    double y = IFace->pos.y;
    double z = IFace->pos.z;
    double csX = IFace->n.x; 
    double csY = IFace->n.y;
    double csZ = IFace->n.z;
    FlowState &fs = *(IFace->fs);
    global_data *gdp = get_global_data_ptr();
    double dt_global = gdp->dt_global;
    size_t t_level = gdp->t_level;
    double area = IFace->area[t_level];
    
    // Call the user-defined function which leaves a table of wall conditions
    // at the top of the stack.
    lua_getglobal(L, "interface");  // function to be called
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, t); lua_setfield(L, -2, "t");
    lua_pushnumber(L, dt_global); lua_setfield(L, -2, "dt");
    lua_pushinteger(L, t_level); lua_setfield(L, -2, "t_level");
    lua_pushnumber(L, x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, z); lua_setfield(L, -2, "z");
    lua_pushnumber(L, area); lua_setfield(L, -2, "area");
    lua_pushnumber(L, csX); lua_setfield(L, -2, "csX");
    lua_pushnumber(L, csY); lua_setfield(L, -2, "csY");
    lua_pushnumber(L, csZ); lua_setfield(L, -2, "csZ");
    lua_pushinteger(L, i); lua_setfield(L, -2, "i");
    lua_pushinteger(L, j); lua_setfield(L, -2, "j");
    lua_pushinteger(L, k); lua_setfield(L, -2, "k");
    lua_pushinteger(L, which_boundary); lua_setfield(L, -2, "which_boundary");
    // Create FlowState table
    Gas_model *gmodel = get_gas_model_ptr();
    create_table_for_fs(L, fs, *gmodel); lua_setfield(L, -2, "fs");
    

    int number_args = 1; // table of
                         // {t st t_level x y z csX csY csZ i j k
                         //  which_boundary fs}
    int number_results = 1; // table of at least {u v w T_wall massf} or nil
    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L, "error running user flow function: %s\n",
			 lua_tostring(L, -1));
    }
    if ( lua_isnil(L, -1) ) {
	// The user-defined function has returned nil so
	// we don't want to interfere with the data that
	// may have bee set up elsewhere.
    } else {
	// Assume that there is a single table at the TOS
	// cout << "Table of mass fractions" << endl;
	lua_getfield(L, -1, "massf"); // put mass-fraction table at TOS 
	for ( size_t isp = 0; isp < nsp; ++isp ) {
	    lua_pushinteger(L, isp);
	    lua_gettable(L, -2);
	    fs.gas->massf[isp] = lua_tonumber(L, -1);
	    lua_pop(L, 1); // remove the number to leave the table at TOS
	}
	lua_pop(L, 1); // remove mf table from top-of-stack
	lua_getfield(L, -1, "T"); // put temperature table at TOS
	for ( size_t imode = 0; imode < nmodes; ++imode ) {
	    lua_pushinteger(L, imode);
	    lua_gettable(L, -2);
	    fs.gas->T[imode] = lua_tonumber(L, -1);
	    lua_pop(L, 1); // remove the number to leave the table at TOS
	}
	lua_pop(L, 1); // remove T table from top-of-stack
	// Now get scalar quantities
	lua_getfield(L, -1, "u"); fs.vel.x = lua_tonumber(L, -1); lua_pop(L, 1);
	lua_getfield(L, -1, "v"); fs.vel.y = lua_tonumber(L, -1); lua_pop(L, 1);
	lua_getfield(L, -1, "w"); fs.vel.z = lua_tonumber(L, -1); lua_pop(L, 1);
	lua_getfield(L, -1, "tke"); fs.tke = lua_tonumber(L, -1); lua_pop(L, 1);
	lua_getfield(L, -1, "omega"); fs.omega = lua_tonumber(L, -1); lua_pop(L, 1);
    }

    lua_settop(L, 0); // clear the stack
    return SUCCESS;
} // end eval_iface_udf()


int UserDefinedBC::eval_visc_flux_udf(double t, size_t i, size_t j, size_t k, FV_Interface *IFace)
{
    // RJG -- 26-Jul-2013
    // NOTE: The user provides viscous flux values directly BUT these
    // need to be added internally to the already present convective fluxes.
    // (If we are doing a separate convective/viscous update, then the
    // main routine should have set flux values to zero and this implementation
    // will still be correct.) Hence, the '+=' to add the viscous fluxes.

    // Call the user-defined function which returns a table 
    // of fluxes of conserved quantities.
    // This give the user complete control over what happens at the boundary.
    // cout << "UserDefinedBC::eval_visc_flux_udf() Begin" << endl;
    double x = IFace->pos.x; 
    double y = IFace->pos.y;
    double z = IFace->pos.z;
    double csX = IFace->n.x; 
    double csY = IFace->n.y;
    double csZ = IFace->n.z;

    lua_getglobal(L, "viscous_flux");  // function to be called
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, t); lua_setfield(L, -2, "t");
    lua_pushnumber(L, x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, z); lua_setfield(L, -2, "z");
    lua_pushnumber(L, csX); lua_setfield(L, -2, "csX");
    lua_pushnumber(L, csY); lua_setfield(L, -2, "csY");
    lua_pushnumber(L, csZ); lua_setfield(L, -2, "csZ");
    lua_pushinteger(L, i); lua_setfield(L, -2, "i");
    lua_pushinteger(L, j); lua_setfield(L, -2, "j");
    lua_pushinteger(L, k); lua_setfield(L, -2, "k");
    lua_pushinteger(L, which_boundary); lua_setfield(L, -2, "which_boundary");
    int number_args = 1; // table of {t x y z csX csY csZ i j k which_boundary}
    int number_results = 1; // one table of results returned on the stack.
    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L, "error running user flow function: %s\n",
			 lua_tostring(L, -1));
    }
    ConservedQuantities &F = *(IFace->F);
    // Assume that there is a single table at the TOS
    // cout << "Table of fluxes of mass fractions" << endl;
    lua_getfield(L, -1, "species"); // put mass-fraction table at TOS 
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_gettable(L, -2);
	F.massf[isp] += lua_tonumber(L, -1);
	lua_pop(L, 1); // remove the number to leave the table at TOS
    }
    lua_pop(L, 1); // remove F_species table from top-of-stack
    lua_getfield(L, -1, "renergies"); // put individual-energies table at TOS 
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	lua_pushinteger(L, imode);
	lua_gettable(L, -2);
	F.energies[imode] += lua_tonumber(L, -1);
	lua_pop(L, 1); // remove the number to leave the table at TOS
    }
    lua_pop(L, 1); // remove F_species table from top-of-stack
    lua_getfield(L, -1, "mass"); F.mass += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_x"); F.momentum.x += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_y"); F.momentum.y += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_z"); F.momentum.z += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "total_energy"); F.total_energy += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "romega"); F.omega += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "rtke"); F.tke += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_settop(L, 0); // clear the stack
    return SUCCESS;
} // end eval_visc_flux_udf()


void UserDefinedBC::handle_lua_error(lua_State *L, const char *fmt, ...)
{
    va_list argp;
    va_start(argp, fmt);
    fprintf(stderr, "Error in UserDefinedBC or AdjacentPlusUDFBC: block %d, face %d\n",
	    static_cast<int>(bdp->id), which_boundary);
    vfprintf(stderr, fmt, argp);
    fprintf(stderr, "\n");
    va_end(argp);
    lua_close(L);
    exit(LUA_ERROR);
}

//------------------------------------------------------------------------

AdjacentPlusUDFBC::AdjacentPlusUDFBC(Block *bdp, int which_boundary, 
				     int other_block, int other_face,
				     int _neighbour_orientation,
				     const std::string filename,
				     bool is_wall, 
				     bool sets_conv_flux, bool sets_visc_flux)
    : UserDefinedBC(bdp, which_boundary, filename, is_wall, sets_conv_flux, sets_visc_flux)
{
    type_code = ADJACENT_PLUS_UDF; // Needs to overwrite UserDefinedBC entry.
    neighbour_block = other_block; 
    neighbour_face = other_face;
    neighbour_orientation = _neighbour_orientation;
}

AdjacentPlusUDFBC::AdjacentPlusUDFBC(const AdjacentPlusUDFBC &bc)
    : UserDefinedBC(bc.bdp, bc.which_boundary, bc.filename, bc.is_wall_flag)
{
    type_code = bc.type_code;
    neighbour_block = bc.neighbour_block; 
    neighbour_face = bc.neighbour_face;
    neighbour_orientation = bc.neighbour_orientation;
    cerr << "AdjacentPlusUDFBC() copy constructor is not implemented FIX-ME." << endl;
    exit( NOT_IMPLEMENTED_ERROR );
}

AdjacentPlusUDFBC::~AdjacentPlusUDFBC() {}

void AdjacentPlusUDFBC::print_info(std::string lead_in)
{
    UserDefinedBC::print_info(lead_in);
    cout << lead_in << "neighbour_block= " << neighbour_block << endl;
    cout << lead_in << "neighbour_face= " << neighbour_face 
	 << " (" << get_face_name(neighbour_face) << ")" << endl;
    cout << lead_in << "neighbour_orientation= " << neighbour_orientation << endl;
    return;
}
