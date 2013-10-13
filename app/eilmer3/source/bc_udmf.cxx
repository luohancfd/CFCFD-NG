
#include "kernel.hh"
#include "bc_udmf.hh"
#include "bc_menter_correction.hh"
#include "lua_service_for_e3.hh"

UserDefinedMassFluxBC::
UserDefinedMassFluxBC(Block *bdp, int which_boundary, const std::string fname)
    : BoundaryCondition(bdp, which_boundary, USER_DEFINED_MASS_FLUX),
      filename_(fname)
{
    is_wall_flag = true;
    sets_conv_flux_flag = false;
    sets_visc_flux_flag = false; // in the general sense
                                 // we specially intercept the
                                 // the diffusive species fluxes 
                                 // in visc.cxx
    start_interpreter();
}

UserDefinedMassFluxBC::
UserDefinedMassFluxBC(const UserDefinedMassFluxBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      filename_(bc.filename_)
{
    is_wall_flag = bc.is_wall_flag;
    sets_conv_flux_flag = bc.sets_conv_flux_flag;
    sets_visc_flux_flag = bc.sets_visc_flux_flag;
    start_interpreter();
}

UserDefinedMassFluxBC::
UserDefinedMassFluxBC()
    : BoundaryCondition(0, 0, USER_DEFINED),
      filename_("")
{
    is_wall_flag = false;
    sets_conv_flux_flag = false;
    sets_visc_flux_flag = false;
    // Cannot do much useful here because we don't have a filename.
}

UserDefinedMassFluxBC&
UserDefinedMassFluxBC::
operator=(const UserDefinedMassFluxBC &bc)
{
    if ( this != &bc ) {
	BoundaryCondition::operator=(bc);
	is_wall_flag = bc.is_wall_flag;
	sets_conv_flux_flag = bc.sets_conv_flux_flag;
	sets_visc_flux_flag = bc.sets_visc_flux_flag;
	filename_ = bc.filename_;
	start_interpreter();
    }
    return *this;
}

UserDefinedMassFluxBC::
~UserDefinedMassFluxBC()
{
    lua_close(L);
}

void
UserDefinedMassFluxBC::
print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "filename= " << filename_ << endl;
    cout << lead_in << "sets_conv_flux= " << sets_conv_flux_flag << endl;
    cout << lead_in << "sets_visc_flux= " << sets_visc_flux_flag << endl;
    return;
}


int
UserDefinedMassFluxBC::
apply_viscous(double t)
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
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		// From user function, set temperature and mass fluxes
		eval_user_fn(t, i, j, k, IFace);
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[EAST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		// From user function, set temperature and mass fluxes
		eval_user_fn(t, i, j, k, IFace);
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[SOUTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		// From user function, set temperature and mass fluxes
		eval_user_fn(t, i, j, k, IFace);
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[WEST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		// From user function, set temperature and mass fluxes
		eval_user_fn(t, i, j, k, IFace);
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[TOP];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		// From user function, set temperature and mass fluxes
		eval_user_fn(t, i, j, k, IFace);
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[BOTTOM];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		// From user function, set temperature and mass fluxes
		eval_user_fn(t, i, j, k, IFace);
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_viscous not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    }
    return SUCCESS;
}

int
UserDefinedMassFluxBC::
start_interpreter()
{
    // COPIED FROM bc_user_defined.cxx -- need to refactor out
    // Set up an instance of the Lua interpreter and run the specified script file
    // to define the flow() and wall() functions.
    gmodel = get_gas_model_ptr();
    nsp = gmodel->get_number_of_species();
    nmodes = gmodel->get_number_of_modes();

#   if ECHO_ALL
    cout << "UserDefinedMassFluxBC(): start a Lua interpreter on file " << filename_ << endl;
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
    if ( luaL_loadfile(L, filename_.c_str()) || lua_pcall(L, 0, 0, 0) ) {
	handle_lua_error(L, "Could not run user file: %s", lua_tostring(L, -1));
    }
    lua_settop(L, 0); // clear the stack
    return SUCCESS;
}

int
UserDefinedMassFluxBC::
eval_user_fn(double t, size_t i, size_t j, size_t k, FV_Interface *IFace)
{
    // We are going to hand the user a table of (possibly) useful
    // values related to the interface we are working on.
    double x = IFace->pos.x; 
    double y = IFace->pos.y;
    double z = IFace->pos.z;
    double csX = IFace->n.x; 
    double csY = IFace->n.y;
    double csZ = IFace->n.z;

    lua_getglobal(L, "mass_flux");
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

    // We expect from the user a table with:
    // t = { Twall = ....,
    //       species = .....
    //      }

    lua_getfield(L, -1, "Twall");
    double Twall = lua_tonumber(L, -1);
    lua_pop(L, 1);
    // Now set all temperatures to Twall value.
    for ( size_t imode = 0; imode < nmodes; ++imode )
	IFace->fs->gas->T[imode] = Twall;

    ConservedQuantities &F = *(IFace->F);
    lua_getfield(L, -1, "species");
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_gettable(L, -2);
	double massf_i = lua_tonumber(L, -1);
	F.massf[isp] += massf_i;
	// And add to total mass (if ablating, for example)
	F.mass += massf_i;
	lua_pop(L, 1); // remove the number to leave the table at TOS
    }
    lua_pop(L, 1); // remove 'species' table from top-of-stack

    lua_settop(L, 0); // clear the stack
    return SUCCESS;
}

int
UserDefinedMassFluxBC::
handle_lua_error(lua_State *L, const char *fmt, ...)
{
    va_list argp;
    va_start(argp, fmt);
    fprintf(stderr, "Error in UserDefinedMassFluxBC: block %d, face %d\n",
	    static_cast<int>(bdp->id), which_boundary);
    vfprintf(stderr, fmt, argp);
    fprintf(stderr, "\n");
    va_end(argp);
    lua_close(L);
    exit(LUA_ERROR);
}



