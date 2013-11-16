#include "../../../lib/util/source/useful.hh"
#include "kernel.hh"
#include "bc_conjugate_ht.hh"

ConjugateHeatTransferBC::
ConjugateHeatTransferBC(Block *bdp, int which_boundary)
    : BoundaryCondition(bdp, which_boundary, CONJUGATE_HT)
{
    sets_conv_flux = false;
    ghost_cell_data_available = false;
    sets_visc_flux = false;
    if ( which_boundary != NORTH ) {
	cout << "Error in boundary condition specification.\n";
	cout << "A ConjugateHT boundary condition can only be applied\n";
	cout << "to a NORTH wall. We detected it being applied to:\n";
	cout << which_boundary << endl;
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
}

ConjugateHeatTransferBC::
ConjugateHeatTransferBC(const ConjugateHeatTransferBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      fname_(bc.fname_)
{
    sets_conv_flux = false;
    ghost_cell_data_available = false;
    sets_visc_flux = false;
}

ConjugateHeatTransferBC::
ConjugateHeatTransferBC()
    : BoundaryCondition(0, 0, CONJUGATE_HT),
      fname_("")
{
    sets_conv_flux = false;
    ghost_cell_data_available = false;
    sets_visc_flux = false;
    // Cannot do much useful here because we don't have a filename.
}

ConjugateHeatTransferBC&
ConjugateHeatTransferBC::
operator=(const ConjugateHeatTransferBC &bc)
{
    if ( this != &bc ) {
	BoundaryCondition::operator=(bc);
	sets_conv_flux_flag = bc.sets_conv_flux_flag;
	ghost_cell_data_available = bc.ghost_cell_data_available
	sets_visc_flux_flag = bc.sets_visc_flux_flag;
	fname_ = bc.fname_;
	T_.resize(bc.T_.size());
	q_.resize(bc.q_.size());
	initialise_wall_model(fname_);
    }
    return *this;
}

ConjugateHeatTransferBC::
~ConjugateHeatTransferBC() {}


void
ConjugateHeatTransferBC::
print_info(string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "fname= " << fname_ << endl;
    return;
}

int
ConjugateHeatTransferBC::
apply_viscous(double t)
{
    // This mostly looks like the fixed-temperature boundary condition.
    // What differs is that the temperature values at the
    // interface are supplied by the wall model.
    size_t i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();
    Block & bd = *bdp;
    global_data *G = get_global_data_ptr();
    int rank = G.my_mpi_rank;
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();
    // ONLY IMPLEMENTED FOR NORTH BOUNDARY

    j = bd.jmax;
    for (k = bd.kmin; k <= bd.kmax; ++k) {
	iT = G.displs[rank];
	for (i = bd.imin; i <= bd.imax; ++i) {
	    cell = bd.get_cell(i,j,k);
	    IFace = cell->iface[NORTH];
	    FlowState &fs = *(IFace->fs);
	    fs.copy_values_from(*(cell->fs));
	    fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
	    fs.tke = 0.0;
	    fs.omega = ideal_omega_at_wall(cell);
	    if (bd.bcp[NORTH]->wc_bc != NON_CATALYTIC) {
		cw->apply(*(cell->fs->gas), fs.gas->massf);
	    }
	    for ( size_t imode = 0; imode < nmodes; ++imode ) {
		fs.gas->T[imode] = G.T_wall[iT];
	    }
	    ++iT;
	} // end i loop
    } // for k
    return SUCCESS;
}
