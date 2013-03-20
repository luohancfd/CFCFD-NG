// bc_extrapolate_out.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_extrapolate_out.hh"
#include "kernel.hh"

//------------------------------------------------------------------------

ExtrapolateOutBC::ExtrapolateOutBC( Block *bdp, int which_boundary, int x_order, int _sponge_flag )
    : BoundaryCondition(bdp, which_boundary, EXTRAPOLATE_OUT, "ExtrapolateOutBC",
			x_order, false, false, -1, -1, 0) 
{}

ExtrapolateOutBC::ExtrapolateOutBC( const ExtrapolateOutBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation) 
{}

ExtrapolateOutBC::ExtrapolateOutBC()
    : BoundaryCondition(0, 0, EXTRAPOLATE_OUT, "ExtrapolateOutBC",
			0, false, false, -1, -1, 0) 
{}

ExtrapolateOutBC & ExtrapolateOutBC::operator=(const ExtrapolateOutBC &bc)
{
    BoundaryCondition::operator=(bc);
    return *this;
}

ExtrapolateOutBC::~ExtrapolateOutBC() {}

int ExtrapolateOutBC::apply_inviscid( double t )
{
    // Fill ghost cells with data from just inside the boundary
    // using zero-order extrapolation (i.e. just copy the data).
    // We assume that this boundary is an outflow boundary.
    size_t i, j, k;
    FV_Cell *src_cell, *dest_cell;
    FV_Cell *cell_1, *cell_2;
    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t nmodes = gmodel->get_number_of_modes();
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		if ( x_order == 1 ) {
		    //  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
		    //      (j-1)        (j)           (j+1)
		    //  dest: ghost cell 1
		    //  [1]: first interior cell
		    //  [2]: second interior cell
		    // This extrapolation assumes that cell-spacing between
		    // cells 1 and 2 continues on in the exterior
		    cell_1 = bd.get_cell(i,j,k);
		    cell_2 = bd.get_cell(i,j-1,k);
		    dest_cell = bd.get_cell(i,j+1,k);
		    // Extrapolate on primitive variables
		    // 1. First exterior point
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		    // 2. Second exterior point
		    //  |---[2]---|||---[1]---|---[dest]------
		    //      (j)        (j+1)       (j+2)
		    cell_2 = cell_1;
		    cell_1 = dest_cell;
		    dest_cell = bd.get_cell(i,j+2,k);
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		}
		else {
		    // Zero-order extrapolation
		    src_cell = bd.get_cell(i,j,k);
		    dest_cell = bd.get_cell(i,j+1,k);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		    dest_cell = bd.get_cell(i,j+2,k);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		} 
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		if ( x_order == 1 ) {
		    //  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
		    //      (i-1)        (i)           (i+1)
		    //  dest: ghost cell 1
		    //  [1]: first interior cell
		    //  [2]: second interior cell
		    // This extrapolation assumes that cell-spacing between
		    // cells 1 and 2 continues on in the exterior
		    cell_1 = bd.get_cell(i,j,k);
		    cell_2 = bd.get_cell(i-1,j,k);
		    dest_cell = bd.get_cell(i+1,j,k);
		    // Extrapolate on primitive variables
		    // 1. First exterior point
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		    // 2. Second exterior point
		    //  |---[2]---|||---[1]---|---[dest]------
		    //      (i)        (i+1)       (i+2)
		    cell_2 = cell_1;
		    cell_1 = dest_cell;
		    dest_cell = bd.get_cell(i+2,j,k);
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		}
		else {
		    src_cell = bd.get_cell(i,j,k);
		    dest_cell = bd.get_cell(i+1,j,k);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		    dest_cell = bd.get_cell(i+2,j,k);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		}
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		if ( x_order == 1 ) {
		    //  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
		    //      (j+1)        (j)           (j-1)
		    //  dest: ghost cell 1
		    //  [1]: first interior cell
		    //  [2]: second interior cell
		    // This extrapolation assumes that cell-spacing between
		    // cells 1 and 2 continues on in the exterior
		    cell_1 = bd.get_cell(i,j,k);
		    cell_2 = bd.get_cell(i,j+1,k);
		    dest_cell = bd.get_cell(i,j-1,k);
		    // Extrapolate on primitive variables
		    // 1. First exterior point
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		    // 2. Second exterior point
		    //  |---[2]---|||---[1]---|---[dest]------
		    //      (j)        (j-1)       (j-2)
		    cell_2 = cell_1;
		    cell_1 = dest_cell;
		    dest_cell = bd.get_cell(i,j-2,k);
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		}
		else {
		    src_cell = bd.get_cell(i,j,k);
		    dest_cell = bd.get_cell(i,j-1,k);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		    dest_cell = bd.get_cell(i,j-2,k);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		}
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		if ( x_order == 1 ) {
		    //  ---[ghost cell 2]---|--- [dest] ---|||--- [1] ---|---[2]----
		    //      (i-2)                 (i-1)           (i)       (i+1)
		    //  dest: ghost cell 1
		    //  [1]: first interior cell
		    //  [2]: second interior cell
		    // This extrapolation assumes that cell-spacing between
		    // cells 1 and 2 continues on in the exterior
		    cell_1 = bd.get_cell(i,j,k);
		    cell_2 = bd.get_cell(i+1,j,k);
		    dest_cell = bd.get_cell(i-1,j,k);
		    // Extrapolate on primitive variables
		    // 1. First exterior point
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		    // 2. Second exterior point
		    //  |---[dest]---|---[1]---|||---[2]---|------|
		    //       (i-2)       (i-1)       (i)
		    cell_2 = cell_1;
		    cell_1 = dest_cell;
		    dest_cell = bd.get_cell(i-2,j,k);
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		}
		else {
		    // Zero-order extrapolation
		    src_cell = bd.get_cell(i,j,k);
		    dest_cell = bd.get_cell(i-1,j,k);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		    dest_cell = bd.get_cell(i-2,j,k);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		}
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		if ( x_order == 1 ) {
		    //  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
		    //      (k-1)        (k)           (k+1)
		    //  dest: ghost cell 1
		    //  [1]: first interior cell
		    //  [2]: second interior cell
		    // This extrapolation assumes that cell-spacing between
		    // cells 1 and 2 continues on in the exterior
		    cell_1 = bd.get_cell(i,j,k);
		    cell_2 = bd.get_cell(i,j,k-1);
		    dest_cell = bd.get_cell(i,j,k+1);
		    // Extrapolate on primitive variables
		    // 1. First exterior point
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		    // 2. Second exterior point
		    //  |---[2]---|||---[1]---|---[dest]------
		    //      (k)        (k+1)       (k+2)
		    cell_2 = cell_1;
		    cell_1 = dest_cell;
		    dest_cell = bd.get_cell(i,j,k+2);
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		}
		else {
		    // Zero-order extrapolation
		    src_cell = bd.get_cell(i,j,k);
		    dest_cell = bd.get_cell(i,j,k+1);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		    dest_cell = bd.get_cell(i,j,k+2);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		}
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		if ( x_order == 1 ) {
		    //  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
		    //      (k+1)        (k)           (k-1)
		    //  dest: ghost cell 1
		    //  [1]: first interior cell
		    //  [2]: second interior cell
		    // This extrapolation assumes that cell-spacing between
		    // cells 1 and 2 continues on in the exterior
		    cell_1 = bd.get_cell(i,j,k);
		    cell_2 = bd.get_cell(i,j,k+2);
		    dest_cell = bd.get_cell(i,j,k-1);
		    // Extrapolate on primitive variables
		    // 1. First exterior point
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		    // 2. Second exterior point
		    //  |---[2]---|||---[1]---|---[dest]------
		    //      (k)        (k-1)       (k-2)
		    cell_2 = cell_1;
		    cell_1 = dest_cell;
		    dest_cell = bd.get_cell(i,j,k-2);
		    dest_cell->fs->gas->rho = 2.0*cell_1->fs->gas->rho - cell_2->fs->gas->rho;
		    for ( size_t imode = 0; imode < nmodes; ++imode ) {
			dest_cell->fs->gas->e[imode] = 2.0*cell_1->fs->gas->e[imode] - cell_2->fs->gas->e[imode];
		    }
		    if ( nsp > 1 ) {
			for ( size_t isp = 0; isp < nsp; ++isp ) {
			    dest_cell->fs->gas->massf[isp] = 2.0*cell_1->fs->gas->massf[isp] - cell_2->fs->gas->massf[isp];
			}
			scale_mass_fractions(dest_cell->fs->gas->massf);
		    }
		    else {
			dest_cell->fs->gas->massf[0] = 1.0;
		    }
		    gmodel->eval_thermo_state_rhoe(*(dest_cell->fs->gas));
		    dest_cell->fs->vel = 2.0*cell_1->fs->vel - cell_2->fs->vel;
		    dest_cell->fs->tke = 2.0*cell_1->fs->tke - cell_2->fs->tke;
		    dest_cell->fs->omega = 2.0*cell_1->fs->omega - cell_2->fs->omega;
		    dest_cell->fs->mu_t = 2.0*cell_1->fs->mu_t - cell_2->fs->mu_t;
		    dest_cell->fs->k_t = 2.0*cell_1->fs->k_t - cell_2->fs->k_t;
		}
		else {
		    src_cell = bd.get_cell(i,j,k);
		    dest_cell = bd.get_cell(i,j,k-1);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		    dest_cell = bd.get_cell(i,j,k-2);
		    dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
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
}
