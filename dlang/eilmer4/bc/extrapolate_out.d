// bc/extrapolate_out.d
//
// Outflow boundary condition that simply extrapolates the flow near the boundary
// into the ghost cells just outside the boundary.
//
// Peter J. 2014-07-26

import std.conv;

import gas;
import fvcore;
import fvinterface;
import fvcell;
import bc;
import block;
import sblock;
import globalconfig;
import globaldata;

class ExtrapolateOutBC: BoundaryCondition {
public:
    int x_order = 0; // default to lowest order

    this(int id, int boundary, int x_order_=0) 
    {
	blk_id = id;
	which_boundary = boundary;
	type_code = BCCode.extrapolate_out;
	is_wall = false;
	ghost_cell_data_available = true;
	sets_conv_flux_directly = false;
	sets_visc_flux_directly = false;
	emissivity = 0.0;
	x_order = x_order_;
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "ExtrapolateOutBC(";
	repr ~= "x_order=" ~ to!string(x_order);
	repr ~= ")";
	return to!string(repr);
    }

    override void apply_convective(double t)
    {
	// Fill ghost cells with data from just inside the boundary
	// using zero-order extrapolation (i.e. just copy the data).
	// We assume that this boundary is an outflow boundary.
	size_t i, j, k;
	FVCell src_cell, dest_cell;
	FVCell cell_1, cell_2;
	auto gmodel = GlobalConfig.gmodel;
	size_t nsp = gmodel.n_species;
	size_t nmodes = gmodel.n_modes;
	auto blk = allBlocks[blk_id];

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    if ( x_order == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (j-1)        (j)           (j+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j-1,k);
			dest_cell = blk.get_cell(i,j+1,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (j)        (j+1)       (j+2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j+2,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    }
		    else {
			// Zero-order extrapolation
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j+1,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j+2,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    } 
		} // end i loop
	    } // for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( x_order == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (i-1)        (i)           (i+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i-1,j,k);
			dest_cell = blk.get_cell(i+1,j,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (i)        (i+1)       (i+2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i+2,j,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    }
		    else {
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i+1,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i+2,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    if ( x_order == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (j+1)        (j)           (j-1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j+1,k);
			dest_cell = blk.get_cell(i,j-1,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (j)        (j-1)       (j-2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j-2,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j-1,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j-2,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end i loop
	    } // for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( x_order == 1 ) {
			//  ---[ghost cell 2]---|--- [dest] ---|||--- [1] ---|---[2]----
			//      (i-2)                 (i-1)           (i)       (i+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i+1,j,k);
			dest_cell = blk.get_cell(i-1,j,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[dest]---|---[1]---|||---[2]---|------|
			//       (i-2)       (i-1)       (i)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i-2,j,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			// Zero-order extrapolation
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i-1,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i-2,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( x_order == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (k-1)        (k)           (k+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j,k-1);
			dest_cell = blk.get_cell(i,j,k+1);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (k)        (k+1)       (k+2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j,k+2);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			// Zero-order extrapolation
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j,k+1);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j,k+2);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( x_order == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (k+1)        (k)           (k-1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j,k+2);
			dest_cell = blk.get_cell(i,j,k-1);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (k)        (k-1)       (k-2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j,k-2);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j,k-1);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j,k-2);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for i
	    break;
	} // end switch
    } // end apply_convective()

    // Let the base class implementations do the work.
    // apply_viscous

} // end class ExtrapolateOutBC
