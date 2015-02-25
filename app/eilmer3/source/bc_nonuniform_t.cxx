// bc_nonuniform_t.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "bc_nonuniform_t.hh"
#include "bc_catalytic.hh"
#include "bc_menter_correction.hh"

//------------------------------------------------------------------------

NonuniformTBC::NonuniformTBC(Block *bdp, int which_boundary, vector<double>& T_non_,
                             int starting_blk, vector<double>& no_blk_, 
                             Vector3 _r_omega, Vector3 _centre, Vector3 _v_trans,
                             double _emissivity)
    : BoundaryCondition(bdp, which_boundary, NONUNIFORM_T),
      T_non(T_non_), starting_blk(starting_blk), no_blk(no_blk_),
      r_omega(_r_omega), centre(_centre), v_trans(_v_trans)
{
    is_wall_flag = true;
    emissivity = _emissivity;
}

NonuniformTBC::NonuniformTBC(const NonuniformTBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      T_non(bc.T_non), starting_blk(bc.starting_blk), no_blk(bc.no_blk),
      r_omega(bc.r_omega), centre(bc.centre), v_trans(bc.v_trans)
{
    is_wall_flag = bc.is_wall_flag;
    emissivity = bc.emissivity;
}

NonuniformTBC::NonuniformTBC()
    : BoundaryCondition(0, 0, NONUNIFORM_T),
      r_omega(Vector3(0.0,0.0,0.0)), centre(Vector3(0.0,0.0,0.0)), v_trans(Vector3(0.0,0.0,0.0))
{
    is_wall_flag = true;
    emissivity = 1.0;
}

NonuniformTBC & NonuniformTBC::operator=(const NonuniformTBC &bc)
{
    BoundaryCondition::operator=(bc);
    T_non = bc.T_non;
    starting_blk = bc.starting_blk;
    no_blk = bc.no_blk;
    r_omega.x = bc.r_omega.x; r_omega.y = bc.r_omega.y; r_omega.z = bc.r_omega.z;
    centre.x = bc.centre.x; centre.y = bc.centre.y; centre.z = bc.centre.z;
    v_trans.x = bc.v_trans.x; v_trans.y = bc.v_trans.y; v_trans.z = bc.v_trans.z;
    return *this;
}

NonuniformTBC::~NonuniformTBC() {}

void NonuniformTBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "T_non= " << T_non << endl;
    cout << lead_in << "starting_blk= " << starting_blk << endl;
    cout << lead_in << "no_blk= " << no_blk << endl;
    cout << lead_in << "r_omega= " << r_omega << endl;
    cout << lead_in << "centre= " << centre << endl;
    cout << lead_in << "v_trans= " << v_trans << endl;
    return;
}

int NonuniformTBC::apply_viscous(double t)
{
    size_t i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();
    Block & bd = *bdp;
    global_data &G = *get_global_data_ptr();     
    
    size_t imax = bd.imax - bd.imin + 1 ;    // grid number in i direction in a single block
    size_t jmax = bd.jmax - bd.jmin + 1 ;    // grid number in j direction in a single block
    //size_t kmax = bd.kmax - bd.kmin + 1 ;  // grid number in k direction in a single block
    size_t blk_id = bdp->id - starting_blk;  // current block id
    size_t imax_all = imax * no_blk[0] ;     // grid number in i direction
    size_t jmax_all = jmax * no_blk[1] ;     // grid number in j direction
    //size_t kmax_all = kmax * no_blk[2] ;   // grid number in k direction
    size_t index = 0;                        // index number for nonuniform wall temperature

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
                index = (k-bd.kmin)*imax + (i-bd.imin) + blk_id*imax + (k-bd.kmin)*(imax_all-imax);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[NORTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                if ( G.moving_grid ) {
		        IFace->ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
		        fs.vel.transform_to_local(IFace->n, IFace->t1, IFace->t2);		    
		        fs.vel.y += IFace->ivel.y;
		        fs.vel.z += IFace->ivel.z;
                        fs.vel.transform_to_global(IFace->n, IFace->t1, IFace->t2);		    
                        IFace->ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);
                }
		for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = T_non[index];
		// [TODO] should we re-evaluate the thermo and transport coeffs?
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[NORTH]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
                index = (k-bd.kmin)*jmax + (j-bd.jmin) + blk_id*jmax + (k-bd.kmin)*(jmax_all-jmax);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[EAST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                if ( G.moving_grid ) {
		        IFace->ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
		        fs.vel.transform_to_local(IFace->n, IFace->t1, IFace->t2);		    
		        fs.vel.y += IFace->ivel.y;
		        fs.vel.z += IFace->ivel.z;
                        fs.vel.transform_to_global(IFace->n, IFace->t1, IFace->t2);		    
                        IFace->ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);
                }               
		for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = T_non[index];
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[EAST]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
                index = (k-bd.kmin)*imax + (i-bd.imin) + blk_id*imax + (k-bd.kmin)*(imax_all-imax);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[SOUTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                if ( G.moving_grid ) {
		        IFace->ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
		        fs.vel.transform_to_local(IFace->n, IFace->t1, IFace->t2);		    
		        fs.vel.y += IFace->ivel.y;
		        fs.vel.z += IFace->ivel.z;
                        fs.vel.transform_to_global(IFace->n, IFace->t1, IFace->t2);		    
                        IFace->ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);
                }            
		for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = T_non[index];
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[SOUTH]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
                index = (k-bd.kmin)*jmax + (j-bd.jmin) + blk_id*jmax + (k-bd.kmin)*(jmax_all-jmax);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[WEST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                if ( G.moving_grid ) {
		        IFace->ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
		        fs.vel.transform_to_local(IFace->n, IFace->t1, IFace->t2);		    
		        fs.vel.y += IFace->ivel.y;
		        fs.vel.z += IFace->ivel.z;
                        fs.vel.transform_to_global(IFace->n, IFace->t1, IFace->t2);		    
                        IFace->ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);
                }              
		for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = T_non[index];
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[WEST]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
                index = (i-bd.imin)*jmax + (j-bd.jmin) + blk_id*jmax + (i-bd.imin)*(jmax_all-jmax);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[TOP];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                if ( G.moving_grid ) {
		        IFace->ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
		        fs.vel.transform_to_local(IFace->n, IFace->t1, IFace->t2);		    
		        fs.vel.y += IFace->ivel.y;
		        fs.vel.z += IFace->ivel.z;
                        fs.vel.transform_to_global(IFace->n, IFace->t1, IFace->t2);		    
                        IFace->ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);
                }               
		for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = T_non[index];
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[TOP]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
                index = (i-bd.imin)*jmax + (j-bd.jmin) + blk_id*jmax + (i-bd.imin)*(jmax_all-jmax);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[BOTTOM];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                if ( G.moving_grid ) {
		        IFace->ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
		        fs.vel.transform_to_local(IFace->n, IFace->t1, IFace->t2);		    
		        fs.vel.y += IFace->ivel.y;
		        fs.vel.z += IFace->ivel.z;
                        fs.vel.transform_to_global(IFace->n, IFace->t1, IFace->t2);		    
                        IFace->ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);
                }               
		for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = T_non[index];
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[BOTTOM]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_viscous not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    }
    return SUCCESS;
} // end FixedTBC::apply_viscous()
