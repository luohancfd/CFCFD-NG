#include <iostream>

#include "bc_adiabatic.hh"
#include "solid_cell.hh"

BC_ADIABATIC::BC_ADIABATIC()
{}

void BC_ADIABATIC::apply_bc(SolidBlock *blk, int which_boundary){
    //Set boundary cell interface fluxes to zero
    //This will not work for both isotropic and anisotropic
    switch (which_boundary){
    case 0: //South boundary
	{
	    int j = 0;
	    for (int i = 0; i<blk->nni-1 ; i++){
		blk->block_cells[i][j].iface[0]->T = blk->block_cells[i][j].T;
		blk->block_cells[i][j].iface[0]->q[0] = 0.0;
		blk->block_cells[i][j].iface[0]->q[1] = 0.0;
	    }
	    break;
	}
    case 1: //East
	{
	    int i = blk->nni-2;
	    for (int j = 0; j<blk->nnj-1 ; j++){
		blk->block_cells[i][j].iface[1]->T = blk->block_cells[i][j].T;
		blk->block_cells[i][j].iface[1]->q[0] =0.0;
		blk->block_cells[i][j].iface[1]->q[1] = 0.0;
	    }
	    break;
	}
    case 2:
	{
	    int j = blk->nnj-2;
	    for (int i = 0; i<blk->nni-1 ; i++){
		blk->block_cells[i][j].iface[2]->T = blk->block_cells[i][j].T;
		blk->block_cells[i][j].iface[2]->q[0] =0.0;
		blk->block_cells[i][j].iface[2]->q[1] = 0.0;
	    }
	    break;
	}
    case 3:
	{
	    int i = 0;
	    for (int j = 0; j<blk->nnj-1 ; j++){
		blk->block_cells[i][j].iface[3]->T = blk->block_cells[i][j].T;
		blk->block_cells[i][j].iface[3]->q[0] =0.0;
		blk->block_cells[i][j].iface[3]->q[1] = 0.0;
	    }
	    break;
	}
    }
}


void BC_ADIABATIC::print_type(){
    std::cout << "Adiabatic BC" <<  std::endl;
}
