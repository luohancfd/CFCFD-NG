#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"

#include "bc_fixed_q.hh"
#include "solid_cell.hh"

BC_FIXED_Q::BC_FIXED_Q(double qwall)
    : qwall(qwall)
{}

void BC_FIXED_Q::apply_bc(SolidBlock *blk, int which_boundary )
{
    std::vector<double> length(2,0.0);
    switch (which_boundary){
    case 0: //South boundary
	{ //Shouldn't have to do anything for the anisotropic cases as incident flux
	    int j = 0;
	    for (int i = 0; i<blk->nni-1 ; i++){
		
		//Set flux direction and scale.
		blk->block_cells[i][j].iface[0]->q = blk->block_cells[i][j].iface[0]->normal;
		scale_vector(blk->block_cells[i][j].iface[0]->q,qwall);
		
		//Solve for interface temperature using conduction eqn.
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[0]->pos);
		blk->block_cells[i][j].iface[0]->T = -1*(-qwall/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) - blk->block_cells[i][j].T );
	    }
	    break;
	}
    case 1:
	{
	    int i = blk->nni-2;
	    for (int j = 0; j<blk->nnj-1 ; j++){
		
		blk->block_cells[i][j].iface[1]->q = blk->block_cells[i][j].iface[1]->normal;
		scale_vector(blk->block_cells[i][j].iface[1]->q,-1 * qwall); //flux inwards
		
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[1]->pos);
		blk->block_cells[i][j].iface[1]->T = (-qwall/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) + blk->block_cells[i][j].T );
	    }
	    break;
	}
    case 2:
	{
	    int j = blk->nnj-2;
	    for (int i = 0; i<blk->nni-1 ; i++){
		
		blk->block_cells[i][j].iface[2]->q = blk->block_cells[i][j].iface[2]->normal;
		scale_vector(blk->block_cells[i][j].iface[2]->q,-1 * qwall); //flux inwards
		
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[2]->pos);
		blk->block_cells[i][j].iface[2]->T = (-qwall/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) + blk->block_cells[i][j].T );
	    }
	    break;
	}
    case 3:
	{
	    int i = 0;
	    for (int j = 0; j<blk->nnj-1 ; j++){
		
		blk->block_cells[i][j].iface[3]->q = blk->block_cells[i][j].iface[3]->normal;
		scale_vector(blk->block_cells[i][j].iface[3]->q,qwall);
		
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[3]->pos);
		blk->block_cells[i][j].iface[3]->T = -1*(-qwall/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) - blk->block_cells[i][j].T );
	    }
        break;
	}
    }
}

void BC_FIXED_Q::print_type(){
    std::cout << "Fixed flux BC: q = " << qwall << "W/m2" <<  std::endl;
}
