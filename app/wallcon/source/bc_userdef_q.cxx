#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"

#include "bc_userdef_q.hh"
#include "solid_cell.hh"
#include "useful_funcs.hh"


BC_USERDEF_Q::BC_USERDEF_Q(){}

BC_USERDEF_Q::BC_USERDEF_Q(std::vector<double> qwall)
    : qwall(qwall)
{}

void BC_USERDEF_Q::apply_bc(SolidBlock *blk, int which_boundary )
{
    std::vector<double> length(2,0.0);
    //Flux defined in vector qwall.
    switch (which_boundary){
    case 0: //South boundary
	{
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc("south", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	this->qwall = terms; //"this" is pointer to the current block
	    }

	    int j = 0;
	    for (int i = 0; i<blk->nni-1 ; i++){

		//Set flux direction and scale
		blk->block_cells[i][j].iface[0]->q = blk->block_cells[i][j].iface[0]->normal;
		scale_vector(blk->block_cells[i][j].iface[0]->q,qwall[i]);
        
		//Solve for interface temperature
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[0]->pos);
		blk->block_cells[i][j].iface[0]->T = -1*(-qwall[i]/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) - blk->block_cells[i][j].T );
	    }
	    break;
	}
    case 1:
	{
	    //Read in time varying terms if flag set
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc("east", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	this->qwall = terms; //"this" is pointer to the current block
	    }

	    int i = blk->nni-2;
	    for (int j = 0; j<blk->nnj-1 ; j++){
  
		blk->block_cells[i][j].iface[1]->q = blk->block_cells[i][j].iface[1]->normal;
		scale_vector(blk->block_cells[i][j].iface[1]->q,-1 * qwall[j]); //flux inwards

		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[1]->pos);
		blk->block_cells[i][j].iface[1]->T = (-qwall[j]/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) + blk->block_cells[i][j].T );
	    }
	    break;
	}
    case 2:
	{

	    //Read in time varying terms if flag set
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc("north", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	this->qwall = terms; //"this" is pointer to the current block
	    }

	    int j = blk->nnj-2;
	    for (int i = 0; i<blk->nni-1 ; i++){

		blk->block_cells[i][j].iface[2]->q = blk->block_cells[i][j].iface[2]->normal;
		scale_vector(blk->block_cells[i][j].iface[2]->q,-1 * qwall[i]); //flux inwards

		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[2]->pos);
		blk->block_cells[i][j].iface[2]->T = (-qwall[i]/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) + blk->block_cells[i][j].T );
	    }
	    break;
	}
    case 3:
	{

	    //Read in time varying terms if flag set
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc("west", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	this->qwall = terms; //"this" is pointer to the current block
	    }

	    int i = 0;
	    for (int j = 0; j<blk->nnj-1 ; j++){
     
		blk->block_cells[i][j].iface[3]->q = blk->block_cells[i][j].iface[3]->normal;
		scale_vector(blk->block_cells[i][j].iface[3]->q,qwall[j]);

		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[3]->pos);
		blk->block_cells[i][j].iface[3]->T = -1*(-qwall[j]/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) - blk->block_cells[i][j].T );
	    }
	    break;
	}
    }
}


void BC_USERDEF_Q::print_type(){
    std::cout << "Fixed flux BC user defined vecotr q [W/m2]" <<  std::endl;
}
