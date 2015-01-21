#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"

#include "bc_timevarying.hh"
#include "solid_cell.hh"
#include "useful_funcs.hh"


BC_TIMEVARYING::BC_TIMEVARYING(){}

BC_TIMEVARYING::BC_TIMEVARYING(std::vector<double> qwall)
    : qwall(qwall)
{}

void BC_TIMEVARYING::apply_bc(SolidBlock *blk, int which_boundary )
{
    std::vector<double> length(2,0.0);
    std::vector<double> Twall;
    //Flux defined in vector qwall.
    switch (which_boundary){
    case 0: //South boundary
	{
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc2("T","south", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	Twall = terms; //"this" is pointer to the current block
	    	terms = read_time_varying_bc2("q","south", blk->time_varying_terms_dir, blk->time_varying_iteration);
		this->qwall = terms;
	    }

	    //This would only work for orthogonal standard grids 
	    int j = 0;
	    for (int i = 0; i<blk->nni-1 ; i++){

		//Set flux direction and scale
		blk->block_cells[i][j].iface[0]->q = blk->block_cells[i][j].iface[0]->normal;
		scale_vector(blk->block_cells[i][j].iface[0]->q,qwall[i]);
		blk->block_cells[i][j].iface[0]->T = Twall[i];

	    }
	    break;
	}
    case 1:
	{
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc2("T","east", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	Twall = terms; //"this" is pointer to the current block
	    	terms = read_time_varying_bc2("q","east", blk->time_varying_terms_dir, blk->time_varying_iteration);
		this->qwall = terms;
	    }
	    int i = blk->nni-2;
	    for (int j = 0; j<blk->nnj-1 ; j++){
		blk->block_cells[i][j].iface[1]->q = blk->block_cells[i][j].iface[1]->normal;
		scale_vector(blk->block_cells[i][j].iface[1]->q,qwall[j]); //flux inwards
		blk->block_cells[i][j].iface[1]->T = Twall[j];
	    }
	    break;
	}
    case 2:
	{
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc2("T","north", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	Twall = terms; //"this" is pointer to the current block
	    	terms = read_time_varying_bc2("q","north", blk->time_varying_terms_dir, blk->time_varying_iteration);
		this->qwall = terms;
	    }

	    int j = blk->nnj-2;
	    for (int i = 0; i<blk->nni-1 ; i++){

		blk->block_cells[i][j].iface[2]->q = blk->block_cells[i][j].iface[2]->normal;
		scale_vector(blk->block_cells[i][j].iface[2]->q, qwall[i]); //flux inwards

		blk->block_cells[i][j].iface[2]->T = Twall[i];
	    }
	    break;
	}
    case 3:
	{
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc2("T","west", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	Twall = terms; //"this" is pointer to the current block
	    	terms = read_time_varying_bc2("q","west", blk->time_varying_terms_dir, blk->time_varying_iteration);
		this->qwall = terms;
	    }

	    int i = 0;
	    for (int j = 0; j<blk->nnj-1 ; j++){
     
		blk->block_cells[i][j].iface[3]->q = blk->block_cells[i][j].iface[3]->normal;
		scale_vector(blk->block_cells[i][j].iface[3]->q,qwall[j]);

		blk->block_cells[i][j].iface[3]->T = Twall[j];
	    }
	    break;
	}
    }
}


void BC_TIMEVARYING::print_type(){
    std::cout << "Time varying BC flux and temp [W/m2]" <<  std::endl;
}
