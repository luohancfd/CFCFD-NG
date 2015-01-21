#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"

#include "bc_userdef_t.hh"
#include "solid_cell.hh"
#include "useful_funcs.hh"


BC_USERDEF_T::BC_USERDEF_T(){}

BC_USERDEF_T::BC_USERDEF_T(std::vector<double> Twall)
    : Twall(Twall)
{}

void BC_USERDEF_T::apply_bc(SolidBlock *blk, int which_boundary )
{
    std::vector<double> length(2,0.0);
    //Temperature given in vector Twall
    switch (which_boundary){
    case 0: //South boundary
	{ 

	    //Read in time varying terms if flag set
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc("south", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	this->Twall = terms; //"this" is pointer to the current block
	    }

	    // std::cout<<"South TBC " << Twall[0] << std::endl;

	    int j = 0;
	    for (int i = 0; i<blk->nni-1 ; i++){

		//Set interface temp
		blk->block_cells[i][j].iface[0]->T = Twall[i];
        
		//Set flux direction and solve analytically
		blk->block_cells[i][j].iface[0]->q = blk->block_cells[i][j].iface[0]->normal;
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[0]->pos);
		scale_vector(blk->block_cells[i][j].iface[0]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].T - blk->block_cells[i][j].iface[0]->T ) / norm2_vector( length,0,length.size() ) );
	    }
	    break;
	}
    case 1:
	{

	    //Read in time varying terms if flag set
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc("east", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	this->Twall = terms; //"this" is pointer to the current block
	    }

	    // std::cout<<"East TBC " << Twall[0] << std::endl;


	    int i = blk->nni-2;
	    for (int j = 0; j<blk->nnj-1 ; j++){

		blk->block_cells[i][j].iface[1]->T = Twall[j];

		blk->block_cells[i][j].iface[1]->q = blk->block_cells[i][j].iface[1]->normal;
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[1]->pos);
		scale_vector(blk->block_cells[i][j].iface[1]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].iface[1]->T - blk->block_cells[i][j].T  ) / norm2_vector( length,0,length.size() ) );
	    }
	    break;
	}
    case 2:
	{

	    //Read in time varying terms if flag set
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc("north", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	this->Twall = terms; //"this" is pointer to the current block
	    }

	    // std::cout<<"North TBC " << Twall[0] << std::endl;


	    int j = blk->nnj-2;
	    for (int i = 0; i<blk->nni-1 ; i++){

		blk->block_cells[i][j].iface[2]->T = Twall[i];

		blk->block_cells[i][j].iface[2]->q = blk->block_cells[i][j].iface[2]->normal;
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[2]->pos);
		scale_vector(blk->block_cells[i][j].iface[2]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].iface[2]->T - blk->block_cells[i][j].T ) / norm2_vector( length,0,length.size() ) );
	    }
	    break;
	}
    case 3:
	{

	    //Read in time varying terms if flag set
	    std::vector<double> terms; //hope the scope works
	    if (blk->time_varying_terms_flag == 1){
	    	terms = read_time_varying_bc("west", blk->time_varying_terms_dir, blk->time_varying_iteration);
	    	this->Twall = terms; //"this" is pointer to the current block
	    }
	    // std::cout<<"West TBC " << Twall[0] << std::endl;

	    int i = 0;
	    for (int j = 0; j<blk->nnj-1 ; j++){

		blk->block_cells[i][j].iface[3]->T = Twall[j];

		blk->block_cells[i][j].iface[3]->q = blk->block_cells[i][j].iface[3]->normal;
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[3]->pos);
		scale_vector(blk->block_cells[i][j].iface[3]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].T - blk->block_cells[i][j].iface[3]->T ) / norm2_vector( length,0,length.size() ) );
	    }
	    break;
	}
    }
}


void BC_USERDEF_T::print_type(){
    std::cout << "Fixed temp BC user defined vecotr T [K]" <<  std::endl;
}
