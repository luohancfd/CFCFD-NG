#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"

#include "bc_convection.hh"
#include "solid_cell.hh"

BC_CONVECTION::BC_CONVECTION(double h, double T_inf)
    : h(h), T_inf(T_inf)
{}

void BC_CONVECTION::print_type(){
    std::cout << "Convection BC: h = "<< h << "WK/m, T_inf = " << T_inf << "K" <<  std::endl;
}

void BC_CONVECTION::apply_bc(SolidBlock *blk, int which_boundary ){
    switch (which_boundary){
    case 0: //South boundary
	{
	    std::vector<double> length(2,0.0);
	    double k;
	    double T;
	    
	    int j = 0;
	    for (int i = 0; i<blk->nni-1 ; i++){
		k = blk->block_cells[i][j].k;
		T = blk->block_cells[i][j].T;
		
		//Constants to be used.
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[0]->pos);
		
		//Set interface temperature by solving analytically over convection-conduction boudary.
		blk->block_cells[i][j].iface[0]->T = ( h*T_inf + ( k/norm2_vector(length,0,length.size()) )*T ) / ( h + k/norm2_vector(length,0,length.size()) );
		
		//Set flux by using temperature above and solving flux eqn over first cell.
		blk->block_cells[i][j].iface[0]->q = blk->block_cells[i][j].iface[0]->normal;
		scale_vector(blk->block_cells[i][j].iface[0]->q, -1 * k * ( T - blk->block_cells[i][j].iface[0]->T ) / norm2_vector( length,0,length.size() ) );
	    }
	    break;
	}
    case 1:
	{
	    std::vector<double> length(2,0.0);
	    double k;
	    double T;

	    int i = blk->nni-2;
	    for (int j = 0; j<blk->nnj-1 ; j++){
		k = blk->block_cells[i][j].k;
		T = blk->block_cells[i][j].T;
		
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[1]->pos);
		
		blk->block_cells[i][j].iface[1]->T = ( h*T_inf + ( k/norm2_vector(length,0,length.size()) )*T ) / ( h + k/norm2_vector(length,0,length.size()) );
		
		blk->block_cells[i][j].iface[1]->q = blk->block_cells[i][j].iface[1]->normal;
		scale_vector(blk->block_cells[i][j].iface[1]->q, -1 * k * ( blk->block_cells[i][j].iface[1]->T - T  ) / norm2_vector( length,0,length.size() ) );
	    }
	    break;
	}
    case 2:
	{ 
	    std::vector<double> length(2,0.0);
	    double k;
	    double T;

	    int j = blk->nnj-2;
	    for (int i = 0; i<blk->nni-1 ; i++){
	        k = blk->block_cells[i][j].k;
		T = blk->block_cells[i][j].T;
		
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[2]->pos);
		
		blk->block_cells[i][j].iface[2]->T = ( h*T_inf + ( k/norm2_vector(length,0,length.size()) )*T ) / ( h + k/norm2_vector(length,0,length.size()) );
		
		blk->block_cells[i][j].iface[2]->q = blk->block_cells[i][j].iface[2]->normal;
		scale_vector(blk->block_cells[i][j].iface[2]->q, -1 * k * ( blk->block_cells[i][j].iface[2]->T - T ) / norm2_vector( length,0,length.size() ) );
	    }
	    break;
	}
    case 3:
	{
	    std::vector<double> length(2,0.0);
	    double k;
	    double T;

	    int i = 0;
	    for (int j = 0; j<blk->nnj-1 ; j++){
		k = blk->block_cells[i][j].k;
		T = blk->block_cells[i][j].T;
		
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[3]->pos);
		
		blk->block_cells[i][j].iface[3]->T = ( h*T_inf + ( k/norm2_vector(length,0,length.size()) )*T ) / ( h + k/norm2_vector(length,0,length.size()) );
		
		blk->block_cells[i][j].iface[3]->q = blk->block_cells[i][j].iface[3]->normal;
		scale_vector(blk->block_cells[i][j].iface[3]->q, -1 * k * ( T - blk->block_cells[i][j].iface[3]->T ) / norm2_vector( length,0,length.size() ) );
	    }
	    break;
	}
    }
}
