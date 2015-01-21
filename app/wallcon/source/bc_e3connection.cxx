#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"

#include "bc_e3connection.hh"
#include "solid_cell.hh"

BC_E3CONNECTION::BC_E3CONNECTION()
{}


BC_E3CONNECTION::BC_E3CONNECTION(std::vector<double> connection_vector)
    : connection_vector(connection_vector)
{}

void BC_E3CONNECTION::apply_bc(SolidBlock *blk, int which_boundary )
{
    std::vector<double> length(2,0.0);
    //set secondary boundary interface temperatures
    switch (which_boundary){
    case 0: //South boundary
       	
	switch(blk->e3connection_type_flag[0]){
	case 0:
	    {
		int j = 0;
		for (int i = 0; i<blk->nni-1 ; i++){

	    //Set flux direction and scale
	    blk->block_cells[i][j].iface[0]->q = blk->block_cells[i][j].iface[0]->normal;
	    scale_vector(blk->block_cells[i][j].iface[0]->q,connection_vector[i]);
	    
	    //Determine interface temperature by solving conduction analytically
	    subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[0]->pos);
	    blk->block_cells[i][j].iface[0]->T = -1*(-connection_vector[i]/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) - blk->block_cells[i][j].T );
		}
		break;
	    }
	
	case 1:
	    {
		int j = 0;
		for (int i = 0; i<blk->nni-1 ; i++){
		    
		    //Set interface temp
		    blk->block_cells[i][j].iface[0]->T = connection_vector[i];
		    
		    //Set flux direction and solve analytically
		    blk->block_cells[i][j].iface[0]->q = blk->block_cells[i][j].iface[0]->normal;
		    subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[0]->pos);
		    scale_vector(blk->block_cells[i][j].iface[0]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].T - blk->block_cells[i][j].iface[0]->T ) / norm2_vector( length,0,length.size() ) );
		}
		break;
	    }
	}
	break;
	
    case 1:
	
	switch (blk->e3connection_type_flag[1]){
	case 0:
	    {
		int i = blk->nni-2;
		for (int j = 0; j<blk->nnj-1 ; j++){
		    
		    blk->block_cells[i][j].iface[1]->q = blk->block_cells[i][j].iface[1]->normal;
		    scale_vector(blk->block_cells[i][j].iface[1]->q,-1 * connection_vector[j]); //flux inwards
		    
		    subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[1]->pos);
		    blk->block_cells[i][j].iface[1]->T = (-connection_vector[j]/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) + blk->block_cells[i][j].T );
		}
		break;
	    }
	
	case 1:
	    {
		int i = blk->nni-2;
		for (int j = 0; j<blk->nnj-1 ; j++){
		    
		    blk->block_cells[i][j].iface[1]->T =  connection_vector[j];
		    
		    blk->block_cells[i][j].iface[1]->q = blk->block_cells[i][j].iface[1]->normal;
		    subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[1]->pos);
		    scale_vector(blk->block_cells[i][j].iface[1]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].iface[1]->T - blk->block_cells[i][j].T  ) / norm2_vector( length,0,length.size() ) );
		}
		break;
	    }
	}
	break;

    case 2:
	
	switch (blk->e3connection_type_flag[2]){
	case 0:
	    {
		int j = blk->nnj-2;
		for (int i = 0; i<blk->nni-1 ; i++){
		    
		    blk->block_cells[i][j].iface[2]->q = blk->block_cells[i][j].iface[2]->normal;
		    scale_vector(blk->block_cells[i][j].iface[2]->q,-1 * connection_vector[i]); //flux inwards
		    
		    subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[2]->pos);
		    blk->block_cells[i][j].iface[2]->T = (-connection_vector[i]/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) + blk->block_cells[i][j].T );
		}
		break;
	    }

	case 1:
	    {
		int j = blk->nnj-2;
		for (int i = 0; i<blk->nni-1 ; i++){
		    
		    blk->block_cells[i][j].iface[2]->T =  connection_vector[i];
		    
		    blk->block_cells[i][j].iface[2]->q = blk->block_cells[i][j].iface[2]->normal;
		    subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[2]->pos);
		    scale_vector(blk->block_cells[i][j].iface[2]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].iface[2]->T - blk->block_cells[i][j].T ) / norm2_vector( length,0,length.size() ) );
		}
		break;
	    }
	}
	break;

    case 3:
	
	switch(blk->e3connection_type_flag[3]){
	case 0:
	    {
		int i = 0;
		for (int j = 0; j<blk->nnj-1 ; j++){
		    
		    blk->block_cells[i][j].iface[3]->q = blk->block_cells[i][j].iface[3]->normal;
		    scale_vector(blk->block_cells[i][j].iface[3]->q,connection_vector[j]);
		    
		    subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[3]->pos);
		    blk->block_cells[i][j].iface[3]->T = -1*(-connection_vector[j]/blk->block_cells[i][j].k * norm2_vector( length,0,length.size() ) - blk->block_cells[i][j].T );
		}
		break;
	    }

	case 1:
	    {
		int i = 0;
		for (int j = 0; j<blk->nnj-1 ; j++){
		    
		    blk->block_cells[i][j].iface[3]->T =  connection_vector[j];
		    
		    blk->block_cells[i][j].iface[3]->q = blk->block_cells[i][j].iface[3]->normal;
		    subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[3]->pos);
		    scale_vector(blk->block_cells[i][j].iface[3]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].T - blk->block_cells[i][j].iface[3]->T ) / norm2_vector( length,0,length.size() ) );
		}
		break;
	    }
	}
	break;
    }
}
    

void BC_E3CONNECTION::print_type() {
    std::cout << "Eilmer3 connection BC" << std::endl;
    std::cout << "connect_vector = [ ";
    for (size_t i =0; i < connection_vector.size() ; i++){
	std::cout << connection_vector[i] << ", ";
    }
    std::cout << "]" <<std::endl;

}
