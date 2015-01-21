#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"
#include <math.h>
#include "bc_fixed_t.hh"
#include "solid_cell.hh"


BC_FIXED_T::BC_FIXED_T(double Twall)
    : Twall(Twall)
{}

void BC_FIXED_T::apply_bc(SolidBlock *blk, int which_boundary){ 
    std::vector<double> length(2,0.0);
    switch (which_boundary){
    case 0: //South boundary
	if (blk->anisotropic_flag == 0){
	    int j = 0;
	    for (int i = 0; i<blk->nni-1 ; i++){

		//Set known temp.
		blk->block_cells[i][j].iface[0]->T = Twall;
	    
		//Set flux direction.
		blk->block_cells[i][j].iface[0]->q = blk->block_cells[i][j].iface[0]->normal;
	    
		//Solve for flux analytically and scale.
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[0]->pos);
		scale_vector(blk->block_cells[i][j].iface[0]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].T - blk->block_cells[i][j].iface[0]->T ) / norm2_vector( length,0,length.size() ) );
	    }
	}
	else if (blk->anisotropic_flag == 1){
	    double dT, dL, q;
	    int j = 0;
	    for (int i = 0; i<blk->nni-1 ; i++){

		dL = (blk->block_cells[i][j].pos[0] - blk->block_cells[i][j].iface[0]->pos[0])*(blk->block_cells[i][j].pos[0] - blk->block_cells[i][j].iface[0]->pos[0]);
		dL += (blk->block_cells[i][j].pos[1] - blk->block_cells[i][j].iface[0]->pos[1])*(blk->block_cells[i][j].pos[1] - blk->block_cells[i][j].iface[0]->pos[1]);
		dL = pow(dL, 0.5);

		//Set known temp.
		blk->block_cells[i][j].iface[0]->T = Twall;
		//Set dT
		dT = blk->block_cells[i][j].T - Twall;

		// k22 coincides with interface normal direction
		q = -blk->block_cells[i][j].k22 * dT/dL;
		    
		//set direction of flux and scale
		blk->block_cells[i][j].iface[0]->q = blk->block_cells[i][j].iface[0]->normal;
		scale_vector( blk->block_cells[i][j].iface[0]->q, q);
	    }
	}
	else {
	    std::cout << "Something went wrong with south boundary condition.\nBailing Out."<<std::endl;
	    exit(1);
	}
	break;
    case 1:
	if (blk->anisotropic_flag == 0){
	    int i = blk->nni-2;
	    for (int j = 0; j<blk->nnj-1 ; j++){

		blk->block_cells[i][j].iface[1]->T = Twall;

		blk->block_cells[i][j].iface[1]->q = blk->block_cells[i][j].iface[1]->normal;

		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[1]->pos);
		scale_vector(blk->block_cells[i][j].iface[1]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].iface[1]->T - blk->block_cells[i][j].T  ) / norm2_vector( length,0,length.size() ) );
	    }
	}
	else if (blk->anisotropic_flag == 1){
	    double dT, dL, q;
	    int i = blk->nni-2;
	    for (int j = 0; j<blk->nnj-1 ; j++){

		dL = (blk->block_cells[i][j].pos[0] - blk->block_cells[i][j].iface[1]->pos[0])*(blk->block_cells[i][j].pos[0] - blk->block_cells[i][j].iface[1]->pos[0]);
		dL += (blk->block_cells[i][j].pos[1] - blk->block_cells[i][j].iface[1]->pos[1])*(blk->block_cells[i][j].pos[1] - blk->block_cells[i][j].iface[1]->pos[1]);
		dL = pow(dL, 0.5);

		//Set known temp.
		blk->block_cells[i][j].iface[1]->T = Twall;
		//Set dT
		dT = Twall - blk->block_cells[i][j].T; //Due to direction of normal

		// k22 coincides with interface normal direction
		q = -blk->block_cells[i][j].k11 * dT/dL;
		    
		//set direction of flux and scale
		blk->block_cells[i][j].iface[1]->q = blk->block_cells[i][j].iface[1]->normal;
		scale_vector( blk->block_cells[i][j].iface[1]->q, q);
	
	    }
	}
	else {
	    std::cout << "Something went wrong with east boundary condition.\nBailing Out."<<std::endl; 
	    exit(1);
	}
	break;
    case 2:
	if (blk->anisotropic_flag == 0){
	    int j = blk->nnj-2;
	    for (int i = 0; i<blk->nni-1 ; i++){

		blk->block_cells[i][j].iface[2]->T = Twall;

		blk->block_cells[i][j].iface[2]->q = blk->block_cells[i][j].iface[2]->normal;

		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[2]->pos);
		scale_vector(blk->block_cells[i][j].iface[2]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].iface[2]->T - blk->block_cells[i][j].T ) / norm2_vector( length,0,length.size() ) );
	    }
	}
	else if (blk-> anisotropic_flag == 1){
	    double dT, dL, q;
	    int j = blk->nnj-2;
	    for (int i = 0; i<blk->nni-1 ; i++){

		dL = (blk->block_cells[i][j].pos[0] - blk->block_cells[i][j].iface[2]->pos[0])*(blk->block_cells[i][j].pos[0] - blk->block_cells[i][j].iface[2]->pos[0]);
		dL += (blk->block_cells[i][j].pos[1] - blk->block_cells[i][j].iface[2]->pos[1])*(blk->block_cells[i][j].pos[1] - blk->block_cells[i][j].iface[2]->pos[1]);
		dL = pow(dL, 0.5);

		//Set known temp.
		blk->block_cells[i][j].iface[2]->T = Twall;
		//Set dT
		dT = Twall - blk->block_cells[i][j].T; //Due to direction of normal

		// k22 coincides with interface normal direction
		q = -blk->block_cells[i][j].k22 * dT/dL;
		    
		//set direction of flux and scale
		blk->block_cells[i][j].iface[2]->q = blk->block_cells[i][j].iface[2]->normal;
		scale_vector( blk->block_cells[i][j].iface[2]->q, q);
	    }
	}
	else{
	    std::cout << "Something went wrong with north boundary condition.\nBailing Out."<<std::endl;
	    exit(1);
	}
	break;
    case 3:
	if (blk->anisotropic_flag == 0) {
	    int i = 0;
	    for (int j = 0; j<blk->nnj-1 ; j++){

		blk->block_cells[i][j].iface[3]->T = Twall;

		blk->block_cells[i][j].iface[3]->q = blk->block_cells[i][j].iface[3]->normal;

		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[3]->pos);
		scale_vector(blk->block_cells[i][j].iface[3]->q, -1 * blk->block_cells[i][j].k * ( blk->block_cells[i][j].T - blk->block_cells[i][j].iface[3]->T ) / norm2_vector( length,0,length.size() ) );
	    }
	}
	else if (blk->anisotropic_flag == 1){
	    double dT, dL, q;
	    int i = 0;
	    for (int j = 0; j<blk->nnj-1 ; j++){
		dL = (blk->block_cells[i][j].pos[0] - blk->block_cells[i][j].iface[3]->pos[0])*(blk->block_cells[i][j].pos[0] - blk->block_cells[i][j].iface[3]->pos[0]);
		dL += (blk->block_cells[i][j].pos[1] - blk->block_cells[i][j].iface[3]->pos[1])*(blk->block_cells[i][j].pos[1] - blk->block_cells[i][j].iface[3]->pos[1]);
		dL = pow(dL, 0.5);

		//Set known temp.
		blk->block_cells[i][j].iface[3]->T = Twall;
		//Set dT
		dT = blk->block_cells[i][j].T - Twall;

		// k22 coincides with interface normal direction
		q = -blk->block_cells[i][j].k11 * dT/dL;
		    
		//set direction of flux and scale
		blk->block_cells[i][j].iface[3]->q = blk->block_cells[i][j].iface[3]->normal;
		scale_vector( blk->block_cells[i][j].iface[3]->q, q);
	    }
	}
	else{
	    std::cout << "Something went wrong with west boundary condition.\nBailing Out."<<std::endl;
	    exit(1);
	}
	break;
    }
}

// void BC_FIXED_T::apply_bc(Solid_Anisotropic_Block *blk, int i, int j, int which_boundary){ 
//     std::vector<double> dx(2,0.0);
//     double dT;

//     // -qx = k11*dTx k12*dTy

//     switch (which_boundary){
//     case 0: //South boundary
// 	//Set known temp.
//         blk->block_cells[i][j].iface[0]->T = Twall;
// 	//Set dT
// 	dT = blk->block_cells[i][j].T - Twall;
// 	// qx
// 	blk->block_cells[i][j].iface[0]->q[0] = -blk->block_cells[i][j].k11*(dT/dx[0]) - blk->block_cells[i][j].k12*(dT/dx[1]);
// 	// qy
// 	blk->block_cells[i][j].iface[0]->q[1] = -blk->block_cells[i][j].k12*(dT/dx[0]) - blk->block_cells[i][j].k22*(dT/dx[1]);
//         break;
//     case 1:
// 	std::cout << "bc_temp_aniso not yet" << std::endl;
//         break;
//     case 2:
// 	std::cout << "bc_temp_aniso not yet" << std::endl;
//         break;
//     case 3:
// 	std::cout << "bc_temp_aniso not yet" << std::endl;
//         break;
//     }
// }

void BC_FIXED_T::print_type(){
    std::cout << "Fixed Temp BC: Twall = " << Twall << "K" <<  std::endl;
}
