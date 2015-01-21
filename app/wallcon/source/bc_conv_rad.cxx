#include <cmath>

#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"
#include "../../../lib/nm/source/ridder.hh"

#include "bc_conv_rad.hh"
#include "solid_cell.hh"

BC_CONV_RAD::BC_CONV_RAD(double eps, double h, double T_inf)
    : eps(eps), h(h), T_inf(T_inf)
{}

void BC_CONV_RAD::apply_bc(SolidBlock *blk, int which_boundary ){
    std::vector<double> length(2,0.0);
    
    double Twall(0),Tcell(0),k(0),distance(0);
    double Tol = 1e-8; //Tolerance for convergence
    double sigma = 5.670373e-8; // W.m^-1.K^-1
    std::function <double (double)> f;
    
    switch (which_boundary){
    case 0: //South boundary
	{
	    int j = 0;
	    for (int i = 0; i<blk->nni-1 ; i++){
		//Constants to be used.
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[0]->pos);
		distance = norm2_vector(length);
		Tcell = blk->block_cells[i][j].T;
		k = blk->block_cells[i][j].k;
		
		//Lambda function for root finding.
		f = [this, sigma, k, Tcell, distance](double T) {
		    return sigma*eps*( pow(T_inf, 4.0) - pow(T, 4) ) + h*(T_inf - T) + k*(Tcell - T)/distance; 
		};
		
		//Solve via ridder search
		Twall = solve(f, Tcell, T_inf, Tol);
		
		//Set properties.
		blk->block_cells[i][j].iface[0]->T = Twall;
		blk->block_cells[i][j].iface[0]->q = blk->block_cells[i][j].iface[0]->normal;
		scale_vector(blk->block_cells[i][j].iface[0]->q, -1*k*(Tcell-Twall)/distance); //Negative here to account reversed normal
	    }
	    break;
	}
    case 1:
	{ 
	    int i = blk->nni-2;
	    for (int j = 0; j<blk->nnj-1 ; j++){
		
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[1]->pos);
		distance = norm2_vector(length);
		Tcell = blk->block_cells[i][j].T;
		k = blk->block_cells[i][j].k;
		
		f = [this, sigma, k, Tcell, distance](double T) {
		    return sigma*eps*( pow(T_inf, 4.0) - pow(T, 4) ) + h* (T_inf - T) + k*(Tcell - T)/distance; 
		};
		
		Twall = solve(f, Tcell, T_inf, Tol);
		blk->block_cells[i][j].iface[1]->T = Twall;
		blk->block_cells[i][j].iface[1]->q = blk->block_cells[i][j].iface[1]->normal;
		scale_vector(blk->block_cells[i][j].iface[1]->q, k*(Tcell-Twall)/distance);
	    }
	    break;
	}
    case 2:
	{
	    int j = blk->nnj-2;
	    for (int i = 0; i<blk->nni-1 ; i++){
		
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[2]->pos);
		distance = norm2_vector(length);
		Tcell = blk->block_cells[i][j].T;
		k = blk->block_cells[i][j].k;
		
		f = [this, sigma, k, Tcell, distance](double T) {
		    return sigma*eps*( pow(T_inf, 4.0) - pow(T, 4) ) + h* (T_inf - T) + k*(Tcell - T)/distance; 
		};
		
		Twall = solve(f, Tcell, T_inf, Tol);
		blk->block_cells[i][j].iface[2]->T = Twall;
		blk->block_cells[i][j].iface[2]->q = blk->block_cells[i][j].iface[2]->normal;
		scale_vector(blk->block_cells[i][j].iface[2]->q, k*(Tcell-Twall)/distance);
	    }	
	    break;
	}
    case 3:
	{
	    int i = 0;
	    for (int j = 0; j<blk->nnj-1 ; j++){
		
		subtract_vectors(length, blk->block_cells[i][j].pos, blk->block_cells[i][j].iface[3]->pos);
		distance = norm2_vector(length);
		Tcell = blk->block_cells[i][j].T;
		k = blk->block_cells[i][j].k;
		
		f = [this, sigma, k, Tcell, distance](double T) {
		    return sigma*eps*( pow(T_inf, 4.0) - pow(T, 4) ) + h* (T_inf - T) + k*(Tcell - T)/distance; 
		};
		
		Twall = solve(f, Tcell, T_inf, Tol);
		blk->block_cells[i][j].iface[3]->T = Twall;
		blk->block_cells[i][j].iface[3]->q = blk->block_cells[i][j].iface[3]->normal;
		scale_vector(blk->block_cells[i][j].iface[3]->q, -1*k*(Tcell-Twall)/distance);
	    }
	    break;
	}
    }
}


void BC_CONV_RAD::print_type(){
    std::cout << "Combined Convection-Radiation BC: emissivity = "<< eps << "h" << h << ", T_inf = " << T_inf << "K" <<  std::endl;
}

