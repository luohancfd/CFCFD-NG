
#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"

#include "solid_block.hh"
#include "solid_bc.hh"
#include <cstdlib>


//Use boundary conditions and cell temps to get fluxes

//This one is obslete as systems of eqns are requried in the bcs
void SolidBlock::OLD_update_boundary_conditions(){ //this may have to change to temperature gradiens for non-isotropic blocks (ie. cell k not unique)
    //South
    { int j = 0;
        for (int i = 0; i<nni-1 ; i++){
            block_bc[0]->apply_bc(this, i,j,0); //"this" is pointer to the current block
        }
    }
    
    //EAST
    { int i = nni-2;
        for (int j = 0; j<nnj-1 ; j++){
            block_bc[1]->apply_bc(this,i, j, 1);
        }
    }
    
    //North
    { int j = nnj-2;
        for (int i = 0; i<nni-1 ; i++){
            block_bc[2]->apply_bc(this, i, j, 2);
        }
    }
    
    //West
    { int i = 0;
        for (int j = 0; j<nnj-1 ; j++){
            block_bc[3]->apply_bc(this, i, j, 3);
        }
    }
}


void SolidBlock::update_boundary_conditions(){ //this may have to change to temperature gradiens for non-isotropic blocks (ie. cell k not unique)

    //South
    block_bc[0]->apply_bc(this, 0); //"this" is pointer to the current block
    
    //EAST
    block_bc[1]->apply_bc(this, 1);
    
    //North
    block_bc[2]->apply_bc(this, 2);
    
    //West
    block_bc[3]->apply_bc(this, 3);
}


void SolidBlock::update_internal_secondary_interface_temperatures(){
    for (int j = 1; j<nnj-1; ++j) {
        for (int i = 1; i<nni-1; ++i){
	    
            //Get interface temperatures for update
            //this does not includes corners
            block_vertices[i][j].iface[0]->T = 0.5 * (block_cells[i-1][j-1].T + block_cells[i][j-1].T); //south
            block_vertices[i][j].iface[1]->T = 0.5 * (block_cells[i][j-1].T + block_cells[i][j].T); //East
            block_vertices[i][j].iface[2]->T = 0.5 * (block_cells[i][j].T + block_cells[i-1][j].T); //North
            block_vertices[i][j].iface[3]->T = 0.5 * (block_cells[i-1][j].T + block_cells[i-1][j-1].T); //West
        }
    }
}

void SolidBlock::update_boundary_secondary_interface_temperatures(){
    //South
    { int j = 0;
        for (int i = 1; i<nni-1 ; i++){
            //perpendicular to boundary
            block_vertices[i][j].iface[1]->T = 0.5 * (block_cells[i][j].iface[0]->T + block_cells[i][j].T); //
            block_vertices[i][j].iface[3]->T = 0.5 * (block_cells[i-1][j].iface[0]->T + block_cells[i-1][j].T); //
	    
            //on boundary
            block_vertices[i][j].iface[0]->T = 0.5 * (block_cells[i-1][j].iface[0]->T + block_cells[i][j].iface[0]->T);
        }
    }
    
    //EAST
    { int i = nni-1;
        for (int j = 1; j<nnj-1 ; j++){
            block_vertices[i][j].iface[0]->T = 0.5 * (block_cells[i-1][j-1].iface[1]->T + block_cells[i-1][j-1].T); //
            block_vertices[i][j].iface[2]->T = 0.5 * (block_cells[i-1][j].iface[1]->T + block_cells[i-1][j].T); //
	    
            //on boundary
            block_vertices[i][j].iface[1]->T = 0.5 * (block_cells[i-1][j-1].iface[1]->T + block_cells[i-1][j].iface[1]->T);
        }
    }
    
    //North
    { int j = nnj-1;
        for (int i = 1; i<nni-1 ; i++){
            block_vertices[i][j].iface[1]->T = 0.5 * (block_cells[i][j-1].iface[2]->T + block_cells[i][j-1].T); //
            block_vertices[i][j].iface[3]->T = 0.5 * (block_cells[i-1][j-1].iface[2]->T + block_cells[i-1][j-1].T); //
	    
            //on boundary
            block_vertices[i][j].iface[2]->T = 0.5 * (block_cells[i-1][j-1].iface[2]->T + block_cells[i][j-1].iface[2]->T);
        }
    }
    
    //West
    { int i = 0;
        for (int j = 1; j<nnj-1 ; j++){
            block_vertices[i][j].iface[0]->T = 0.5 * (block_cells[i][j-1].iface[3]->T + block_cells[i][j-1].T); //
            block_vertices[i][j].iface[2]->T = 0.5 * (block_cells[i][j].iface[3]->T + block_cells[i][j].T); //
	    
            //on boundary
            block_vertices[i][j].iface[3]->T = 0.5 * (block_cells[i][j-1].iface[3]->T + block_cells[i][j].iface[3]->T);
        }
    }
}

void SolidBlock::update_vertex_derivatives(){
    for (int j = 0; j<nnj; ++j) {
        for (int i = 0; i<nni; ++i){
	    
            //T . n * dS
            //Store in secondary iface fluxes
            //block corners don't matter see below
	    //Uses vertex interface (ie secondary interface) flux pointer to store grad(T).
            block_vertices[i][j].iface[0]->q = block_vertices[i][j].iface[0]->normal; //set vecotr
            scale_vector(block_vertices[i][j].iface[0]->q, -1.0*block_vertices[i][j].iface[0]->T * block_vertices[i][j].iface[0]->length); // T*n as T scalar

            block_vertices[i][j].iface[1]->q = block_vertices[i][j].iface[1]->normal; //set vecotr
            scale_vector(block_vertices[i][j].iface[1]->q, block_vertices[i][j].iface[1]->T * block_vertices[i][j].iface[1]->length); // T*n as T scalar

            //make normal vector negative so to point outwards
            block_vertices[i][j].iface[2]->q = block_vertices[i][j].iface[2]->normal; //set vecotr
            scale_vector(block_vertices[i][j].iface[2]->q, block_vertices[i][j].iface[2]->T * block_vertices[i][j].iface[2]->length); // T*n as T scalar

            //make normal vector negative so to point outwards
            block_vertices[i][j].iface[3]->q = block_vertices[i][j].iface[3]->normal; //set vecotr
            scale_vector(block_vertices[i][j].iface[3]->q, -1.0*block_vertices[i][j].iface[3]->T * block_vertices[i][j].iface[3]->length); // T*n as T scalar

            //Sum ( T . n * dS )
            block_vertices[i][j].grad_T[0] = 0.0;
            block_vertices[i][j].grad_T[1] = 0.0;
            add_vectors(block_vertices[i][j].grad_T, block_vertices[i][j].iface[0]->q, block_vertices[i][j].iface[1]->q);
            add_vectors(block_vertices[i][j].grad_T, block_vertices[i][j].grad_T, block_vertices[i][j].iface[2]->q);
            add_vectors(block_vertices[i][j].grad_T, block_vertices[i][j].grad_T, block_vertices[i][j].iface[3]->q);

            // 1/V * ( Sum ( T . n * dS ) )
            scale_vector(block_vertices[i][j].grad_T, 1.0/block_vertices[i][j].volume);

        }
    }
}

void SolidBlock::update_interface_fluxes(){ 
    //Primary interfaces using average of vertex temperature derivatives
    //Don't update boundary interfaces only internal and perpendicular to boundary
    //Puts q in direction  of temperature gradient
    //Not this uses block properties. If want to use cell properties might have to average them.

    //Eastward facing normals

    if (anisotropic_flag == 0){
	for (int j=0; j<(nnj-1); j++){
	    for (int i=1; i<nni-1; i++){
		add_vectors(block_iface_i[i][j].q, block_vertices[i][j].grad_T, block_vertices[i][j+1].grad_T);
		scale_vector(block_iface_i[i][j].q, -k*0.5); //average
	    }
	}
	
	//Northward facing normals
	for (int j=1;j<(nnj-1);j++){
	    for (int i=0;i<(nni-1);i++){
		add_vectors(block_iface_j[i][j].q, block_vertices[i][j].grad_T, block_vertices[i+1][j].grad_T);
		scale_vector(block_iface_j[i][j].q, -k*0.5); //average
	    }
	}
    }
    else if (anisotropic_flag == 1){
	std::vector<double> gradT(2,0.0);
	std::vector<double> fluxi(2,0.0), fluxj(2,0.0);
	double gradTn, gradTt;

	//Eastward Normals
	for (int j=0; j<(nnj-1); j++){
	    for (int i=1; i<nni-1; i++){
		//Average grad(T)
		add_vectors(gradT, block_vertices[i][j].grad_T, block_vertices[i][j+1].grad_T);
		scale_vector(gradT, 0.5); //average

		//Resolve gradT to interface n t directions direction
		gradTn = gradT[0]*block_iface_i[i][j].normal[0] + gradT[1]*block_iface_i[i][j].normal[1];
		gradTt = gradT[0]*block_iface_i[i][j].tangent[0] + gradT[1]*block_iface_i[i][j].tangent[1];

		//normal is the cell i direction. tangent is cell j direction
		// -qi = k11*dT/di + k12*dT/dj    -qj = k12*dT/di + k22*dT/dj
		fluxi = block_iface_i[i][j].normal;
		scale_vector( fluxi, -k11*gradTn -k12*gradTt); //normal flux in x,y convention
		fluxj = block_iface_i[i][j].tangent;
		scale_vector( fluxj, -k12*gradTn -k22*gradTt); //tangential flux in x,y convention

		//Resolve back to
		block_iface_i[i][j].q[0] = fluxi[0] + fluxj[0];
		block_iface_i[i][j].q[1] = fluxi[1] + fluxj[1];

		//Resolve back
	    }
	}
	
	//Northward facing normals
	for (int j=1;j<(nnj-1);j++){
	    for (int i=0;i<(nni-1);i++){
		add_vectors(gradT, block_vertices[i][j].grad_T, block_vertices[i+1][j].grad_T);
		scale_vector(gradT, 0.5); //average

		//Resolve gradT to interface n t directions direction
		gradTn = gradT[0]*block_iface_j[i][j].normal[0] + gradT[1]*block_iface_j[i][j].normal[1];
		gradTt = gradT[0]*block_iface_j[i][j].tangent[0] + gradT[1]*block_iface_j[i][j].tangent[1];

		//normal is the cell j direction. tangent is cell i direction
		// -qi = k11*dT/di + k12*dT/dj    -qj = k12*dT/di + k22*dT/dj
		fluxi = block_iface_j[i][j].tangent;
		scale_vector( fluxi, -k11*gradTt -k12*gradTn); //tangential flux in x,y convention
		fluxj = block_iface_j[i][j].normal;
		scale_vector( fluxj, -k12*gradTt -k22*gradTn); //normal flux in x,y convention

		//Resolve back to
		block_iface_j[i][j].q[0] = fluxi[0] + fluxj[0];
		block_iface_j[i][j].q[1] = fluxi[1] + fluxj[1];
	    }
	}
    }
    else {
	std::cout << "Something went wrong updating intrface fluxes. Bailing out." << std::endl;
	exit(1);
    }
    
}


void SolidBlock::update_cell_energy_derivative(){
    //dondt(e) = (1/V * sum fluxes . n * dS ) + sources
    //flux in is positive => if  q.n < 0 inwards if q.n>0 outwards
    
    for (int j=0; j<(nnj-1); j++){
        for (int i=0; i<(nni-1); i++){
	    
            block_cells[i][j].deondt = block_cells[i][j].source + ( 1.0/block_cells[i][j].volume ) * (
                    -1*block_cells[i][j].iface[0]->length * ( block_cells[i][j].iface[0]->q[0] * -1*block_cells[i][j].iface[0]->normal[0] + block_cells[i][j].iface[0]->q[1] * -1*block_cells[i][j].iface[0]->normal[1])
                    +
                    -1*block_cells[i][j].iface[1]->length * ( block_cells[i][j].iface[1]->q[0] * block_cells[i][j].iface[1]->normal[0] + block_cells[i][j].iface[1]->q[1] * block_cells[i][j].iface[1]->normal[1])
                    +
                    -1*block_cells[i][j].iface[2]->length * ( block_cells[i][j].iface[2]->q[0] * block_cells[i][j].iface[2]->normal[0] + block_cells[i][j].iface[2]->q[1] * block_cells[i][j].iface[2]->normal[1])
                    +
                    -1*block_cells[i][j].iface[3]->length * ( block_cells[i][j].iface[3]->q[0] * -1*block_cells[i][j].iface[3]->normal[0] + block_cells[i][j].iface[3]->q[1] * -1*block_cells[i][j].iface[3]->normal[1])
                    ); //reversed normals for faces 0 and 3.
	    
        }
    }
}

void SolidBlock::time_update(double &dt, int update_scheme){
    //dondt e = sum fluxes . n * dS
    //flux in is positive => if  q.n < 0 inwards if q.n>0 outwards
    
    if (update_scheme == 0){ //Euler
	// std::cout << "Wallcon update is Euler" << std::endl;
	for (int j=0; j<(nnj-1); j++){
	    for (int i=0; i<(nni-1); i++){
	    
		block_cells[i][j].one_stage_update(dt); //Euler update to begine with change this later to be more general
	    }
        }
    }
    else if (update_scheme == 1){ //Euler
	// std::cout << "Wallcon update is Predictor Corrector" << std::endl;
	for (int j=0; j<(nnj-1); j++){
	    for (int i=0; i<(nni-1); i++){	
		block_cells[i][j].e0 = block_cells[i][j].e;
		block_cells[i][j].de0ondt = block_cells[i][j].deondt;
		block_cells[i][j].one_stage_update(dt); 
	    }
        }
	this->update_boundary_conditions();
        this->update_boundary_secondary_interface_temperatures();
        this->update_internal_secondary_interface_temperatures();
        this->update_vertex_derivatives();
        this->update_interface_fluxes();
        this->update_cell_energy_derivative();
	for (int j=0; j<(nnj-1); j++){
	    for (int i=0; i<(nni-1); i++){	
		block_cells[i][j].two_stage_update(dt); 
	    }
        }

	
    }
    else{
	std::cout<<"Something went wrong in SolidBlock::time_update()" << std::endl;
    }
}

