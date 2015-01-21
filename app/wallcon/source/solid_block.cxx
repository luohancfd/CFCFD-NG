#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"

#include <cmath>
#include "solid_block.hh"

using namespace std;

SolidBlock::SolidBlock()
    :anisotropic_flag(0), time_varying_terms_flag(0),time_varying_iteration(0),
     e3connection_type_flag(4,0)
     
{}

//allocate memory to arrays in block
//block nni is number of nodes starting at 1 to n in i direction etc //chech memorey allocations as vectors are dynamic and could increas later on
void SolidBlock::allocate_memory(){
    //Verticies
    for (int i=0 ; i<nni ; i++){
        block_vertices.push_back(vector<Solid_FV_Vertex>(nnj));
    }
    
    //Interfaces
    //Eastward facing normals
    for (int i=0 ; i<nni ; i++){
        block_iface_i.push_back(vector<Solid_FV_Interface>(nnj-1));
    }
    //Northward facing normals
    for (int i=0 ; i<nni-1 ; i++){
        block_iface_j.push_back(vector<Solid_FV_Interface>(nnj));
    }
    
    //Cells
    for (int i=0 ; i<nni-1 ; i++){
        block_cells.push_back(vector<Solid_FV_Cell>(nnj-1));
    }
    
    //Secondary Objects
    //Interfaces
    //Eastward facing normals
    for (int i=0 ; i<nni+1 ; i++){
        block_secondary_iface_i.push_back(vector<Solid_FV_Interface>(nnj));
    }
    //Northward facing normals
    for (int i=0 ; i<nni ; i++){
        block_secondary_iface_j.push_back(vector<Solid_FV_Interface>(nnj+1));
    }

}


//Assign cells adding verticies
void SolidBlock::assign_cells_to_block(){ //assigns cells to block after verticies have been created
    std::vector<double> mid1(2,0.0), mid2(2,0.0), centroid(2,0.0);

    for (int j=0; j<(nnj-1); j++){
        for (int i=0; i<(nni-1); i++){
            //verticies
            block_cells[i][j].vrtx[0] = &block_vertices[i][j]; //vertex numbering start in lower left and move anti-clockwise
            block_cells[i][j].vrtx[1] = &block_vertices[i+1][j];
            block_cells[i][j].vrtx[2] = &block_vertices[i+1][j+1];
            block_cells[i][j].vrtx[3] = &block_vertices[i][j+1];
	    
            //Centroid position
            add_vectors(mid1,block_cells[i][j].vrtx[0]->pos,block_cells[i][j].vrtx[2]->pos);
            scale_vector(mid1,0.5); //midpoint diagonal 1
            add_vectors(mid2,block_cells[i][j].vrtx[1]->pos,block_cells[i][j].vrtx[3]->pos);
            scale_vector(mid2,0.5); //midpoint diagonal 2
            add_vectors(centroid,mid1,mid2);
            scale_vector(centroid,0.5); //centroid of shape
            block_cells[i][j].pos = centroid;

            //Area aka volume
            block_cells[i][j].volume = 0.5*fabs( (block_cells[i][j].vrtx[0]->pos[0]-block_cells[i][j].vrtx[2]->pos[0]) * (block_cells[i][j].vrtx[1]->pos[1]-block_cells[i][j].vrtx[3]->pos[1]) - (block_cells[i][j].vrtx[1]->pos[0]-block_cells[i][j].vrtx[3]->pos[0]) * (block_cells[i][j].vrtx[0]->pos[1]-block_cells[i][j].vrtx[2]->pos[1]) );
        }
    }
}


//create interfaces
//horizontal 1st then vertical, normals such that south and west need to be reversed
void SolidBlock::assign_ifaces_to_block(){
    
    std::vector<double> temp_vector(2,0.0),temp_vector1(2,0.0);
    
    //Eastward facing normals
    for (int j=0; j<(nnj-1); j++){
        for (int i=0; i<nni; i++){
	    
            //iface position
            add_vectors(block_iface_i[i][j].pos,block_vertices[i][j].pos,block_vertices[i][j+1].pos);
            scale_vector(block_iface_i[i][j].pos, 1/2.0);

            //iface length
            subtract_vectors(temp_vector,block_vertices[i][j+1].pos,block_vertices[i][j].pos); //Causes negative zero. This may affect error?
            block_iface_i[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());

            //iface normal
            temp_vector1[0] = temp_vector[1];
            temp_vector1[1] = -1*temp_vector[0];
            scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
            block_iface_i[i][j].normal = temp_vector1;

	    //iface tangent
	    block_iface_i[i][j].tangent[0] = -1*block_iface_i[i][j].normal[1];
	    block_iface_i[i][j].tangent[1] = block_iface_i[i][j].normal[0];

        }
    }

    //Northward facing normals
    for (int j=0;j<(nnj);j++){
        for (int i=0;i<(nni-1);i++){
	    
            //iface position
            add_vectors(block_iface_j[i][j].pos,block_vertices[i][j].pos,block_vertices[i+1][j].pos);
            scale_vector(block_iface_j[i][j].pos, 1/2.0);
	    
            //iface length
            subtract_vectors(temp_vector,block_vertices[i][j].pos,block_vertices[i+1][j].pos); //use negative vector so cross product results in +ve normal north
            block_iface_j[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	    
            //iface normal
            temp_vector1[0] = temp_vector[1];
            temp_vector1[1] = -1*temp_vector[0];
            scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
            block_iface_j[i][j].normal = temp_vector1;

	    //iface tangent
	    block_iface_j[i][j].tangent[0] = block_iface_j[i][j].normal[1];
	    block_iface_j[i][j].tangent[1] = -1*block_iface_j[i][j].normal[0];

        }
    }
}

//Assign interfaces
//South=0 East=1 North=2 West=3
void SolidBlock::assign_ifaces_to_cells(){
    
    for (int j=0; j<(nnj-1); j++){ 
        for (int i=0; i<(nni-1); i++){
	    
            //Vert ifaces
            block_cells[i][j].iface[3] = &block_iface_i[i][j];
            block_cells[i][j].iface[1] = &block_iface_i[i+1][j];
	    
            //Hoz ifaces
            block_cells[i][j].iface[0] = &block_iface_j[i][j];
            block_cells[i][j].iface[2] = &block_iface_j[i][j+1];
	    //block_iface_j[i][j].print_pos();
        }
    }
}



void SolidBlock::assign_secondary_ifaces(){
    
    std::vector<double> temp_vector(2,0.0),temp_vector1(2,0.0);
    
    //Internal secondary i faces
    //Eastward facing normals
    for (int j = 1; j<nnj-1; ++j) {
        for (int i = 1; i<nni; ++i){
	    //iface position //not sure I care about this?
	    add_vectors(block_secondary_iface_i[i][j].pos,block_cells[i-1][j].pos,block_cells[i-1][j-1].pos);
	    scale_vector(block_secondary_iface_i[i][j].pos, 1/2.0);
	    
	    //iface length
	    subtract_vectors(temp_vector,block_cells[i-1][j].pos,block_cells[i-1][j-1].pos); //Causes negative zero. This may affect error?
	    block_secondary_iface_i[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	    
	    //iface normal
	    temp_vector1[0] = temp_vector[1];
	    temp_vector1[1] = -1*temp_vector[0];
	    scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
	    block_secondary_iface_i[i][j].normal = temp_vector1;
	}
    }
    //Northward facing normals
    for (int j = 1; j< nnj; ++j) {
        for (int i = 1; i<nni-1; ++i){
            //iface position //not sure I care about this?
            add_vectors(block_secondary_iface_j[i][j].pos,block_cells[i-1][j-1].pos,block_cells[i][j-1].pos);
            scale_vector(block_secondary_iface_j[i][j].pos, 1/2.0);
	    
            //iface length
            subtract_vectors(temp_vector,block_cells[i-1][j-1].pos,block_cells[i][j-1].pos); //Causes negative zero. This may affect error?
            block_secondary_iface_j[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	    
            //iface normal
            temp_vector1[0] = temp_vector[1];
            temp_vector1[1] = -1*temp_vector[0];
            scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
            block_secondary_iface_j[i][j].normal = temp_vector1;
        }
    }
    
    //West boundary
    {int i = 0;
        for (int j = 1; j<nnj-1; ++j) {
	    
	    //iface position //not sure I care about this?
	    add_vectors(block_secondary_iface_i[i][j].pos,block_iface_i[i][j].pos,block_iface_i[i][j-1].pos);
	    scale_vector(block_secondary_iface_i[i][j].pos, 1/2.0);
	    
	    //iface length
	    subtract_vectors(temp_vector,block_iface_i[i][j].pos,block_iface_i[i][j-1].pos); //Causes negative zero. This may affect error?
	    block_secondary_iface_i[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	    
	    //iface normal
	    temp_vector1[0] = temp_vector[1];
	    temp_vector1[1] = -1*temp_vector[0];
	    scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
	    block_secondary_iface_i[i][j].normal = temp_vector1;
        }
	//Northward facing normals
        for (int j = 1; j< nnj; ++j) {
	    
            //iface position //not sure I care about this?
            add_vectors(block_secondary_iface_j[i][j].pos,block_iface_i[i][j-1].pos,block_cells[i][j-1].pos);
            scale_vector(block_secondary_iface_j[i][j].pos, 1/2.0);
	    
            //iface length
            subtract_vectors(temp_vector,block_iface_i[i][j-1].pos,block_cells[i][j-1].pos); //Causes negative zero. This may affect error?
            block_secondary_iface_j[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	    
            //iface normal
            temp_vector1[0] = temp_vector[1];
            temp_vector1[1] = -1*temp_vector[0];
            scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
            block_secondary_iface_j[i][j].normal = temp_vector1;
        }
    }
    
    //EAST boundary
    {int i = nni;
        for (int j = 1; j<nnj-1; ++j) {
	    
	    //iface position //not sure I care about this?
	    add_vectors(block_secondary_iface_i[i][j].pos,block_iface_i[i-1][j].pos,block_iface_i[i-1][j-1].pos);
	    scale_vector(block_secondary_iface_i[i][j].pos, 1/2.0);
	    
	    //iface length
	    subtract_vectors(temp_vector,block_iface_i[i-1][j].pos,block_iface_i[i-1][j-1].pos); //Causes negative zero. This may affect error?
	    block_secondary_iface_i[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	    
	    //iface normal
	    temp_vector1[0] = temp_vector[1];
	    temp_vector1[1] = -1*temp_vector[0];
	    scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
	    block_secondary_iface_i[i][j].normal = temp_vector1;
        }
	//Northward facing normals
        for (int j = 1; j< nnj; ++j) {
	    
            //iface position //not sure I care about this?
            add_vectors(block_secondary_iface_j[i-1][j].pos,block_cells[i-2][j-1].pos,block_iface_i[i-1][j-1].pos);
            scale_vector(block_secondary_iface_j[i-1][j].pos, 1/2.0);
	    
            //iface length
            subtract_vectors(temp_vector,block_cells[i-2][j-1].pos,block_iface_i[i-1][j-1].pos); //Causes negative zero. This may affect error?
            block_secondary_iface_j[i-1][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	    
            //iface normal
            temp_vector1[0] = temp_vector[1];
            temp_vector1[1] = -1*temp_vector[0];
            scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
            block_secondary_iface_j[i-1][j].normal = temp_vector1;
        }
    }
    
    //South Boundary
    //Eastward facing normals
    {int j=0;
        for (int i = 1; i<nni; ++i){
	    //iface position //not sure I care about this?
	    add_vectors(block_secondary_iface_i[i][j].pos,block_cells[i-1][j].pos,block_iface_j[i-1][j].pos);
	    scale_vector(block_secondary_iface_i[i][j].pos, 1/2.0);
	    
	    //iface length
	    subtract_vectors(temp_vector,block_cells[i-1][j].pos,block_iface_j[i-1][j].pos); //Causes negative zero. This may affect error?
	    block_secondary_iface_i[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	    
	    //iface normal
	    temp_vector1[0] = temp_vector[1];
	    temp_vector1[1] = -1*temp_vector[0];
	    scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
	    block_secondary_iface_i[i][j].normal = temp_vector1;
	}
	//Northward facing normals
	
        for (int i = 1; i<nni-1; ++i){
            //iface position //not sure I care about this?
            add_vectors(block_secondary_iface_j[i][j].pos,block_iface_j[i-1][j].pos,block_iface_j[i][j].pos);
            scale_vector(block_secondary_iface_j[i][j].pos, 1/2.0);
	    
            //iface length
            subtract_vectors(temp_vector,block_iface_j[i-1][j].pos,block_iface_j[i][j].pos); //Causes negative zero. This may affect error?
            block_secondary_iface_j[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	    
            //iface normal
            temp_vector1[0] = temp_vector[1];
            temp_vector1[1] = -1*temp_vector[0];
            scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
            block_secondary_iface_j[i][j].normal = temp_vector1;
        }
    }
    
    //North Boundary
    //Eastward facing normals
    {int j=nnj;
        for (int i = 1; i<nni; ++i){
	    //iface position //not sure I care about this?
	    add_vectors(block_secondary_iface_i[i][j-1].pos,block_iface_j[i-1][j-1].pos,block_cells[i-1][j-2].pos);
	    scale_vector(block_secondary_iface_i[i][j-1].pos, 1/2.0);
	    
	    //iface length
	    subtract_vectors(temp_vector,block_iface_j[i-1][j-1].pos,block_cells[i-1][j-2].pos); //Causes negative zero. This may affect error?
	    block_secondary_iface_i[i][j-1].length = norm2_vector(temp_vector,0,temp_vector.size());
	    
	    //iface normal
	    temp_vector1[0] = temp_vector[1];
	    temp_vector1[1] = -1*temp_vector[0];
	    scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
	    block_secondary_iface_i[i][j-1].normal = temp_vector1;
	}
	//Northward facing normals
	
        for (int i = 1; i<nni-1; ++i){
            //iface position //not sure I care about this?
            add_vectors(block_secondary_iface_j[i][j].pos,block_iface_j[i-1][j-1].pos,block_iface_j[i][j-1].pos);
            scale_vector(block_secondary_iface_j[i][j].pos, 1/2.0);
	    
            //iface length
            subtract_vectors(temp_vector,block_iface_j[i-1][j-1].pos,block_iface_j[i][j-1].pos); //Causes negative zero. This may affect error?
            block_secondary_iface_j[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	    
            //iface normal
            temp_vector1[0] = temp_vector[1];
            temp_vector1[1] = -1*temp_vector[0];
            scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
            block_secondary_iface_j[i][j].normal = temp_vector1;
        }
    }
    
    //SW Corner
    {int i = 0; int j=0;
	//Eastward facing normals
        //iface position //not sure I care about this?
        add_vectors(block_secondary_iface_i[i][j].pos,block_iface_i[i][j].pos,block_vertices[i][j].pos);
        scale_vector(block_secondary_iface_i[i][j].pos, 1/2.0);
	
        //iface length
        subtract_vectors(temp_vector,block_iface_i[i][j].pos,block_vertices[i][j].pos); //Causes negative zero. This may affect error?
        block_secondary_iface_i[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	
        //iface normal
        temp_vector1[0] = temp_vector[1];
        temp_vector1[1] = -1*temp_vector[0];
        scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
        block_secondary_iface_i[i][j].normal = temp_vector1;
	
	//Northward facing normals
        //iface position //not sure I care about this?
        add_vectors(block_secondary_iface_j[i][j].pos,block_vertices[i][j].pos,block_iface_j[i][j].pos);
        scale_vector(block_secondary_iface_j[i][j].pos, 1/2.0);
	
        //iface length
        subtract_vectors(temp_vector,block_vertices[i][j].pos,block_iface_j[i][j].pos); //Causes negative zero. This may affect error?
        block_secondary_iface_j[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());

        //iface normal
        temp_vector1[0] = temp_vector[1];
        temp_vector1[1] = -1*temp_vector[0];
        scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
        block_secondary_iface_j[i][j].normal = temp_vector1;
    }

//SE Corner
    {int i = nni; int j=0;
	//Eastward facing normals
        //iface position //not sure I care about this?
        add_vectors(block_secondary_iface_i[i][j].pos,block_iface_i[i-1][j].pos,block_vertices[i-1][j].pos);
        scale_vector(block_secondary_iface_i[i][j].pos, 1/2.0);
	
        //iface length
        subtract_vectors(temp_vector,block_iface_i[i-1][j].pos,block_vertices[i-1][j].pos); //Causes negative zero. This may affect error?
        block_secondary_iface_i[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	
        //iface normal
        temp_vector1[0] = temp_vector[1];
        temp_vector1[1] = -1*temp_vector[0];
        scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
        block_secondary_iface_i[i][j].normal = temp_vector1;
	
	//Northward facing normals
        //iface position //not sure I care about this?
        add_vectors(block_secondary_iface_j[i-1][j].pos,block_iface_j[i-2][j].pos,block_vertices[i-1][j].pos);
        scale_vector(block_secondary_iface_j[i-1][j].pos, 1/2.0);
	
        //iface length
        subtract_vectors(temp_vector,block_iface_j[i-2][j].pos,block_vertices[i-1][j].pos); //Causes negative zero. This may affect error?
        block_secondary_iface_j[i-1][j].length = norm2_vector(temp_vector,0,temp_vector.size());
	
        //iface normal
        temp_vector1[0] = temp_vector[1];
        temp_vector1[1] = -1*temp_vector[0];
        scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
        block_secondary_iface_j[i-1][j].normal = temp_vector1;
    }
    
    //NE Corner
    {int i = nni; int j=nnj;
	//Eastward facing normals
        //iface position //not sure I care about this?
        add_vectors(block_secondary_iface_i[i][j-1].pos,block_vertices[i-1][j-1].pos,block_iface_i[i-1][j-2].pos);
        scale_vector(block_secondary_iface_i[i][j-1].pos, 1/2.0);
	
        //iface length
        subtract_vectors(temp_vector,block_vertices[i-1][j-1].pos,block_iface_i[i-1][j-2].pos); //Causes negative zero. This may affect error?
        block_secondary_iface_i[i][j-1].length = norm2_vector(temp_vector,0,temp_vector.size());

        //iface normal
        temp_vector1[0] = temp_vector[1];
        temp_vector1[1] = -1*temp_vector[0];
        scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
        block_secondary_iface_i[i][j-1].normal = temp_vector1;

	//Northward facing normals
        //iface position //not sure I care about this?
        add_vectors(block_secondary_iface_j[i-1][j].pos,block_iface_j[i-2][j-1].pos,block_vertices[i-1][j-1].pos);
        scale_vector(block_secondary_iface_j[i-1][j].pos, 1/2.0);

        //iface length
        subtract_vectors(temp_vector,block_iface_j[i-2][j-1].pos,block_vertices[i-1][j-1].pos); //Causes negative zero. This may affect error?
        block_secondary_iface_j[i-1][j].length = norm2_vector(temp_vector,0,temp_vector.size());

        //iface normal
        temp_vector1[0] = temp_vector[1];
        temp_vector1[1] = -1*temp_vector[0];
        scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
        block_secondary_iface_j[i-1][j].normal = temp_vector1;
    }
    
    //NW Corner
    {int i = 0; int j=nnj;
	//Eastward facing normals
        //iface position //not sure I care about this?
        add_vectors(block_secondary_iface_i[i][j-1].pos,block_vertices[i][j-1].pos,block_iface_i[i][j-2].pos);
        scale_vector(block_secondary_iface_i[i][j-1].pos, 1/2.0);
	
        //iface length
        subtract_vectors(temp_vector,block_vertices[i][j-1].pos,block_iface_i[i][j-2].pos); //Causes negative zero. This may affect error?
        block_secondary_iface_i[i][j-1].length = norm2_vector(temp_vector,0,temp_vector.size());

        //iface normal
        temp_vector1[0] = temp_vector[1];
        temp_vector1[1] = -1*temp_vector[0];
        scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
        block_secondary_iface_i[i][j-1].normal = temp_vector1;

	//Northward facing normals
        //iface position //not sure I care about this?
        add_vectors(block_secondary_iface_j[i][j].pos,block_vertices[i][j-1].pos,block_iface_j[i][j-1].pos);
        scale_vector(block_secondary_iface_j[i][j].pos, 1/2.0);

        //iface length
        subtract_vectors(temp_vector,block_vertices[i][j-1].pos,block_iface_j[i][j-1].pos); //Causes negative zero. This may affect error?
        block_secondary_iface_j[i][j].length = norm2_vector(temp_vector,0,temp_vector.size());

        //iface normal
        temp_vector1[0] = temp_vector[1];
        temp_vector1[1] = -1*temp_vector[0];
        scale_vector(temp_vector1,1/norm2_vector(temp_vector1,0,temp_vector1.size()));
        block_secondary_iface_j[i][j].normal = temp_vector1;
    }
    
}

//Using primary cell centres as vertices for secondary cell
void SolidBlock::assign_internal_secondary_geometry_to_vertices(){ //maybe just write an update function so don have to waste memory etc passing info
    
    std::vector<double> v0(2,0.0),v1(2,0.0),v2(2,0.0),v3(2,0.0),p1(2,0.0),p2(2,0.0);
    
    for (int j = 0; j<nnj; ++j) {
        for (int i = 0; i<nni; ++i){
	    
            //Interfaces
            //Hoz (north bdy)
            block_vertices[i][j].iface[0] = &block_secondary_iface_j[i][j];
            block_vertices[i][j].iface[2] = &block_secondary_iface_j[i][j+1];
            //Vert (east bdy)
            block_vertices[i][j].iface[1] = &block_secondary_iface_i[i+1][j];
            block_vertices[i][j].iface[3] = &block_secondary_iface_i[i][j];
	    
            //Area this might be slow or cause round off errors
            //vectors parrallel to interfaces should be +ve i direction
            p1[0] = block_vertices[i][j].iface[0] -> normal[1];
            p1[1] = -1.0 * block_vertices[i][j].iface[0] -> normal[0];
            p2[0] = block_vertices[i][j].iface[2] -> normal[1];
            p2[1] = -1.0 * block_vertices[i][j].iface[2] -> normal[0];
            //cout << p1[0] << p1[1] <<endl;
            scale_vector(p1, 0.5 * block_vertices[i][j].iface[0]->length);
            scale_vector(p2, 0.5 * block_vertices[i][j].iface[0]->length);
            //vertices defining secondary cells
            subtract_vectors(v0, block_vertices[i][j].iface[0]->pos, p1);
            add_vectors(v1, block_vertices[i][j].iface[0]->pos, p1);
            subtract_vectors(v3, block_vertices[i][j].iface[2]->pos, p2);
            add_vectors(v2, block_vertices[i][j].iface[2]->pos, p2);
            block_vertices[i][j].volume = 0.5*abs( (v0[0]-v2[0]) * (v1[1]-v3[1]) - (v1[0]-v3[0]) * (v0[1]-v2[1]) );
        }
    }
}

//Thermal properties
// TODO: Non-constant, non-homogenous and isotropic
// TODO: in update stage. currently using block properties. should use cell properties but still valid for homgenity. Note this will suck up more memory.
void SolidBlock::assign_properties_to_cells(){
    if (anisotropic_flag == 0){
	for (int j=0; j<(nnj-1); j++){
	    for (int i=0; i<(nni-1); i++){
		block_cells[i][j].rho = rho;
		block_cells[i][j].cp = cp;
		block_cells[i][j].k = k;
	    }
	}
    }
    else if (anisotropic_flag == 1){
	for (int j=0; j<(nnj-1); j++){
	    for (int i=0; i<(nni-1); i++){
		block_cells[i][j].rho = rho;
		block_cells[i][j].cp = cp;
		block_cells[i][j].k11 = k11;
		block_cells[i][j].k12 = k12;
		block_cells[i][j].k22 = k22;
	    }
	}
    }
    else {
	std::cout << "Something went wrong assigning properties to cells.\n" \
		  << "Bailing out!" << std::endl;
	exit(1);
    }
}

void SolidBlock::initialise_cells(double T_init){ //probably already done in creating objects but can now set specifics
    //Cell temp
    for (int j=0; j<(nnj-1); j++){
	
	
        for (int i=0; i<(nni-1); i++){
	    
            //T_init = 100 + 100*i;//To test gradient implementation
	    
            block_cells[i][j].T = T_init;
            block_cells[i][j].e = block_cells[i][j].rho * block_cells[i][j].cp * block_cells[i][j].T;
        }
    }
}
