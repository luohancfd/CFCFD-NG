/// \file block_moving_grid.cxx
/// \ingroup eilmer3
/// \brief Andrew Pastrello's functions that implement the moving-grid behaviour.
///
/// \version 23-Mar-2013 extracted from block.cxx
///

#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdexcept>
extern "C" {
#include <zlib.h>
}
#include "cell.hh"
#include "kernel.hh"
#include "block.hh"
#include "bc.hh"

//-----------------------------------------------------------------------------


/// \brief Perform Euler step using calculated vertex
/// velocity and set the vertex position for time-level 1.
///
/// First of two stages.
int Block::predict_vertex_positions(size_t dimensions, double dt)
{
    size_t krangemax = ( dimensions == 2 ) ? kmax : kmax+1;
    const double gamma_1 = 1.0;
    for ( size_t k = kmin; k <= krangemax; ++k ) {
	for ( size_t j = jmin; j <= jmax+1; ++j ) {
	    for ( size_t i = imin; i <= imax+1; ++i ) {
		FV_Vertex *vtx = get_vtx(i,j,k);
		vtx->pos[1] = vtx->pos[0] + dt * gamma_1 * vtx->vel[0];
	    }
	}
    }
    return SUCCESS;
}

/// \brief Perform corrector step using re-calculated vertex
/// velocity and set the vertex position for time-level 2.
///
/// Stage 2 of 2.
int Block::correct_vertex_positions(size_t dimensions, double dt)
{
    size_t krangemax = ( dimensions == 2 ) ? kmax : kmax+1;
    for ( size_t k = kmin; k <= krangemax; ++k ) {
	for ( size_t j = jmin; j <= jmax+1; ++j ) {
	    for ( size_t i = imin; i <= imax+1; ++i ) {
		FV_Vertex *vtx = get_vtx(i,j,k);
		vtx->pos[2] = vtx->pos[1];
	    }
	}
    }
    return SUCCESS;
}

double velocity_weighting_factor(FV_Interface &IFace, Vector3 &vp)
// Smooth, upwind-biased weighting function.
// Equation 4.20 in Ian Johnston's thesis.
{
    double M = dot(IFace.fs->vel, unit(vp - IFace.pos)) / IFace.fs->gas->a;
    double w = 0.0;
    if ( M > 1.0 ) {
        w = M;
    } else if ( M > 0.0 ) {
	double Mp1 = M + 1.0;
        w = 0.125*(Mp1*Mp1 + (Mp1)*fabs(Mp1));
    }
    return w;
}

/// \brief Calculate shock speed at interface for 3D. 
/// See Ian Johnston's thesis for an explanation.
int calc_boundary_vertex_velocity(std::vector<FV_Interface *> &IFaceList,
				  FV_Vertex &vtx, Vector3 trv, size_t gtl)
{   
    std::vector<double> w;
    Vector3 wv(0.0,0.0,0.0);
    Vector3 vp = vtx.pos[gtl];
    if ( IFaceList.size() == 0 ) {
	throw std::runtime_error("In calc_boundary_vertex_velocity, No interfaces in list");
    }
    for ( FV_Interface *facep : IFaceList ) {
	w.push_back(velocity_weighting_factor(*facep, vp));
    }
    double sum_w = 0.0; for ( double wi : w ) sum_w += wi;
    if ( sum_w >= 1.0e-3 ) {
	for ( size_t i =0; i < w.size(); ++i ) wv += w[i] * IFaceList[i]->ivel;
	wv /= sum_w;
    } else {
	for ( size_t i =0; i < w.size(); ++i ) wv += IFaceList[i]->ivel;
	wv /= w.size();
    }
    // Finally, constrain vertex velocity to be directed along 
    // the grid-line, toward the body.
    vtx.vel[gtl] = dot(wv, trv) * trv; 
    return SUCCESS;					       			           
}

int Block::set_geometry_velocities(size_t dimensions, size_t gtl)
{
    if ( dimensions == 2 ) { 
	set_vertex_velocities2D(gtl);
	// Probably doesn't need to be here if using GCL
	//set_interface_velocities2D(gtl);
    } else {
	set_vertex_velocities3D(gtl);
	// Probably doesn't need to be here if using GCL
	//set_interface_velocities3D(gtl);
    }
    return SUCCESS;
}

/// \brief Set vertex velocities based on previously calculated boundary interface velocities in 2D. 
///
///  Based on Ian Johnston's thesis, see for explanation.
///  Assumes inflow at west boundary and wall at east boundary.
///
int Block::set_vertex_velocities2D(size_t gtl)
{
    // Only works with one block in the i-direction.
    // Supports multiple blocks in the j-direction.
    std::vector<FV_Interface *> IFaceList;
    FV_Vertex *svtx, *wvtx, *vtx;
    Vector3 trv;
    double length;

    size_t i = imin;
    size_t k = kmin;
    // Set boundary vertex velocities.
    for ( size_t j = jmin; j <= jmax+1; ++j ) {
	IFaceList.clear();
	if ( j > jmin ) IFaceList.push_back(get_ifi(i,j-1,k));
	if ( j < jmax+1 ) IFaceList.push_back(get_ifi(i,j,k));
	vtx = get_vtx(i,j,k);
	wvtx = get_vtx(imax,j,k);
	// Direction vector from vertex to body.
	trv = unit(wvtx->pos[gtl] - vtx->pos[gtl]);
	calc_boundary_vertex_velocity(IFaceList, *vtx, trv, gtl);
    } // for j
    // Set interior vertex velocities.
    // Velocities are set as linear functions of position between
    // the shock boundary and the wall.
    for ( size_t j = jmin; j <= jmax+1; ++j ) {
	svtx = get_vtx(imin,j,k); // Shock boundary vertex.
	wvtx = get_vtx(imax+1,j,k); // Wall vertex.
	length = vabs(svtx->pos[gtl] - wvtx->pos[gtl]);
	for ( size_t i = imin; i <= imax; ++i ) {
	    vtx = get_vtx(i,j,k);
	    vtx->vel[gtl] = (vabs(vtx->pos[gtl] - wvtx->pos[gtl])/length) * svtx->vel[gtl];
	}
    }
    return SUCCESS;
} // end Block::set_vertex_velocities2D()

/// \brief Set vertex velocities based on previously calculated boundary interface velocities in 3D.
///
///  Based on Ian Johnston's thesis, see for explanation.
///  Assumes inflow at west boundary and wall at east boundary.
int Block::set_vertex_velocities3D(size_t gtl)
{
    // Only works with one block in the i-direction. Supports multiple blocks
    // in the j-direction.
    std::vector<FV_Interface *> IFaceList;
    FV_Vertex *svtx, *wvtx, *vtx;
    Vector3 trv;
    double length;
    size_t i = imin;
    // Set boundary vertex velocities.
    // Ghost cell geometry will be invalid, but NaNs will be caught by the weighting function.
    for ( size_t k = kmin; k <= kmax+1; ++k ) {
	for ( size_t j = jmin; j <= jmax+1; ++j ) {
	    IFaceList.clear();
	    IFaceList.push_back(get_ifi(i,j,k));
	    if ( j > jmin && k < kmax+1 ) IFaceList.push_back(get_ifi(i,j-1,k));
	    if ( j < jmax+1 && k < kmax+1 ) IFaceList.push_back(get_ifi(i,j,k));
	    if ( j < jmax+1 && k > kmin ) IFaceList.push_back(get_ifi(i,j,k-1));
	    if ( j > jmin && k > kmin ) IFaceList.push_back(get_ifi(i,j-1,k-1));
	    vtx = get_vtx(i,j,k);
	    wvtx = get_vtx(imax,j,k);
	    // Direction vector from vertex to body.
	    trv = unit(wvtx->pos[gtl] - vtx->pos[gtl]); 
	    calc_boundary_vertex_velocity(IFaceList, *vtx, trv, gtl);
	} // for j
    } // for k
    // Set interior vertex velocities.
    // Velocities are set as linear functions of position between
    // the shock boundary and the wall.
    for ( size_t k = kmin; k <= kmax+1; ++k ) {
	for ( size_t j = jmin; j <= jmax+1; ++j ) {
	    svtx = get_vtx(imin,j,k); // Shock boundary vertex.
	    wvtx = get_vtx(imax+1,j,k); // Wall vertex.
	    length = vabs(svtx->pos[gtl] - wvtx->pos[gtl]);
	    for (i = imin; i <= imax; ++i) {
		vtx = get_vtx(i,j,k);
		vtx->vel[gtl] = (vabs(vtx->pos[gtl] - wvtx->pos[gtl])/length) * svtx->vel[gtl];
	    }
	}
    }
    return SUCCESS;
} // end Block::set_vertex_velocities3D()

int Block::set_gcl_interface_properties(size_t dimensions, size_t gtl, double dt)
{
    if ( dimensions == 2 ) { 
	set_gcl_interface_properties2D(gtl, dt);

    } else {
	set_gcl_interface_properties3D(gtl, dt);
    }
    return SUCCESS;
}

/// \brief Set interface velocities and area to average value over timestep for GCL adherence. 
///
int Block::set_gcl_interface_properties2D( size_t gtl, double dt )
{
    size_t i, j, k;
    FV_Vertex *vtx1, *vtx2;
    FV_Interface *IFace;
    Vector3 pos1, pos2, temp;
    Vector3 averaged_ivel;
    
    k = kmin;
    for (j = jmin; j <= jmax; ++j) {
	for (i = imin; i <= imax+1; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i,j+1,k);
	    IFace = get_ifi(i,j,k);   	
	    pos1 = vtx1->pos[gtl] - vtx2->pos[0];
	    pos2 = vtx2->pos[gtl] - vtx1->pos[0];
	    averaged_ivel = (vtx1->vel[0] + vtx2->vel[0]) / 2.0;
	    // Use effective edge velocity
	    // Reference: D. Ambrosi, L. Gasparini and L. Vigenano
	    // Full Potential and Euler solutios for transonic unsteady flow
	    // Aeronautical Journal November 1994
	    // Eqn 25
	    temp = 0.5 * cross( pos1, pos2 ) / ( dt * IFace->area[gtl] );
	    IFace->ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
	    averaged_ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
	    IFace->ivel.x = temp.z;
	    IFace->ivel.y = averaged_ivel.y;
	    IFace->ivel.z = averaged_ivel.z;
	    averaged_ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);	    
	    IFace->ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);	    				
	}
    }
    for (j = jmin; j <= jmax+1; ++j) {
	for (i = imin; i <= imax; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i+1,j,k);
	    IFace = get_ifj(i,j,k);
	    pos1 = vtx2->pos[gtl] - vtx1->pos[0];
	    pos2 = vtx1->pos[gtl] - vtx2->pos[0];
	    averaged_ivel = (vtx1->vel[0] + vtx2->vel[0]) / 2.0;	    
	    temp = 0.5 * cross( pos1, pos2 ) / ( dt * IFace->area[gtl] );
	    IFace->ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
	    averaged_ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);	    
	    IFace->ivel.x = temp.z;
	    IFace->ivel.y = averaged_ivel.y;
	    IFace->ivel.z = averaged_ivel.z;	    
	    averaged_ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);	    
	    IFace->ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);	  						
	}
    }
    return SUCCESS;
}

/// \brief Set interface velocities and area to average value over timestep for GCL adherence. 
///
int Block::set_gcl_interface_properties3D(size_t gtl, double dt)
{
    size_t i, j, k;
    Vector3 p0, p1, p2, p3, p4, p5, p6, p7;
    Vector3 pmN, pmE, pmS, pmW, pmT, pmB;
    Vector3 centroid;
    double facial_vi;    
    FV_Interface *IFace;
    Vector3 averaged_ivel;
    FV_Vertex *vtx1, *vtx2, *vtx3, *vtx4;
        
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
	    	vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i,j+1,k);
		vtx3 = get_vtx(i,j,k+1);
		vtx4 = get_vtx(i,j+1,k+1);
                averaged_ivel = (vtx1->vel[0] + vtx2->vel[0] + vtx3->vel[0] + vtx4->vel[0]) / 4.0;		
		p0 = get_vtx(i,j,k)->pos[0];
		p1 = get_vtx(i,j+1,k)->pos[0];
		p2 = get_vtx(i,j+1,k+1)->pos[0];
		p3 = get_vtx(i,j,k+1)->pos[0];
		p4 = get_vtx(i,j,k)->pos[gtl];
		p5 = get_vtx(i,j+1,k)->pos[gtl];
		p6 = get_vtx(i,j+1,k+1)->pos[gtl];
		p7 = get_vtx(i,j,k+1)->pos[gtl];		
		centroid = 0.125 * (p0+p1+p2+p3+p4+p5+p6+p7);
                pmN = 0.25*(p3+p2+p6+p7);
                pmE = 0.25*(p1+p2+p6+p5);
                pmS = 0.25*(p0+p1+p5+p4);
                pmW = 0.25*(p0+p3+p7+p4);
                pmT = 0.25*(p4+p5+p6+p7);
                pmB = 0.25*(p0+p1+p2+p3);
                facial_vi = 0.0;
                facial_vi += tetragonal_dipyramid(p6, p7, p3, p2, pmN, centroid);
                facial_vi += tetragonal_dipyramid(p5, p6, p2, p1, pmE, centroid);
                facial_vi += tetragonal_dipyramid(p4, p5, p1, p0, pmS, centroid);
                facial_vi += tetragonal_dipyramid(p7, p4, p0, p3, pmW, centroid);
                facial_vi += tetragonal_dipyramid(p7, p6, p5, p4, pmT, centroid);
                facial_vi += tetragonal_dipyramid(p0, p1, p2, p3, pmB, centroid);               		
		IFace = get_ifi(i,j,k);
                IFace->ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
                averaged_ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
	        IFace->ivel.x = facial_vi / ( dt * IFace->area[gtl] );	
	        IFace->ivel.y = averaged_ivel.y;
	        IFace->ivel.z = averaged_ivel.z;	        
                averaged_ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);	        
	        IFace->ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);
	    }
	}
    }
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    for (i = imin; i <= imax; ++i) {
	    	vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i,j,k+1);
		vtx3 = get_vtx(i+1,j,k);
		vtx4 = get_vtx(i+1,j,k+1);
                averaged_ivel = (vtx1->vel[0] + vtx2->vel[0] + vtx3->vel[0] + vtx4->vel[0]) / 4.0;		
		p0 = get_vtx(i,j,k)->pos[0];
		p1 = get_vtx(i,j,k+1)->pos[0];
		p2 = get_vtx(i+1,j,k+1)->pos[0];
		p3 = get_vtx(i+1,j,k)->pos[0];
		p4 = get_vtx(i,j,k)->pos[gtl];
		p5 = get_vtx(i,j,k+1)->pos[gtl];
		p6 = get_vtx(i+1,j,k+1)->pos[gtl];
		p7 = get_vtx(i+1,j,k)->pos[gtl];
		centroid = 0.125 * (p0+p1+p2+p3+p4+p5+p6+p7);
                Vector3 pmN = 0.25*(p3+p2+p6+p7);
                Vector3 pmE = 0.25*(p1+p2+p6+p5);
                Vector3 pmS = 0.25*(p0+p1+p5+p4);
                Vector3 pmW = 0.25*(p0+p3+p7+p4);
                Vector3 pmT = 0.25*(p4+p5+p6+p7);
                Vector3 pmB = 0.25*(p0+p1+p2+p3);
                facial_vi = 0.0;
                facial_vi += tetragonal_dipyramid(p6, p7, p3, p2, pmN, centroid);
                facial_vi += tetragonal_dipyramid(p5, p6, p2, p1, pmE, centroid);
                facial_vi += tetragonal_dipyramid(p4, p5, p1, p0, pmS, centroid);
                facial_vi += tetragonal_dipyramid(p7, p4, p0, p3, pmW, centroid);
                facial_vi += tetragonal_dipyramid(p7, p6, p5, p4, pmT, centroid);
                facial_vi += tetragonal_dipyramid(p0, p1, p2, p3, pmB, centroid);                   		
		IFace = get_ifj(i,j,k);
                IFace->ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
                averaged_ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
	        IFace->ivel.x = facial_vi / ( dt * IFace->area[gtl] );	
	        IFace->ivel.y = averaged_ivel.y;
	        IFace->ivel.z = averaged_ivel.z;	        
                averaged_ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);	        
	        IFace->ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);
	    }
	}
    }
    for (k = kmin; k <= kmax+1; ++k) {
	for (j = jmin; j <= jmax; ++j) {
	    for (i = imin; i <= imax; ++i) {
	    	vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i+1,j,k);
		vtx3 = get_vtx(i,j+1,k);
		vtx4 = get_vtx(i+1,j+1,k);
                averaged_ivel = (vtx1->vel[0] + vtx2->vel[0] + vtx3->vel[0] + vtx4->vel[0]) / 4.0;		
		p0 = get_vtx(i,j,k)->pos[0];
		p1 = get_vtx(i+1,j,k)->pos[0];
		p2 = get_vtx(i+1,j+1,k)->pos[0];
		p3 = get_vtx(i,j+1,k)->pos[0];
		p4 = get_vtx(i,j,k)->pos[gtl];
		p5 = get_vtx(i+1,j,k)->pos[gtl];
		p6 = get_vtx(i+1,j+1,k)->pos[gtl];
		p7 = get_vtx(i,j+1,k)->pos[gtl];
		centroid = 0.125 * (p0+p1+p2+p3+p4+p5+p6+p7);
                pmN = 0.25*(p3+p2+p6+p7);
                pmE = 0.25*(p1+p2+p6+p5);
                pmS = 0.25*(p0+p1+p5+p4);
                pmW = 0.25*(p0+p3+p7+p4);
                pmT = 0.25*(p4+p5+p6+p7);
                pmB = 0.25*(p0+p1+p2+p3);
                facial_vi = 0.0;
                facial_vi += tetragonal_dipyramid(p6, p7, p3, p2, pmN, centroid);
                facial_vi += tetragonal_dipyramid(p5, p6, p2, p1, pmE, centroid);
                facial_vi += tetragonal_dipyramid(p4, p5, p1, p0, pmS, centroid);
                facial_vi += tetragonal_dipyramid(p7, p4, p0, p3, pmW, centroid);
                facial_vi += tetragonal_dipyramid(p7, p6, p5, p4, pmT, centroid);
                facial_vi += tetragonal_dipyramid(p0, p1, p2, p3, pmB, centroid);              		
		IFace = get_ifk(i,j,k);
                IFace->ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
                averaged_ivel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
	        IFace->ivel.x = facial_vi / ( dt * IFace->area[gtl] );	
	        IFace->ivel.y = averaged_ivel.y;
	        IFace->ivel.z = averaged_ivel.z;	        
                averaged_ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);	        
	        IFace->ivel.transform_to_global(IFace->n, IFace->t1, IFace->t2);
	    }
	}
    }
    return SUCCESS;
}

/// \brief Set interface velocities as average of adjacent vertex velocities. 
///
int Block::set_interface_velocities2D(size_t gtl)
{
    size_t i, j, k;
    FV_Vertex *vtx1, *vtx2;
    FV_Interface *IFace;
    k = kmin;
    for (j = jmin; j <= jmax; ++j) {
	for (i = imin; i <= imax+1; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i,j+1,k);
	    IFace = get_ifi(i,j,k);
	    IFace->ivel = (vtx1->vel[gtl] + vtx2->vel[gtl]) / 2.0;
	}
    }
    for (j = jmin; j <= jmax+1; ++j) {
	for (i = imin; i <= imax; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i+1,j,k);
	    IFace = get_ifj(i,j,k);
	    IFace->ivel = (vtx1->vel[gtl] + vtx2->vel[gtl]) / 2.0;
	}
    }
    return SUCCESS;
}

/// \brief Set interface velocities as average of adjacent vertex velocities. 
///
int Block::set_interface_velocities3D(size_t gtl)
{
    size_t i, j, k;
    FV_Vertex *vtx1, *vtx2, *vtx3, *vtx4;
    FV_Interface *IFace;
    
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i,j+1,k);
		vtx3 = get_vtx(i,j,k+1);
		vtx4 = get_vtx(i,j+1,k+1);
		IFace = get_ifi(i,j,k);
		IFace->ivel = (vtx1->vel[gtl] + vtx2->vel[gtl] + vtx3->vel[gtl] + vtx4->vel[gtl]) / 4.0;
	    }
	}
    }
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    for (i = imin; i <= imax; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i,j,k+1);
		vtx3 = get_vtx(i+1,j,k);
		vtx4 = get_vtx(i+1,j,k+1);
		IFace = get_ifj(i,j,k);
		IFace->ivel = (vtx1->vel[gtl] + vtx2->vel[gtl] + vtx3->vel[gtl] + vtx4->vel[gtl]) / 4.0;
	    }
	}
    }
    for (k = kmin; k <= kmax+1; ++k) {
	for (j = jmin; j <= jmax; ++j) {
	    for (i = imin; i <= imax; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i+1,j,k);
		vtx3 = get_vtx(i,j+1,k);
		vtx4 = get_vtx(i+1,j+1,k);
		IFace = get_ifk(i,j,k);
		IFace->ivel = (vtx1->vel[gtl] + vtx2->vel[gtl] + vtx3->vel[gtl] + vtx4->vel[gtl]) / 4.0;
	    }
	}
    }
    return SUCCESS;
}

int Block::determine_moving_time_step_size()
/// \brief Compute the local moving time step limit for all cells in the block.
///
{
    global_data *gdp = get_global_data_ptr();
    bool first;
    double dt_moving_local, signal;

    first = true;
    for ( FV_Cell *cp: active_cells ) {
	signal = cp->moving_signal_frequency(gdp->dimensions);
	dt_moving_local = gdp->cfl_moving_target / signal; // Recommend a time step size.
	if ( first ) {
	    dt_moving_allow = dt_moving_local;
	    first = false;
	} else {
	    if (dt_moving_local < dt_moving_allow) dt_moving_allow = dt_moving_local;
	}
    } // for cp
    cout << " moving time step for this block is " << dt_moving_allow << endl;
    
    return SUCCESS;
} // end of determine_moving_time_step_size()

// Helper's function
double tetragonal_dipyramid(const Vector3 &p0, const Vector3 &p1, 
				   const Vector3 &p2, const Vector3 &p3, 
				   const Vector3 &pb, const Vector3 &pc)
// Directly copied from geom.cxx by ignoring negative volume
// J. Grandy (1997) Efficient Computation of Volume of Hexahedral Cells UCRL-ID-128886.
// Base of each dipyramid is specified clockwise from the outside.
// pc is apex
// pb is barycentre of base quad.
// base quad p0->p1->p2->p3->p0 counterclockwise when looking from pc toward base.
{
    double volume = dot(pc-pb, cross(p1-p0+p2-p3, p3-p0+p2-p1)) / 12.0;
    return volume;
}

int Block::clear_vertex_velocities(size_t dimensions)
{
    size_t krangemax = ( dimensions == 2 ) ? kmax : kmax+1;
        for ( size_t k = kmin; k <= krangemax; ++k ) {
	    for ( size_t j = jmin; j <= jmax+1; ++j ) {
	        for ( size_t i = imin; i <= imax+1; ++i ) {
	            FV_Vertex *vtx = get_vtx(i,j,k);
	            vtx->vel[0].x = 0.0;
	            vtx->vel[0].y = 0.0;
	            vtx->vel[0].z = 0.0;
	        }
	    }
        }
    return SUCCESS;       
}
