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
    const double th = 0.5;
    const double th_inv = 0.5;
    size_t tl_old = 0;
    size_t krangemax = ( dimensions == 2 ) ? kmax : kmax+1;
    for ( size_t k = kmin; k <= krangemax; ++k ) {
	for ( size_t j = jmin; j <= jmax+1; ++j ) {
	    for ( size_t i = imin; i <= imax+1; ++i ) {
		FV_Vertex *vtx = get_vtx(i,j,k);
		vtx->pos[2] = vtx->pos[tl_old] + dt * (th_inv * vtx->vel[0] + th * vtx->vel[1]);
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
	for ( size_t i =0; i < w.size(); ++i ) wv += w[i] * IFaceList[i]->vel;
	wv /= sum_w;
    } else {
	for ( size_t i =0; i < w.size(); ++i ) wv += IFaceList[i]->vel;
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
    global_data &G = *get_global_data_ptr();
    size_t i, j, k;
    FV_Vertex *vtx1, *vtx2;
    FV_Interface *IFace;
    Vector3 vpm1, vpm2;
    double xA, xB, yA, yB;
    size_t tl_old = gtl - 1;
    k = kmin;
    for (j = jmin; j <= jmax; ++j) {
	for (i = imin; i <= imax+1; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i,j+1,k);
	    IFace = get_ifi(i,j,k);   
	    vpm1 = 0.5 * ( vtx1->pos[tl_old] + vtx1->pos[gtl] );
	    vpm2 = 0.5 * ( vtx2->pos[tl_old] + vtx2->pos[gtl] );
	    IFace->pos = 0.5 * (vpm1 + vpm2);
	    IFace->ivel = 0.5 * (vtx1->pos[gtl] + vtx2->pos[gtl] - 
	    			vtx1->pos[tl_old] - vtx2->pos[tl_old]) / dt;		
            xA = vpm1.x;
	    yA = vpm1.y;
            xB = vpm2.x;
	    yB = vpm2.y;	 
	    // Interface area at midpoint.   
	    IFace->area[gtl] = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA)); 
	    if ( G.axisymmetric ) {
		IFace->Ybar = 0.5 * (yA + yB);
                IFace->area[gtl] *= IFace->Ybar;
            }
	}
    }
    for (j = jmin; j <= jmax+1; ++j) {
	for (i = imin; i <= imax; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i+1,j,k);
	    IFace = get_ifj(i,j,k);
	    vpm1 = 0.5 * ( vtx1->pos[tl_old] + vtx1->pos[gtl] );
	    vpm2 = 0.5 * ( vtx2->pos[tl_old] + vtx2->pos[gtl] );
	    IFace->pos = 0.5 * (vpm1 + vpm2);
	    IFace->ivel = 0.5 * (vtx1->pos[gtl] + vtx2->pos[gtl] - 
	    			vtx1->pos[tl_old] - vtx2->pos[tl_old]) / dt;			
            xA = vpm2.x;
	    yA = vpm2.y;
            xB = vpm1.x;
	    yB = vpm1.y;
	    // Interface area at midpoint.   
	    IFace->area[gtl] = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA)); 
	    if ( G.axisymmetric ) {
		IFace->Ybar = 0.5 * (yA + yB);
                IFace->area[gtl] *= IFace->Ybar;
            }
	}
    }
    return SUCCESS;
}

/// \brief Set interface velocities and area to average value over timestep for GCL adherence. 
///
int Block::set_gcl_interface_properties3D(size_t gtl, double dt)
{
    size_t i, j, k;
    FV_Vertex *vtx1, *vtx2, *vtx3, *vtx4;
    FV_Interface *IFace;
    Vector3 vpm1, vpm2, vpm3, vpm4, p1, p2, p3, p4;
    size_t tl_old = gtl - 1;
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i,j+1,k);
		vtx3 = get_vtx(i,j,k+1);
		vtx4 = get_vtx(i,j+1,k+1);
		IFace = get_ifi(i,j,k);   
		vpm1 = 0.5 * ( vtx1->pos[tl_old] + vtx1->pos[gtl] );
		vpm2 = 0.5 * ( vtx2->pos[tl_old] + vtx2->pos[gtl] );
		vpm3 = 0.5 * ( vtx3->pos[tl_old] + vtx3->pos[gtl] );
		vpm4 = 0.5 * ( vtx4->pos[tl_old] + vtx4->pos[gtl] );
		IFace->pos = 0.25 * (vpm1 + vpm2 + vpm3 + vpm4);
		IFace->ivel = 0.25 * (vtx1->pos[gtl] + vtx2->pos[gtl] +
				     vtx3->pos[gtl] + vtx4->pos[gtl] - 
				     vtx1->pos[tl_old] - vtx2->pos[tl_old] - 
				     vtx3->pos[tl_old] - vtx4->pos[tl_old]) / dt;
		p1 = vpm1;
		p4 = vpm2;
		p2 = vpm3;
		p3 = vpm4;
		// Interface area at midpoint.
		IFace->area[gtl] = vabs(0.25 * cross(p2-p1+p3-p4, p4-p1+p3-p2)); 
	    }
	}
    }
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    for (i = imin; i <= imax; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i+1,j,k);
		vtx3 = get_vtx(i,j,k+1);
		vtx4 = get_vtx(i+1,j,k+1);
		IFace = get_ifj(i,j,k);
		vpm1 = 0.5 * ( vtx1->pos[tl_old] + vtx1->pos[gtl] );
		vpm2 = 0.5 * ( vtx2->pos[tl_old] + vtx2->pos[gtl] );
		vpm3 = 0.5 * ( vtx3->pos[tl_old] + vtx3->pos[gtl] );
		vpm4 = 0.5 * ( vtx4->pos[tl_old] + vtx4->pos[gtl] );
		IFace->pos = 0.25 * (vpm1 + vpm2 + vpm3 + vpm4);
		IFace->ivel = 0.25 * (vtx1->pos[gtl] + vtx2->pos[gtl] +
				     vtx3->pos[gtl] + vtx4->pos[gtl] - 
				     vtx1->pos[tl_old] - vtx2->pos[tl_old] - 
				     vtx3->pos[tl_old] - vtx4->pos[tl_old]) / dt;
		p1 = vpm1;
		p4 = vpm2;
		p2 = vpm3;
		p3 = vpm4;
		// Interface area at midpoint.		
		IFace->area[gtl] = vabs(0.25 * cross(p2-p1+p3-p4, p4-p1+p3-p2)); 
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
		vpm1 = 0.5 * ( vtx1->pos[tl_old] + vtx1->pos[gtl] );
		vpm2 = 0.5 * ( vtx2->pos[tl_old] + vtx2->pos[gtl] );
		vpm3 = 0.5 * ( vtx3->pos[tl_old] + vtx3->pos[gtl] );
		vpm4 = 0.5 * ( vtx4->pos[tl_old] + vtx4->pos[gtl] );
		IFace->pos = 0.25 * (vpm1 + vpm2 + vpm3 + vpm4);
		IFace->ivel = 0.25 * (vtx1->pos[gtl] + vtx2->pos[gtl] +
				     vtx3->pos[gtl] + vtx4->pos[gtl] - 
				     vtx1->pos[tl_old] - vtx2->pos[tl_old] - 
				     vtx3->pos[tl_old] - vtx4->pos[tl_old]) / dt;
		p1 = vpm1;
		p4 = vpm2;
		p2 = vpm3;
		p3 = vpm4;
		// Interface area at midpoint.		
		IFace->area[gtl] = vabs(0.25 * cross(p2-p1+p3-p4, p4-p1+p3-p2)); 
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
