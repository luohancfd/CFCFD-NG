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

/// \brief Calculate shock speed at interface for 2D. 
/// See Ian Johnston's thesis for an explanation.
///
int Block::calc_boundary_vertex_velocity(FV_Interface &IFace1, FV_Interface &IFace2,     
                                         FV_Vertex &vtx, Vector3 trv, size_t gtl )
{   
    double w1, w2;
    Vector3 ws1, ws2, vp;
    vp = vtx.pos[gtl];
    velocity_weighting_factor(IFace1, vp, w1, ws1);
    velocity_weighting_factor(IFace2, vp, w2, ws2);
    if ( (w1 + w2) < 1.0e-3 ) {
	w1 = w2 = 1.0;
    }
    Vector3 wv = (w1*ws1 + w2*ws2) / (w1 + w2);
    // Finally, constrain vertex velocity to body radial direction.
    vtx.vel[gtl] = dot(wv, trv) * trv; 
    return SUCCESS;				       			           
}

/// \brief Calculate shock speed at interface for 3D. 
/// See Ian Johnston's thesis for an explanation.
/// 
int Block::calc_boundary_vertex_velocity(FV_Interface &IFace1, FV_Interface &IFace2,     
                                         FV_Interface &IFace3, FV_Interface &IFace4,
                                         FV_Vertex &vtx, Vector3 trv, size_t gtl)
{   
    double w1, w2, w3, w4;
    Vector3 ws1, ws2, ws3, ws4, vp;
    vp = vtx.pos[gtl];
    velocity_weighting_factor(IFace1, vp, w1, ws1);
    velocity_weighting_factor(IFace2, vp, w2, ws2);
    velocity_weighting_factor(IFace3, vp, w3, ws3);
    velocity_weighting_factor(IFace4, vp, w4, ws4);
    if ( (w1 + w2 + w3 + w4) < 1.0e-3 ) {
	w1 = w2 = w3 = w4 = 1.0;
    }
    Vector3 wv = (w1*ws1 + w2*ws2 + w3*ws3 + w4*ws4) / (w1 + w2 + w3 + w4);
    // Finally, constrain vertex velocity to body radial direction.
    vtx.vel[gtl] = dot(wv, trv) * trv; 
    return SUCCESS;					       			           
}

int Block::velocity_weighting_factor(FV_Interface &IFace, Vector3 vp, double &w, Vector3 &ws)
{
    double M = dot(IFace.fs->vel, unit(vp - IFace.pos)) / IFace.fs->gas->a;
    if ( vabs(IFace.vel) == vabs(IFace.vel) ) // If not NaN.
	ws = IFace.vel;
    else {
	ws = 0.0;
    }
    if ( M > 1 ) {
        w = M;
    } else if ( M != M ) { // Check for NaN.
	w = 0.0;
    } else {
        w = 0.125*( pow(M+1, 2) + (M+1)*fabs(M+1) );
    }
    return SUCCESS;
}

int Block::set_geometry_velocities(size_t dimensions, size_t gtl)
{
    if ( dimensions == 2 ) { 
	set_vertex_velocities2D(gtl);
	set_interface_velocities2D(gtl);
    } else {
	set_vertex_velocities3D(gtl);
	set_interface_velocities3D(gtl);
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
    // Only works with one block in the i-direction. Supports multiple blocks
    // in the j-direction.
    size_t i, j, k;
    FV_Interface *IFaceU, *IFaceD;
    FV_Vertex *svtx, *wvtx, *vtx;
    Vector3 trv;
    double length;

    i = imin;
    k = kmin;
    // Set boundary vertex velocities.
    // Ghost cell geometry will be invalid, but NaNs will be caught by the weighting function.
    for (j = jmin; j <= jmax+1; ++j) {
	IFaceD = get_ifi(i,j-1,k);
	IFaceU = get_ifi(i,j,k);
	vtx = get_vtx(i,j,k);
	wvtx = get_vtx(imax,j,k);
	// Direction vector from vertex to body.
	trv = unit(wvtx->pos[gtl] - vtx->pos[gtl]); 
	calc_boundary_vertex_velocity(*IFaceD, *IFaceU, *vtx, trv, gtl);
    } // for j
    // // Set first and last two boundary vertex velocities
    // if ( bcp[SOUTH]->type_code != ADJACENT ) { // If not adjacent to another block on south side.
    // 	// First
    // 	vtx = get_vtx(imin,jmin,k);
    // 	wvtx = get_vtx(imax,jmin,k);
    // 	IFaceU = get_ifi(imin,jmin,k);
    // 	trv = unit(wvtx->pos[gtl] - vtx->pos[gtl]); 	    
    // 	vtx->vel[gtl] = dot(IFaceU->vel, trv) * trv; // Constrain vertex velocity to body radial direction.
    // }
    // if ( bcp[NORTH]->type_code != ADJACENT ) { // If not adjacent to another block on north side.
    // 	// Last
    // 	vtx = get_vtx(imin,jmax+1,k);
    // 	wvtx = get_vtx(imax,jmax+1,k);
    // 	IFaceD = get_ifi(imin,jmax,k);
    // 	trv = unit(wvtx->pos[gtl] - vtx->pos[gtl]);
    // 	vtx->vel[gtl] = dot(IFaceD->vel, trv) * trv; // Constrain vertex velocity to body radial direction.
    // }
    // Set interior vertex velocities.
    // Velocities are set as linear functions of position between
    // the shock boundary and the wall.
    for (j = jmin; j <= jmax+1; ++j) {
	svtx = get_vtx(imin,j,k); // Shock boundary vertex.
	wvtx = get_vtx(imax+1,j,k); // Wall vertex.
	length = vabs(svtx->pos[gtl] - wvtx->pos[gtl]);
	for (i = imin; i <= imax; ++i) {
	    vtx = get_vtx(i,j,k);
	    vtx->vel[gtl] = (vabs(vtx->pos[gtl] - wvtx->pos[gtl])/length) * svtx->vel[gtl];
	}
    }
    return SUCCESS;
}

/// \brief Set vertex velocities based on previously calculated boundary interface velocities in 3D.
///
///  Based on Ian Johnston's thesis, see for explanation.
///  Assumes inflow at west boundary and wall at east boundary.
///
int Block::set_vertex_velocities3D(size_t gtl)
{
    // Only works with one block in the i-direction. Supports multiple blocks
    // in the j-direction.
    size_t i, j, k;
    FV_Interface *IFace1, *IFace2, *IFace3, *IFace4;
    FV_Vertex *svtx, *wvtx, *vtx;
    Vector3 trv;
    double length;
    i = imin;
    // Set boundary vertex velocities.
    // Ghost cell geometry will be invalid, but NaNs will be caught by the weighting function.
    for (k = kmin; k <= kmax+1; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    IFace1 = get_ifi(i,j,k);
	    IFace2 = get_ifi(i,j-1,k);
	    IFace3 = get_ifi(i,j,k-1);
	    IFace4 = get_ifi(i,j-1,k-1);
	    vtx = get_vtx(i,j,k);
	    wvtx = get_vtx(imax,j,k);
	    // Direction vector from vertex to body.
	    trv = unit(wvtx->pos[gtl] - vtx->pos[gtl]); 
	    calc_boundary_vertex_velocity(*IFace1, *IFace2, *IFace3, *IFace4, *vtx, trv, gtl);
	} // for j
    } // for k
   // // Set vertex velocities on edges of domain as ghost cell will not have valid data.
   //  if ( bcp[SOUTH]->type_code != ADJACENT ) { // If not adjacent to another block on south side.
   // 	for (k = kmin; k <= kmax+1; ++k) {
   // 	    vtx = get_vtx(imin,jmin,kmin);
   // 	    wvtx = get_vtx(imax,jmin,kmin);
   // 	    IFace1 = get_ifi(imin,jmin,kmin);
   // 	    IFace1 = get_ifi(imin,jmin,kmin);
   // 	    trv = unit(wvtx->pos[gtl] - vtx->pos[gtl]); 	    
   // 	    calc_boundary_vertex_velocity(*IFace1, *IFace2, *vtx, trv, gtl);
   // 	}
   //  }
   //  if ( bcp[NORTH]->type_code != ADJACENT ) { // If not adjacent to another block on north side.
   // 	for (k = kmin; k <= kmax+1; ++k) {
   // 	    vtx = get_vtx(imin,jmax+1,kmin);
   // 	    wvtx = get_vtx(imax,jmax+1,kmin);
   // 	    IFace1 = get_ifi(imin,jmax,kmin);
   // 	    trv = unit(wvtx->pos[gtl] - vtx->pos[gtl]);
   // 	    calc_boundary_vertex_velocity(*IFace1, *IFace2, *vtx, trv, gtl);
   // 	}
   //  }
    // Set interior vertex velocities.
    // Velocities are set as linear functions of position between
    // the shock boundary and the wall.
    for (k = kmin; k <= kmax+1; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
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
}


/// \brief Function used to test GCL adherence.
///
int Block::set_gcl_test_vertex_velocities2D( size_t gtl )
{
    // Only works with one block in the i-direction. Supports multiple blocks
    // in the j-direction.
    size_t i, j, k;
    FV_Vertex *vtx;
    k = kmin;
    // Set boundary vertex velocities.
    // Ghost cell geometry will be invalid, but NaNs will be caught by the weighting function.
    for (j = jmin; j <= jmax+1; ++j) {
	for (i = imin; i <= imax+1; ++i) {
	    vtx = get_vtx(i,j,k);
	    vtx->vel[gtl] = 1.0 * vtx->pos[gtl];
	}
    } // for j
    return SUCCESS;
}

/// \brief Function used to test GCL adherence. 
///
int Block::set_gcl_test_vertex_velocities3D(size_t gtl)
{
    // Only works with one block in the i-direction. Supports multiple blocks
    // in the j-direction.
    size_t i, j, k;
    FV_Vertex *vtx;
    // Set boundary vertex velocities.
    // Ghost cell geometry will be invalid, but NaNs will be caught by the weighting function.
    for (k = kmin; k <= kmax+1; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
		vtx = get_vtx(i,j,k);
		vtx->vel[gtl] = 1.0 * vtx->pos[gtl];
	    }
	}
    }
    return SUCCESS;
}

/// \brief Set vertex velocities based on previously calculated boundary interface velocities in 2D. 
///
///  Based on Ian Johnston's thesis, see for explanation.
///  Assumes inflow at west boundary and wall at east boundary.
///
int Block::set_gcl_test_random_vertex_velocities2D(size_t gtl)
{
    // Only works with one block in the i-direction. Supports multiple blocks
    // in the j-direction.
    size_t i, j, k;
    FV_Vertex *vtx;
    FV_Cell *cell;
    k = kmin;
    // Set boundary vertex velocities.
    // Ghost cell geometry will be invalid, but NaNs will be caught by the weighting function.
    srand ( 1 );
    for (j = jmin+1; j <= jmax; ++j) {
	for (i = imin+1; i <= imax; ++i) {
	    vtx = get_vtx(i,j,k);
	    cell = get_cell(i,j,k);
	    vtx->vel[gtl].x = static_cast<double>(rand() % static_cast<int>(cell->fs->gas->a / 50.0));
	    vtx->vel[gtl].y = static_cast<double>(rand() % static_cast<int>(cell->fs->gas->a / 50.0));
	    vtx->vel[gtl].z = 0.0;
	    cout << vtx->vel[gtl] << endl;
	    cout << vabs(vtx->vel[gtl]) << endl;
	    cout << cell->fs->gas->a / 50.0 << endl;
	}
    } // for j
    return SUCCESS;
}


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
    Vector3 vpm1, vpm2;
    double xA, xB, yA, yB;
    size_t tl_old = 0;
    k = kmin;
    for (j = jmin; j <= jmax; ++j) {
	for (i = imin; i <= imax+1; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i,j+1,k);
	    IFace = get_ifi(i,j,k);   
	    vpm1 = 0.5 * ( vtx1->pos[tl_old] + vtx1->pos[gtl+1] );
	    vpm2 = 0.5 * ( vtx2->pos[tl_old] + vtx2->pos[gtl+1] );
	    IFace->pos = 0.5 * (vpm1 + vpm2);
	    IFace->vel = 0.5 * (vtx1->pos[gtl+1] + vtx2->pos[gtl+1] - 
	    			vtx1->pos[tl_old] - vtx2->pos[tl_old]) / dt;
            xA = vpm1.x;
	    yA = vpm1.y;
            xB = vpm2.x;
	    yB = vpm2.y;	 
	    // Interface area at midpoint.   
	    IFace->area[gtl] = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA)); 
	    if (get_axisymmetric_flag() == 1) {
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
	    vpm1 = 0.5 * ( vtx1->pos[tl_old] + vtx1->pos[gtl+1] );
	    vpm2 = 0.5 * ( vtx2->pos[tl_old] + vtx2->pos[gtl+1] );
	    IFace->pos = 0.5 * (vpm1 + vpm2);
	    IFace->vel = 0.5 * (vtx1->pos[gtl+1] + vtx2->pos[gtl+1] - 
	    			vtx1->pos[tl_old] - vtx2->pos[tl_old]) / dt;
            xA = vpm2.x;
	    yA = vpm2.y;
            xB = vpm1.x;
	    yB = vpm1.y;
	    // Interface area at midpoint.   
	    IFace->area[gtl] = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA)); 
	    if (get_axisymmetric_flag() == 1) {
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
    size_t tl_old = 0;
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i,j+1,k);
		vtx3 = get_vtx(i,j,k+1);
		vtx4 = get_vtx(i,j+1,k+1);
		IFace = get_ifi(i,j,k);   
		vpm1 = 0.5 * ( vtx1->pos[tl_old] + vtx1->pos[gtl+1] );
		vpm2 = 0.5 * ( vtx2->pos[tl_old] + vtx2->pos[gtl+1] );
		vpm3 = 0.5 * ( vtx3->pos[tl_old] + vtx3->pos[gtl+1] );
		vpm4 = 0.5 * ( vtx4->pos[tl_old] + vtx4->pos[gtl+1] );
		IFace->pos = 0.25 * (vpm1 + vpm2 + vpm3 + vpm4);
		IFace->vel = 0.25 * (vtx1->pos[gtl+1] + vtx2->pos[gtl+1] +
				     vtx3->pos[gtl+1] + vtx4->pos[gtl+1] - 
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
		vpm1 = 0.5 * ( vtx1->pos[tl_old] + vtx1->pos[gtl+1] );
		vpm2 = 0.5 * ( vtx2->pos[tl_old] + vtx2->pos[gtl+1] );
		vpm3 = 0.5 * ( vtx3->pos[tl_old] + vtx3->pos[gtl+1] );
		vpm4 = 0.5 * ( vtx4->pos[tl_old] + vtx4->pos[gtl+1] );
		IFace->pos = 0.25 * (vpm1 + vpm2 + vpm3 + vpm4);
		IFace->vel = 0.25 * (vtx1->pos[gtl+1] + vtx2->pos[gtl+1] +
				     vtx3->pos[gtl+1] + vtx4->pos[gtl+1] - 
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
		vpm1 = 0.5 * ( vtx1->pos[tl_old] + vtx1->pos[gtl+1] );
		vpm2 = 0.5 * ( vtx2->pos[tl_old] + vtx2->pos[gtl+1] );
		vpm3 = 0.5 * ( vtx3->pos[tl_old] + vtx3->pos[gtl+1] );
		vpm4 = 0.5 * ( vtx4->pos[tl_old] + vtx4->pos[gtl+1] );
		IFace->pos = 0.25 * (vpm1 + vpm2 + vpm3 + vpm4);
		IFace->vel = 0.25 * (vtx1->pos[gtl+1] + vtx2->pos[gtl+1] +
				     vtx3->pos[gtl+1] + vtx4->pos[gtl+1] - 
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
	    IFace->vel = (vtx1->vel[gtl] + vtx2->vel[gtl]) / 2.0;
	}
    }
    for (j = jmin; j <= jmax+1; ++j) {
	for (i = imin; i <= imax; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i+1,j,k);
	    IFace = get_ifj(i,j,k);
	    IFace->vel = (vtx1->vel[gtl] + vtx2->vel[gtl]) / 2.0;
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
		IFace->vel = (vtx1->vel[gtl] + vtx2->vel[gtl] + vtx3->vel[gtl] + vtx4->vel[gtl]) / 4.0;
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
		IFace->vel = (vtx1->vel[gtl] + vtx2->vel[gtl] + vtx3->vel[gtl] + vtx4->vel[gtl]) / 4.0;
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
		IFace->vel = (vtx1->vel[gtl] + vtx2->vel[gtl] + vtx3->vel[gtl] + vtx4->vel[gtl]) / 4.0;
	    }
	}
    }
    return SUCCESS;
}

// FIX-ME moving grid : Andrew, if the following functions really are not wanted, can we delete them?

int Block::diffuse_vertex_velocities(double mu, size_t npass, size_t dimensions, size_t gtl)
/// \brief Filter the cell-centred primary variables.
///
/// This filtering is done on a block-by-block basis.
/// Valid flow data are needed in the ghost cells at the edges.
/// \param alpha : filter coefficient (closer to 1.0, more fudging)
/// \param npass : this many passes of the simple averager
//
{
    FV_Vertex *vtx, *vtxN, *vtxS;
    size_t i,j;
    std::vector<Vector3> diffuse(nnj+5);
    i = imin;
    if ( bcp[SOUTH]->type_code != ADJACENT ) { // If not adjacent to another block on south side.
	// First
	get_vtx(i,jmin-1)->vel[gtl] = get_vtx(i,jmin)->vel[gtl];
	// Second
	get_vtx(i,jmin-2)->vel[gtl] = get_vtx(i,jmin+1)->vel[gtl];
    }
    if ( bcp[NORTH]->type_code != ADJACENT ) { // If not adjacent to another block on south side.
	// Last
	get_vtx(i,jmax+2)->vel[gtl] = get_vtx(i,jmax+1)->vel[gtl];
	// Second last
	get_vtx(i,jmax+3)->vel[gtl] = get_vtx(i,jmax)->vel[gtl];
    }
    for (j = jmin; j <= jmax+1; ++j) {
	vtx = get_vtx(i,j);
	vtxN = get_vtx(i,j+1);
	vtxS = get_vtx(i,j-1);
	diffuse[j] = ( 1 - 4 * mu ) * (vtx)->vel[gtl] +
	    mu * ( (vtxN)->vel[gtl] + (vtxS)->vel[gtl]);
    }
     for (j = jmin; j <= jmax+1; ++j) {
	vtx = get_vtx(i,j);
	vtx->vel[gtl] = diffuse[j];
    }
    return SUCCESS;
} // end of diffuse_vertex_velocities()

int Block::anti_diffuse_vertex_velocities(double mu, size_t npass, size_t dimensions, size_t gtl)
/// \brief Filter the cell-centred primary variables.
///
/// This filtering is done on a block-by-block basis.
/// Valid flow data are needed in the ghost cells at the edges.
/// \param alpha : filter coefficient
/// \param npass : this many passes
//
{
    Vector3 m2, m1, p1, p2;
    Vector3 vel;
    std::vector<Vector3> adflux(nnj+5);
    // Apply the anti-diffusion.
    size_t i = imin;
    for ( size_t j = jmin; j <= jmax+2; ++j ) {
	p2 = get_vtx(i,j+1)->vel[gtl];
	p1 = get_vtx(i,j)->vel[gtl];
	m1 = get_vtx(i,j-1)->vel[gtl];
	m2 = get_vtx(i,j-2)->vel[gtl];
	vel.x = calc_anti_diffusive_flux(m2.x, m1.x, p1.x, p2.x, mu);
	vel.y = calc_anti_diffusive_flux(m2.y, m1.y, p1.y, p2.y, mu);
	vel.z = calc_anti_diffusive_flux(m2.z, m1.z, p1.z, p2.z, mu);
	adflux[j] = vel;
    }
    for ( size_t j = jmin; j <= jmax+1; ++j ) {
	FV_Vertex *vtx = get_vtx(i,j);
	vtx->vel[gtl] += adflux[j] - adflux[j+1];
    }
    return SUCCESS;
} // end of anti_diffuse_vertex_velocities()
