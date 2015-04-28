/**
 * ssolidblock.d
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-22-04
 */

module ssolidblock;

import geom;
import globalconfig;
import solidblock;
import solidfvcell;
import solidfvinterface;
import solidfvvertex;
import solidbc;


class SSolidBlock : SolidBlock {
public:
    size_t nicell;
    size_t njcell;
    size_t nkcell;
    size_t imin, imax;
    size_t jmin, jmax;
    size_t kmin, kmax;
    size_t[] hicell, hjcell, hkcell; // locations of sample cells for history record
    SolidBoundaryCondition[6] bc;

private:
    size_t _nidim;
    size_t _njdim;
    size_t _nkdim;

    SolidFVCell[] _ctr;
    SolidFVInterface[] _ifi;
    SolidFVInterface[] _ifj;
    SolidFVInterface[] _ifk;
    SolidFVVertex[] _vtx;
    //    SolidFVInterface[] _sifi;
    //    SolidFVInterface[] _sifj;
    //    SolidFVInterface[] _sifk;

public:
    // [TODO] Constructors

    // [TODO] Add functionality
    override void assembleArrays() {}
    override void bindFacesAndVerticesToCells() {}

    override void applyPreSpatialDerivAction(double t, int tLevel)
    {
	bc[Face.north].applyPreSpatialDerivAction(t, tLevel);
	bc[Face.east].applyPreSpatialDerivAction(t, tLevel);
	bc[Face.south].applyPreSpatialDerivAction(t, tLevel);
	bc[Face.west].applyPreSpatialDerivAction(t, tLevel);
	if ( GlobalConfig.dimensions == 3 ) {
	    bc[Face.top].applyPreSpatialDerivAction(t, tLevel);
	    bc[Face.bottom].applyPreSpatialDerivAction(t, tLevel);
	}
    }

    override void applyPostFluxAction(double t, int tLevel)
    {
	bc[Face.north].applyPostFluxAction(t, tLevel);
	bc[Face.east].applyPostFluxAction(t, tLevel);
	bc[Face.south].applyPostFluxAction(t, tLevel);
	bc[Face.west].applyPostFluxAction(t, tLevel);
	if ( GlobalConfig.dimensions == 3 ) {
	    bc[Face.top].applyPostFluxAction(t, tLevel);
	    bc[Face.bottom].applyPostFluxAction(t, tLevel);
	}
    }

    override void computeSpatialDerivatives(int ftl)
    {
	// NOTE: This presently uses the Eilmer3 2D formulation for
	// computing spatial derivatives. We might move this to the
	// 3D formulation at some point in the future.
	// [2015-25-04]

	if ( GlobalConfig.dimensions == 3 ) {
	    throw new Error("computeSpatialDerivatives() not implemented for 3D yet.");
	}
	size_t i, j;
	SolidFVVertex vtx;
	SolidFVCell cell;
	SolidFVInterface a, b;
	double xA, xB, xC, xD;
	double yA, yB, yC, yD;
	double TA, TB, TC, TD;
	double areaInv;

	// Work on all internal secondary cells.
	for ( j = jmin+1; j <= jmax; ++j ) {
	    for ( i = imin+1; i <= imax; ++i ) {
		vtx = getVtx(i,j);
		areaInv = 1.0/vtx.area;
		// Corners of secondary cells
		xA = getCell(i,j-1).pos.x;
		yA = getCell(i,j-1).pos.y;
		xB = getCell(i,j).pos.x;
		yB = getCell(i,j).pos.y;
		xC = getCell(i-1,j).pos.x;
		yC = getCell(i-1,j).pos.y;
		xD = getCell(i-1,j-1).pos.x;
		yD = getCell(i-1,j-1).pos.y;
		// Temperature at the corners of the secondary cells
		TA = getCell(i,j-1).T[ftl];
		TB = getCell(i,j).T[ftl];
		TC = getCell(i-1,j).T[ftl];
		TD = getCell(i-1,j-1).T[ftl];
		// Compute derivative using Gauss-Green theorem
		vtx.dTdx = 0.5 * areaInv * 
		    ((TB + TA) * (yB - yA) + (TC + TB) * (yC - yB) +
		     (TD + TC) * (yD - yC) + (TA + TD) * (yA - yD)); 
		vtx.dTdy = -0.5 * areaInv *
		    ((TB + TA) * (xB - xA) + (TC + TB) * (xC - xB) +
		     (TD + TC) * (xD - xC) + (TA + TD) * (xA - xD));
	    } // end for j
	} // end for i

	// Next, work on the edges.
	// EAST boundary.
	i = imax + 1;
	for ( j = jmin+1; j <= jmax; ++j ) {
	    vtx = getVtx(i,j);
	    areaInv = 1.0 / vtx.area;
	    // Corners for the secondary cell
	    xA = getIfi(i, j-1).pos.x;
	    yA = getIfi(i, j-1).pos.y;
	    xB = getIfi(i,j).pos.x;
	    yB = getIfi(i,j).pos.y;
	    xC = getCell(i-1,j).pos.x;
	    yC = getCell(i-1,j).pos.y;
	    xD = getCell(i-1,j-1).pos.x;
	    yD = getCell(i-1,j-1).pos.y;
	    // Temperatures
	    TA = getIfi(i, j-1).T;
	    TB = getIfi(i, j).T;
	    TC = getCell(i-1,j).T[ftl];
	    TD = getCell(i-1,j-1).T[ftl];
	    // Compute derivative using Gauss-Green theorem
	    vtx.dTdx = 0.5 * areaInv * 
		((TB + TA) * (yB - yA) + (TC + TB) * (yC - yB) +
		 (TD + TC) * (yD - yC) + (TA + TD) * (yA - yD)); 
	    vtx.dTdy = -0.5 * areaInv *
		((TB + TA) * (xB - xA) + (TC + TB) * (xC - xB) +
		 (TD + TC) * (xD - xC) + (TA + TD) * (xA - xD));

	}
	// WEST boundary
	i = imin;
	for ( j = jmin+1; j <= jmax; ++j ) {
	    vtx = getVtx(i, j);
	    areaInv = 1.0 / vtx.area;
	    // Corners of the secondary cell
	    xA = getCell(i, j-1).pos.x;
	    yA = getCell(i, j-1).pos.y;
	    xB = getCell(i, j).pos.x;
	    yB = getCell(i, j).pos.y;
	    xC = getIfi(i, j).pos.x;
	    yC = getIfi(i, j).pos.y;
	    xD = getIfi(i, j-1).pos.x;
	    yD = getIfi(i, j-1).pos.y;
	    // Temperatures
	    TA = getCell(i, j-1).T[ftl];
	    TB = getCell(i, j).T[ftl];
	    TC = getIfi(i, j).T;
	    TD = getIfi(i, j-1).T;
	    // Compute derivative using Gauss-Green theorem
	    vtx.dTdx = 0.5 * areaInv * 
		    ((TB + TA) * (yB - yA) + (TC + TB) * (yC - yB) +
		     (TD + TC) * (yD - yC) + (TA + TD) * (yA - yD)); 
	    vtx.dTdy = -0.5 * areaInv *
		((TB + TA) * (xB - xA) + (TC + TB) * (xC - xB) +
		 (TD + TC) * (xD - xC) + (TA + TD) * (xA - xD));

	}
	// NORTH boundary
	j = jmax + 1;
	for ( i = imin+1; i <= imax; ++i ) {
	    vtx = getVtx(i, j);
	    areaInv = 1.0 / vtx.area;
	    // Corners of the secondary cell
	    xA = getCell(i, j-1).pos.x;
	    yA = getCell(i, j-1).pos.y;
	    xB = getIfj(i, j).pos.x;
	    yB = getIfj(i, j).pos.y;
	    xC = getIfj(i-1, j).pos.x;
	    yC = getIfj(i-1, j).pos.y;
	    xD = getCell(i-1, j-1).pos.x;
	    yD = getCell(i-1, j-1).pos.y;
	    // Temperatures
	    TA = getCell(i, j-1).T[ftl];
	    TB = getIfj(i, j).T;
	    TC = getIfj(i-1, j).T;
	    TD = getCell(i-1, j-1).T[ftl];
	    // Compute derivative using Gauss-Green theorem
	    vtx.dTdx = 0.5 * areaInv * 
		((TB + TA) * (yB - yA) + (TC + TB) * (yC - yB) +
		 (TD + TC) * (yD - yC) + (TA + TD) * (yA - yD)); 
	    vtx.dTdy = -0.5 * areaInv *
		((TB + TA) * (xB - xA) + (TC + TB) * (xC - xB) +
		 (TD + TC) * (xD - xC) + (TA + TD) * (xA - xD));
	}
	// SOUTH boundary
	j = jmin;
	for ( i = imin+1; i <= imax; ++i ) {
	    vtx = getVtx(i, j);
	    areaInv = 1.0 / vtx.area;
	    // Corners of the secondary cell
	    xA = getIfj(i, j).pos.x;
	    yA = getIfj(i, j).pos.y;
	    xB = getCell(i, j).pos.x;
	    yB = getCell(i, j).pos.y;
	    xC = getCell(i-1, j).pos.x;
	    yC = getCell(i-1, j).pos.y;
	    xD = getIfj(i-1, j).pos.x;
	    yD = getIfj(i-1, j).pos.y;
	    // Temperatures
	    TA = getIfj(i, j).T;
	    TB = getCell(i, j).T[ftl];
	    TC = getCell(i-1, j).T[ftl];
	    TD = getIfj(i-1, j).T;
	    // Compute derivative using Gauss-Green theorem
	    vtx.dTdx = 0.5 * areaInv * 
		    ((TB + TA) * (yB - yA) + (TC + TB) * (yC - yB) +
		     (TD + TC) * (yD - yC) + (TA + TD) * (yA - yD));
	    vtx.dTdy = -0.5 * areaInv *
		((TB + TA) * (xB - xA) + (TC + TB) * (xC - xB) +
		 (TD + TC) * (xD - xC) + (TA + TD) * (xA - xD));
	}

	// Finally, derivatives at corners
	// NORTH-EAST corner
	i = imax;
	j = jmax;
	vtx = getVtx(i+1, j+1);
	cell = getCell(i, j);
	a = getIfj(i, j+1);
	b = getIfi(i+1, j);
	xA = a.pos.x; yA = a.pos.y;
	xB = b.pos.x; yB = b.pos.y;
	xC = cell.pos.x; yC = cell.pos.y;
	double denom = (xC - xA)*(yB - yA) - (xB - xA)*(yC - yA);
	TA = a.T;
	TB = b.T;
	TC = cell.T[ftl];
	vtx.dTdx = ((TC-TA)*(yB-yA) - (TB-TA)*(yC-yA))/denom;
	vtx.dTdy = ((TB-TA)*(xC-xA) - (TC-TA)*(xB-xA))/denom;
	// SOUTH-EAST corner
	i = imax;
	j = jmin;
	vtx = getVtx(i+1, j);
	cell = getCell(i, j);
	a = getIfj(i, j);
	b = getIfi(i+1, j);
	xA = a.pos.x; yA = a.pos.y;
	xB = b.pos.x; yB = b.pos.y;
	xC = cell.pos.x; yC = cell.pos.y;
	denom = (xC - xA)*(yB - yA) - (xB - xA)*(yC - yA);
	TA = a.T;
	TB = b.T;
	TC = cell.T[ftl];
	vtx.dTdx = ((TC-TA)*(yB-yA) - (TB-TA)*(yC-yA))/denom;
	vtx.dTdy = ((TB-TA)*(xC-xA) - (TC-TA)*(xB-xA))/denom;
	// SOUTH-WEST corner
	i = imin;
	j = jmin;
	vtx = getVtx(i, j);
	cell = getCell(i, j);
	a = getIfj(i, j);
	b = getIfi(i, j);
	xA = a.pos.x; yA = a.pos.y;
	xB = b.pos.x; yB = b.pos.y;
	xC = cell.pos.x; yC = cell.pos.y;
	denom = (xC - xA)*(yB - yA) - (xB - xA)*(yC - yA);
	TA = a.T;
	TB = b.T;
	TC = cell.T[ftl];
	vtx.dTdx = ((TC-TA)*(yB-yA) - (TB-TA)*(yC-yA))/denom;
	vtx.dTdy = ((TB-TA)*(xC-xA) - (TC-TA)*(xB-xA))/denom;
	// NORTH-WEST corner
	i = imin;
	j = jmax;
	vtx = getVtx(i, j+1);
	cell = getCell(i, j);
	a = getIfj(i, j+1);
	b = getIfi(i, j);
	xA = a.pos.x; yA = a.pos.y;
	xB = b.pos.x; yB = b.pos.y;
	xC = cell.pos.x; yC = cell.pos.y;
	denom = (xC - xA)*(yB - yA) - (xB - xA)*(yC - yA);
	TA = a.T;
	TB = b.T;
	TC = cell.T[ftl];
	vtx.dTdx = ((TC-TA)*(yB-yA) - (TB-TA)*(yC-yA))/denom;
	vtx.dTdy = ((TB-TA)*(xC-xA) - (TC-TA)*(xB-xA))/denom;
    }

    override void computeFluxes()
    {
	size_t i, j, k;
	SolidFVInterface IFace;
	SolidFVVertex vtx1, vtx2;
	double dTdx, dTdy;
	double qx, qy;
	double k_eff;
	// East-facing interfaces
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( i = imin; i <= imax + 1; ++i ) {
		IFace = getIfi(i, j);
		vtx1 = getVtx(i, j+1);
		vtx2 = getVtx(i, j);
		dTdx = 0.5*(vtx1.dTdx + vtx2.dTdx);
		dTdy = 0.5*(vtx1.dTdy + vtx2.dTdy);
		qx = -sp.k * dTdx;
		qy = -sp.k * dTdy;
		IFace.flux = qx * IFace.n.x + qy * IFace.n.y;
	    }
	}
	// North-facing interfaces
	for ( j = jmin; j <= jmax + 1; ++j ) {
	    for ( i = imin; i <= imax; ++i ) {
		IFace = getIfj(i, j);
		vtx1 = getVtx(i, j);
		vtx2 = getVtx(i+1, j);
		dTdx = 0.5*(vtx1.dTdx + vtx2.dTdx);
		dTdy = 0.5*(vtx1.dTdy + vtx2.dTdy);
		qx = -sp.k * dTdx;
		qy = -sp.k * dTdy;
		IFace.flux = qx * IFace.n.x + qy * IFace.n.y;
	    }
	}
    }

    // -- Service methods
    @nogc
    size_t toGlobalIndex(size_t i, size_t j, size_t k) const
    in {
	assert(i < _nidim && j < _njdim && k < _nkdim, "Index out of bounds.");
    }
    body {
	return k * (_njdim * _nidim) + j * _nidim + i; 
    }

    size_t[] toIJKIndices(size_t gid) const
    {
	size_t k = gid / (_njdim * _nidim);
	size_t j = (gid - k * (_njdim * _nidim)) / _nidim;
	size_t i = gid - k * (_njdim * _nidim) - j * _nidim;
	return [i, j, k];
    }

    @nogc ref SolidFVCell getCell(size_t i, size_t j, size_t k=0) { return _ctr[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getIfi(size_t i, size_t j, size_t k=0) { return _ifi[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getIfj(size_t i, size_t j, size_t k=0) { return _ifj[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getIfk(size_t i, size_t j, size_t k=0) { return _ifk[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVVertex getVtx(size_t i, size_t j, size_t k=0) { return _vtx[toGlobalIndex(i,j,k)]; }
    /*
    @nogc ref SolidFVInterface getSifi(size_t i, size_t j, size_t k=0) { return _sifi[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getSifj(size_t i, size_t j, size_t k=0) { return _sifj[toGlobalIndex(i,j,k)]; }
    @nogc ref SolidFVInterface getSifk(size_t i, size_t j, size_t k=0) { return _sifk[toGlobalIndex(i,j,k)]; }
    */


}
