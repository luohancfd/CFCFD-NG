/** \file surface.cxx
 *  \ingroup libgeom2
 *  \brief Implementation of the geometric-surface classes. 
 *  \author PJ
 *  \version 09-Jan-2006 initial coding
 *  \version 23-Jan-2006 Added AOPatch
 *  \version 26-Jan-2005 -- Added TrianglePatch class.
 *  \version 02-Mar-2006 -- Extended TrianglePatch class.
 *  \version 2007 -- Melrose (ADFA): MeshPatch class adapted from MeshVolume class.
 *  \version 07-Feb-2008: PolarSurface added for Hannes and Paul's turbomachinery grids.
 *
 * This module is a C++ replacement for the combined C+Python geometry code that
 * was built for mb_cns and Elmer.
 */
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>

#include "geom.hh"
#include "gpath.hh"
#include "../../nm/source/secant.hh"
#include "surface.hh"
#include "../../nm/source/nelmin.hh"
#include "volume.hh"

using namespace std;

//----------------------------------------------------------------------------------

// Base class for geometric surfaces.
ParametricSurface::ParametricSurface( const string label, 
				      double r0, double r1,
				      double s0, double s1 )
    : label(label), r0(r0), r1(r1), s0(s0), s1(s1) {}
ParametricSurface::ParametricSurface( const ParametricSurface &surf )
    : label(surf.label), r0(surf.r0), r1(surf.r1), s0(surf.s0), s1(surf.s1) {}
ParametricSurface::~ParametricSurface() {}
ParametricSurface* ParametricSurface::clone() const
{
    return new ParametricSurface(*this);
}
ParametricSurface* ParametricSurface::copy() const
{
    return new ParametricSurface(*this);
}
Vector3 ParametricSurface::eval( double r, double s ) const 
{ 
    cout << "ParametricSurface::eval() does nothing." << endl;
    return Vector3(0.0, 0.0, 0.0); 
}
Vector3 ParametricSurface::dpdr( double r, double s ) const
{
    // Obtain the derivative approximately, via a finite-difference.
    double dr = 0.001;
    Vector3 p0 = eval(r, s);
    Vector3 derivative = Vector3();
    if ( r+dr > 1.0 ) {
	// r is close to the r=1.0 boundary, use a one-sided difference.
	Vector3 pminus1 = eval(r-dr, s);
	derivative = (p0 - pminus1) / dr;
    } else if ( r-dr < 0.0 ) {
	// r is close to the r=0 boundary, use a one-sided difference.
	Vector3 pplus1 = eval(r+dr, s);
	derivative = (pplus1 - p0) / dr;
    } else {
	// Not near a boundary, use central-difference.
	Vector3 pminus1 = eval(r-dr, s);
	Vector3 pplus1 = eval(r+dr, s);
	derivative = (pplus1 - pminus1) / (2.0 * dr);
    }
    return derivative;
}
Vector3 ParametricSurface::dpds( double r, double s ) const
{
    // Obtain the derivative approximately, via a finite-difference.
    double ds = 0.001;
    Vector3 p0 = eval(r, s);
    Vector3 derivative = Vector3();
    if ( s+ds > 1.0 ) {
	// s is close to the s=1.0 boundary, use a one-sided difference.
	Vector3 pminus1 = eval(r, s-ds);
	derivative = (p0 - pminus1) / ds;
    } else if ( s-ds < 0.0 ) {
	// s is close to the s=0 boundary, use a one-sided difference.
	Vector3 pplus1 = eval(r, s+ds);
	derivative = (pplus1 - p0) / ds;
    } else {
	// Not near a boundary, use central-difference.
	Vector3 pminus1 = eval(r, s-ds);
	Vector3 pplus1 = eval(r, s+ds);
	derivative = (pplus1 - pminus1) / (2.0 * ds);
    }
    return derivative;
}
string ParametricSurface::str() const
{
    ostringstream ost;
    ost << "ParametricSurface(" << ", \"" << label << "\")";
    return ost.str();
}
ParametricSurface* ParametricSurface::translate( const Vector3 &v ) 
{
    cout << "ParametricSurface::translate() does nothing." << endl 
	 << "We should never have called this code." << endl;
    return this;
}
ParametricSurface* ParametricSurface::translate( double vx, double vy, double vz ) 
{
    translate(Vector3(vx, vy, vz));
    return this;
}
ParametricSurface* ParametricSurface::mirror_image( const Vector3 &point, const Vector3 &normal ) 
{
    cout << "ParametricSurface::mirror_image(): does nothing." << endl 
	 << "We should never have called this code." << endl;
    return this;
}
ParametricSurface* ParametricSurface::rotate_about_zaxis( double dtheta ) 
{
    cout << "ParametricSurface::rotate_about_zaxis(): does nothing." << endl 
	 << "We should never have called this code." << endl;
    return this;
}

// Overload stream output for ParametricSurface objects
ostream& operator<<( ostream &os, const ParametricSurface &surf )
{
    os << surf.str();
    return os;
}

//----------------------------------------------------------------------------------

CoonsPatch::CoonsPatch( const Path &_cA, const Path &_cB,
			const Path &_cC, const Path &_cD,
			string label, 
			double r0, double r1, 
			double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), 
      cA(_cA.clone()), cB(_cB.clone()),
      cC(_cC.clone()), cD(_cD.clone())
{
    // We have cloned the boundary paths to ensure their persistence.
    // Compute the corners (for later use in TFI).
    p00 = cA->eval(0.0);
    p10 = cA->eval(1.0);
    p01 = cB->eval(0.0);
    p11 = cB->eval(1.0);
    // Check that the corners close.
    Vector3 p00_alt = cC->eval(0.0);
    Vector3 p10_alt = cD->eval(0.0);
    Vector3 p01_alt = cC->eval(1.0);
    Vector3 p11_alt = cD->eval(1.0);
    double average_length = 0.25 * (vabs(p10-p00) + vabs(p01-p00) +
				    vabs(p11-p10) + vabs(p11-p01));
    double tolerance = 1.0e-6 * (average_length + 1.0);
    if ( vabs(p00 - p00_alt) > tolerance ) {
	cout << "p00 corner appears open: " << p00 << ", " << p00_alt << endl;
    }
    if ( vabs(p10 - p10_alt) > tolerance ) {
	cout << "p10 corner appears open: " << p10 << ", " << p10_alt << endl;
    }
    if ( vabs(p01 - p01_alt) > tolerance ) {
	cout << "p01 corner appears open: " << p01 << ", " << p01_alt << endl;
    }
    if ( vabs(p11 - p11_alt) > tolerance ) {
	cout << "p11 corner appears open: " << p11 << ", " << p11_alt << endl;
    }
}
CoonsPatch::CoonsPatch( const Vector3 p00, const Vector3 p10, 
			const Vector3 p11, const Vector3 p01,
			string label, double r0, double r1,
			double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1),
      cA(0), cB(0), cC(0), cD(0),
      p00(p00), p10(p10), p11(p11), p01(p01)
{
    cA = new Line(p00, p10);
    cB = new Line(p01, p11);
    cC = new Line(p00, p01);
    cD = new Line(p10, p11);
}
CoonsPatch::CoonsPatch( const CoonsPatch &surf )
    : ParametricSurface(surf.label, surf.r0, surf.r1, surf.s0, surf.s1), 
      cA(0), cB(0), cC(0), cD(0), 
      p00(surf.p00), p10(surf.p10), p11(surf.p11), p01(surf.p01) 
{
    cA = surf.cA->clone();
    cB = surf.cB->clone();
    cC = surf.cC->clone();
    cD = surf.cD->clone();
}
CoonsPatch::~CoonsPatch() 
{
    delete cA; delete cB; delete cC; delete cD;
}
CoonsPatch* CoonsPatch::clone() const
{
    return new CoonsPatch(*this);
}
CoonsPatch* CoonsPatch::copy() const
{
    return new CoonsPatch(*this);
}
Vector3 CoonsPatch::eval( double r, double s ) const 
{ 
    // Locate a point on the CoonsPatch surface by blended linear interpolation.
    // @param r: interpolation parameter along curves cA and cB, 0.0<=r<=1.0
    // @param s: interpolation parameter along curves cC and cD, 0.0<=s<=1.0
    // @returns: a L{Vector3} value for the point.
    r = r0 + (r1 - r0) * r;
    s = s0 + (s1 - s0) * s;
    Vector3 cAr = cA->eval(r); 
    Vector3 cBr = cB->eval(r);
    Vector3 cCs = cC->eval(s); 
    Vector3 cDs = cD->eval(s);
    Vector3 p = (1.0-s) * cAr + s * cBr + (1.0-r) * cCs + r * cDs - 
	( (1.0-r)*(1.0-s)*p00 + (1.0-r)*s*p01 + r*(1.0-s)*p10 + r*s*p11 );
    return p; 
}
string CoonsPatch::str() const
{
    ostringstream ost;
    ost << "CoonsPatch("
	<< "[ " << cA->str()
	<< ", " << cB->str()
	<< ", " << cC->str()
	<< ", " << cD->str()
	<< "], \"" << label << ", ";
    ost << r0 << ", " << r1 << ", " << s0 << ", " << s1 << ")";
    return ost.str();
}
CoonsPatch* CoonsPatch::translate( const Vector3 &v ) 
{
    cA->translate(v);
    cB->translate(v);
    cC->translate(v);
    cD->translate(v);
    // Recompute the corners (for later use in TFI).
    p00 = cA->eval(0.0);
    p10 = cA->eval(1.0);
    p01 = cB->eval(0.0);
    p11 = cB->eval(1.0);
    return this;
}
CoonsPatch* CoonsPatch::translate( double vx, double vy, double vz ) 
{
    translate(Vector3(vx, vy, vz));
    return this;
}
CoonsPatch* CoonsPatch::mirror_image( const Vector3 &point, const Vector3 &normal ) 
{
    cA->mirror_image(point, normal);
    cB->mirror_image(point, normal);
    cC->mirror_image(point, normal);
    cD->mirror_image(point, normal);
    // Recompute the corners (for later use in TFI).
    p00 = cA->eval(0.0);
    p10 = cA->eval(1.0);
    p01 = cB->eval(0.0);
    p11 = cB->eval(1.0);
    return this;
}
CoonsPatch* CoonsPatch::rotate_about_zaxis( double dtheta ) 
{
    cA->rotate_about_zaxis(dtheta);
    cB->rotate_about_zaxis(dtheta);
    cC->rotate_about_zaxis(dtheta);
    cD->rotate_about_zaxis(dtheta);
    // Recompute the corners (for later use in TFI).
    p00 = cA->eval(0.0);
    p10 = cA->eval(1.0);
    p01 = cB->eval(0.0);
    p11 = cB->eval(1.0);
    return this;
}

//---------------------------------------------------------------------------------

AOPatch::AOPatch( const Path &cA, const Path &cB,
		  const Path &cC, const Path &cD,
		  string label, int nx, int ny,
		  double r0, double r1, double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), 
      tfi_surface(CoonsPatch(cA, cB, cC, cD)), nx(nx), ny(ny)
{
    allocate_background_mesh();
    compute_background_mesh();
}
AOPatch::AOPatch( const Vector3 p00, const Vector3 p10, 
		  const Vector3 p11, const Vector3 p01,
		  string label, int nx, int ny,
		  double r0, double r1, double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), 
      tfi_surface(CoonsPatch(p00, p10, p11, p01)), nx(nx), ny(ny)
{
    allocate_background_mesh();
    compute_background_mesh();
}
AOPatch::AOPatch( const CoonsPatch &tfi_surf, string label, int nx, int ny,
		  double r0, double r1, double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), 
      tfi_surface(tfi_surf), nx(nx), ny(ny)
{
    allocate_background_mesh();
    compute_background_mesh();
}
AOPatch::AOPatch( const AOPatch &surf )
    : ParametricSurface(surf.label, surf.r0, surf.r1, surf.s0, surf.s1), 
      tfi_surface(surf.tfi_surface), 
      nx(surf.nx), ny(surf.ny) 
{
    allocate_background_mesh();
    compute_background_mesh();
}
AOPatch::~AOPatch() 
{
    free_background_mesh();
}
AOPatch* AOPatch::clone() const
{
    return new AOPatch(*this);
}
AOPatch* AOPatch::copy() const
{
    return new AOPatch(*this);
}
Vector3 AOPatch::eval( double r, double s ) const 
{ 
    r = r0 + (r1 - r0) * r;
    s = s0 + (s1 - s0) * s;

    // Locate a point on the AOPatch surface.
    Vector3 p = tfi_surface.eval(r, s); // TFI is a fall-back

    // Use TFI close to the boundaries.  
    // The background mesh is pretty coarse.
    double tol = 1.0e-4;
    if ( r < tol || s < tol || r > 1.0-tol || s > 1.0-tol ) return p;

    if ( mesh_ok ) {
	// Interpolate within the AO background mesh.
	// This involves finding the relevant coarse cell and
	// using a bilinear interpolation within that cell.
	double dr = 1.0 / nx;
	double ds = 1.0 / ny;
	int ix_coarse = (int) (r / dr);
	int iy_coarse = (int) (s / ds);
	// Parametric coordinate within the coarse cell.
	double local_r = (r - dr * ix_coarse) / dr;
	double local_s = (s - ds * iy_coarse) / ds;
	// BiLinear interpolation for each component.
	p = (1.0 - local_r) * (1.0 - local_s) * bgmesh[ix_coarse][iy_coarse] +
	    (1.0 - local_r) * local_s * bgmesh[ix_coarse][iy_coarse + 1] +
	    local_r * (1.0 - local_s) * bgmesh[ix_coarse + 1][iy_coarse] +
	    local_r * local_s * bgmesh[ix_coarse + 1][iy_coarse + 1];
    }
    return p; 
}
string AOPatch::str() const
{
    ostringstream ost;
    ost << "AOPatch(" << tfi_surface.str() << ", \"" << label << "\"" 
	<< ", " << nx << ", " << ny << ", ";
    ost << r0 << ", " << r1 << ", " << s0 << ", " << s1 << ")";
    return ost.str();
}
AOPatch* AOPatch::translate( const Vector3 &v ) 
{
    tfi_surface.translate(v);
    compute_background_mesh();
    return this;
}
AOPatch* AOPatch::translate( double vx, double vy, double vz ) 
{
    translate(Vector3(vx, vy, vz));
    return this;
}
AOPatch* AOPatch::mirror_image( const Vector3 &point, const Vector3 &normal ) 
{
    tfi_surface.mirror_image(point, normal);
    compute_background_mesh();
    return this;
}
AOPatch* AOPatch::rotate_about_zaxis( double dtheta ) 
{
    tfi_surface.rotate_about_zaxis(dtheta);
    compute_background_mesh();
    return this;
}
void AOPatch::allocate_background_mesh()
{
    mesh_ok = 1;
    bgmesh = (Vector3 **) calloc(nx+1, sizeof(Vector3*));
    if ( bgmesh == NULL ) {
	cout << "AOPatch::allocate_background_mesh() could not allocate memory." << endl;
	mesh_ok = 0;
    }
    for ( int ix = 0; ix <= nx; ++ix ) {
	bgmesh[ix] = (Vector3 *) calloc(ny+1, sizeof(Vector3));
	if ( bgmesh[ix] == NULL ) {
	    cout << "AOPatch::allocate_background_mesh() could not allocate memory." << endl;
	    mesh_ok = 0;
	}
    }
    return;
}
void AOPatch::free_background_mesh()
{
    for ( int ix = 0; ix <= nx; ++ix ) {
	if ( bgmesh[ix] != NULL ) free(bgmesh[ix]);
    }
    if ( bgmesh != NULL ) free(bgmesh);
    mesh_ok = 0;
    return;
}
void AOPatch::compute_background_mesh()
{
    if ( !mesh_ok ) return;  // fall back to TFI surface

    // Initial positions of the mesh points are just TFI locations.
    double dXi = 1.0 / ((double) nx);
    double dEta = 1.0 / ((double) nx);
    double dXi2 = dXi * dXi;
    double dEta2 = dEta * dEta;
    for ( int ix = 0; ix <= nx; ++ix ) {
        for ( int iy = 0; iy <= ny; ++iy ) {
            double r = dXi * ((double) ix);
            double s = dEta * ((double) iy);
            bgmesh[ix][iy] = tfi_surface.eval(r, s);
	}
    }

    // Now, adjust the mesh point locations.
    // Note that the current formulation is for the (x,y)-plane only.
    // z-components remain as the TFI value.
    double x_tol = 1.0e-6;
    double y_tol = 1.0e-6;
    double largest_x_move;
    double largest_y_move;
    int max_count = 500;
    int count;
    for (count = 1; count <= max_count; ++count) {
        largest_x_move = 0.0;
        largest_y_move = 0.0;
        /* Adjust the internal points only. */
        for ( int ix = 1; ix < nx; ++ix ) {
            for ( int iy = 1; iy < ny; ++iy ) {
                /* Save the old position. */
                double x_old = bgmesh[ix][iy].x;
                double y_old = bgmesh[ix][iy].y;

                /* Calculate the partial derivatives. */
                double dxdXi = (bgmesh[ix+1][iy].x - bgmesh[ix-1][iy].x) / (2.0 * dXi);
                double dxdEta = (bgmesh[ix][iy+1].x - bgmesh[ix][iy-1].x) / (2.0 * dEta);
		double d2xdXidEta = ((bgmesh[ix+1][iy+1].x - bgmesh[ix-1][iy+1].x) -
				     (bgmesh[ix+1][iy-1].x - bgmesh[ix-1][iy-1].x)) / 
		    (4.0 * dXi * dEta);
                double dydXi = (bgmesh[ix+1][iy].y - bgmesh[ix-1][iy].y) / (2.0 * dXi);
                double dydEta = (bgmesh[ix][iy+1].y - bgmesh[ix][iy-1].y) / (2.0 * dEta);
                double d2ydXidEta = ((bgmesh[ix+1][iy+1].y - bgmesh[ix-1][iy+1].y) -
				     (bgmesh[ix+1][iy-1].y - bgmesh[ix-1][iy-1].y)) / 
		    (4.0 * dXi * dEta);

                /* Calculate intermediate quantities */
                double B = dxdXi * dydEta + dxdEta * dydXi;
                double Ax = (4.0 * dxdXi * dxdEta * d2xdXidEta +
			     2.0 * B * d2ydXidEta) * (dXi2 * dEta2);
                double Ay = (4.0 * dydXi * dydEta * d2ydXidEta +
			     2.0 * B * d2xdXidEta) * (dXi2 * dEta2);
                double g11 = dxdXi * dxdXi + dydXi * dydXi;
                double g22 = dxdEta * dxdEta + dydEta * dydEta;

                /* Update the node position. */
                double numer = Ax + g22 * dEta2 * (bgmesh[ix+1][iy].x + bgmesh[ix-1][iy].x) +
                    g11 * dXi2 * (bgmesh[ix][iy+1].x + bgmesh[ix][iy-1].x);
                double denom = 2.0 * (g22 * dEta2 + g11 * dXi2);
                bgmesh[ix][iy].x = numer / denom;

                numer = Ay + g22 * dEta2 * (bgmesh[ix+1][iy].y + bgmesh[ix-1][iy].y) +
                    g11 * dXi2 * (bgmesh[ix][iy+1].y + bgmesh[ix][iy-1].y);
                denom = 2.0 * (g22 * dEta2 + g11 * dXi2);
                bgmesh[ix][iy].y = numer / denom;

                double dx = fabs(bgmesh[ix][iy].x - x_old);
                double dy = fabs(bgmesh[ix][iy].y - y_old);
                if ( dx > largest_x_move ) largest_x_move = dx;
                if ( dy > largest_y_move ) largest_y_move = dy;

                // printf("Iteration %d node[%d, %d] x,y(%f, %f)\n",
                //        count, ix, iy, bgmesh[ix][iy].x, bgmesh[ix][iy].y);
            }   /* iy loop */
        }   /* ix loop */

        /* Check for convergence */
        if (largest_x_move <= x_tol && largest_y_move <= y_tol) break;

    }   /* count loop */

    // printf("AO iteration number of iterations = %d.\n", count);
    // printf("Final node movements dx,dy=(%e,%e).\n",
    //        largest_x_move, largest_y_move);
    if ( count > max_count ) {
        printf("AO iteration did not converge.\n");
    }
    return;
}

//----------------------------------------------------------------------------------

MeshPatch::MeshPatch( const vector<Vector3*> &_p, int _ni, int _nj,
		      string label, double r0, double r1, double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1),
      p(vector<Vector3>()), ni(_ni), nj(_nj)
{
    if ( (size_t)(ni*nj) != _p.size() ) {
	cout << "MeshPatch Warning: ni=" << ni << " nj=" << nj 
	     << " p.size()=" << _p.size() << endl;
    }
    for ( size_t i = 0; i < _p.size(); ++i ) p.push_back(*(_p[i]));
}
MeshPatch::MeshPatch( const vector<Vector3> &_p, int _ni, int _nj,
		      string label, double r0, double r1, double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1),
      p(vector<Vector3>()), ni(_ni), nj(_nj)
{
    if ( (size_t)(ni*nj) != _p.size() ) {
	cout << "MeshPatch Warning: ni=" << ni << " nj=" << nj 
	     << " p.size()=" << _p.size() << endl;
    }
    for ( size_t i = 0; i < _p.size(); ++i ) p.push_back(_p[i]);
}
MeshPatch::MeshPatch( const MeshPatch &surf )
    : ParametricSurface(surf.label, surf.r0, surf.r1, surf.s0, surf.s1), 
      p(surf.p), ni(surf.ni), nj(surf.nj) {}
MeshPatch::~MeshPatch() {}
MeshPatch* MeshPatch::clone() const
{
    return new MeshPatch(*this);
} 
MeshPatch* MeshPatch::copy() const
{
    return new MeshPatch(*this);
} 
MeshPatch::MeshPatch( const string file_name, int vtk_format,
			const string label, double r0, double r1, double s0, double s1)
    : ParametricSurface(label, r0, r1, s0, s1),
      p(vector<Vector3>()), ni(0), nj(0)
{
    cout << "MeshPatch constructor: read from file " << file_name << endl;
    FILE* fp = fopen(file_name.c_str(), "r");
    if ( fp == NULL ) {
	cout << "MeshPatch constructor: could not open file " << file_name << endl;
    }

    char buffer[256];
    if ( vtk_format == 1 ) {
	// Discard the first 5 lines of the VTK-format file.
	// line 1 expect  # vtk DataFile Version 2.0
	// line 2 expect  some title string
	// line 3 expect  ASCII
	// line 4 expect  blank line
	// line 5 expect  DATASET STRUCTURED_GRID
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cout << "MeshPatch constructor: problem reading VTK header for grid" << endl;
	}
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cout << "MeshPatch constructor: problem reading VTK header for grid" << endl;
	}
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cout << "MeshPatch constructor: problem reading VTK header for grid" << endl;
	}
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cout << "MeshPatch constructor: problem reading VTK header for grid" << endl;
	}
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cout << "MeshPatch constructor: problem reading VTK header for grid" << endl;
	}
	// Now we get to some interesting data.
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cout << "MeshPatch constructor: problem reading VTK header for grid" << endl;
	}
	sscanf(buffer, "DIMENSIONS %d %d 1", &ni, &nj);
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cout << "MeshPatch constructor: problem reading VTK header for grid" << endl;
	}
	int total_points;
	sscanf(buffer, "POINTS %d float\n", &total_points);
	if ( total_points != ni * nj) {
	    cout << "MeshPatch constructor: wrong number of total points=" 
		 << total_points << " ni*nj=" << ni*nj << endl;
	}
    } else {
	// Assume TECPLOT format.
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cout << "MeshPatch constructor: problem reading TECPLOT header for grid" << endl;
	}
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cout << "MeshPatch constructor: problem reading TECPLOT header for grid" << endl;
	}
	// Now, we get to interesting data.
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cout << "MeshPatch constructor: problem reading TECPLOT header for grid" << endl;
	}
	sscanf(buffer, "ZONE I=%d, J=%d, K=1, F=POINT", &ni, &nj);
    }
    // Read the vertex coordinates.
    // Note the order of the loops: i index runs fastest.
    double x, y, z;
    for ( int j = 0; j < nj; ++j ) {
	for ( int i = 0; i < ni; ++i) {
	    if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
		cout << "MeshPatch constructor: problem reading grid" 
		     << " i=" << i << ", j=" << j << endl;
	    }
	    sscanf(buffer, "%lf %lf %lf", &x, &y, &z);
	    p.push_back( Vector3(x, y, z) );
	}
    }
    fclose(fp);
    cout << "Finished MeshPatch construction." << endl;
}
Vector3 MeshPatch::eval( double r, double s ) const
{
    r = r0 + (r1 - r0) * r;
    s = s0 + (s1 - s0) * s;

    // Interpolate within the background mesh.
    // This involves finding the relevant cell and
    // using a bilinear interpolation within that cell.
    double dr = 1.0 / (ni-1);
    double ds = 1.0 / (nj-1);
    int i = (int) (r / dr); i = (i < 0) ? 0 : i; i = (i > ni-2) ? ni-2 : i;
    int j = (int) (s / ds); j = (j < 0) ? 0 : j; j = (j > nj-2) ? nj-2 : j;
    // Identify the four corners of the cell.
    // This indexing essentially defines the order in which points are
    // stored in the vector p.  The i-index varies most quickly.
    int indx00 = j*ni + i;
    int indx10 = j*ni + i+1;
    int indx11 = (j+1)*ni + i+1;
    int indx01 = (j+1)*ni + i;
    // Parametric coordinate within the coarse cell.
    double local_r = (r - dr * i) / dr;
    double local_s = (s - ds * j) / ds;
    // BiLinear interpolation.
    Vector3 point = (1.0 - local_r) * (1.0 - local_s) * p[indx00] +
	(1.0 - local_r) * local_s * p[indx01] +
	local_r * (1.0 - local_s) * p[indx10] +
	local_r * local_s * p[indx11];
    return point;
}
string MeshPatch::str() const
{
    ostringstream ost;
    size_t i;
    ost << "MeshPatch(";
    ost << "p=[";
    for ( i = 0; i < p.size()-1; ++i ) ost << p[i] << ", ";
    ost << p[p.size()-1] << "], ";
    ost << "ni=" << ni << ", nj=" << nj << ", ";
    ost << "\"" << label << "\"" << ", ";
    ost << r0 << ", " << r1 << ", " << s0 << ", " << s1 << ")";
    return ost.str();
}
MeshPatch* MeshPatch::translate( const Vector3 &v )
{
    for ( size_t i = 0; i < p.size(); ++i ) p[i] += v;
    return this;
}
MeshPatch* MeshPatch::translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    return this;
}
MeshPatch* MeshPatch::mirror_image( const Vector3 &point, const Vector3 &normal )
{
    for ( size_t i = 0; i < p.size(); ++i ) p[i].mirror_image(point, normal);
    return this;
}
MeshPatch* MeshPatch::rotate_about_zaxis( double dtheta )
{
    for ( size_t i = 0; i < p.size(); ++i ) p[i].rotate_about_zaxis(dtheta);
    return this;
}

//---------------------------------------------------------------------------------

TrianglePatch::TrianglePatch( const vector<Vector3*> &_p, const vector<int> &_itri,
			      const vector<int> &_iA, const vector<int> &_iB, 
			      const vector<int> &_iC, const vector<int> &_iD, 
			      string label,
			      double r0, double r1, double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), 
      p(vector<Vector3>()), itri(_itri),
      iA(_iA), iB(_iB), iC(_iC), iD(_iD), tfi_surface(0), ntri(0)
{
    size_t i;
    ntri = itri.size() / 3;
    for ( i = 0; i < _p.size(); ++i ) p.push_back(*(_p[i]));
    // Check user-supplied index values.
    int max_index = (int) p.size() - 1;
    for ( i = 0; i < itri.size(); ++i ) {
	if ( itri[i] > max_index ) {
	    cout << "TrianglePatch constructor: bad index value: " << itri[i] << endl;
	    itri[i] = max_index; // The mesh is buggered anyway.
	}
    }
    for ( i = 0; i < iA.size(); ++i ) {
	if ( iA[i] > max_index ) {
	    cout << "TrianglePatch constructor: bad index value: " << iA[i] << endl;
	    iA[i] = max_index;
	}
    }
    for ( i = 0; i < iB.size(); ++i ) {
	if ( iB[i] > max_index ) {
	    cout << "TrianglePatch constructor: bad index value: " << iB[i] << endl;
	    iB[i] = max_index;
	}
    }
    for ( i = 0; i < iC.size(); ++i ) {
	if ( iC[i] > max_index ) {
	    cout << "TrianglePatch constructor: bad index value: " << iC[i] << endl;
	    iC[i] = max_index;
	}
    }
    for ( i = 0; i < iD.size(); ++i ) {
	if ( iD[i] > max_index ) {
	    cout << "TrianglePatch constructor: bad index value: " << iD[i] << endl;
	    iD[i] = max_index;
	}
    }
    setup_tfi_surface();
}
TrianglePatch::TrianglePatch( const TrianglePatch &surf )
    : ParametricSurface(surf.label, surf.r0, surf.r1, surf.s0, surf.s1), 
      p(surf.p), itri(surf.itri),
      iA(surf.iA), iB(surf.iB), iC(surf.iC), iD(surf.iD), 
      tfi_surface(surf.tfi_surface->clone()), ntri(surf.ntri) {}
TrianglePatch::TrianglePatch( const ParametricSurface *surf, int Nr, int Ns, 
			      string label, double r0, double r1,
			      double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), 
      p(vector<Vector3>()), itri(vector<int>()),
      iA(vector<int>()), iB(vector<int>()), iC(vector<int>()), 
      iD(vector<int>()), tfi_surface(0), ntri(0)
{
    int i, j, i00, i10, i11, i01;
    double r, s;
    double dr = 1.0 / Nr;
    double ds = 1.0 / Ns;
    // Sample the underlying surface to make the vertices of the triangles.
    for ( j = 0; j <= Ns; ++j ) {
	s = j * ds;
	for ( i = 0; i <= Nr; ++i ) {
	    r = i * dr;
	    p.push_back(surf->eval(r,s));
	}
    }
    // Produce the index list for the individual triangles.
    for ( j = 0; j < Ns; ++j ) {
	for ( i = 0; i < Nr; ++i ) {
	    // Each quad panel will be composed of two triangles.
	    //   i01---i11
	    //    |   / |
            //    | /   |
            //   i00---i10
	    i00 = i + j * (Nr+1);
	    i10 = i00 + 1;
	    i01 = i + (j+1) * (Nr+1);
	    i11 = i01 + 1;
	    itri.push_back(i00); itri.push_back(i10); itri.push_back(i11);
	    itri.push_back(i01); itri.push_back(i00); itri.push_back(i11);
	    ntri += 2;
	}
    }
    // Produce index lists for each of the edges.
    for ( i = 0; i <= Nr; ++i ) {
	iA.push_back(i);
	iB.push_back(i + Ns * (Nr+1));
    }
    for ( j = 0; j <= Ns; ++j ) {
	iC.push_back(j * (Nr+1));
	iD.push_back(Nr + j * (Nr+1));
    }
    // Finally, build the query surface (to be used internally in eval)
    setup_tfi_surface();
}
TrianglePatch::~TrianglePatch() 
{
    delete tfi_surface;
}
TrianglePatch* TrianglePatch::clone() const
{
    return new TrianglePatch(*this);
} 
TrianglePatch* TrianglePatch::copy() const
{
    return new TrianglePatch(*this);
} 
Vector3 TrianglePatch::eval( double r, double s ) const
{
    r = r0 + (r1 - r0) * r;
    s = s0 + (s1 - s0) * s;

    // Locate a query point on the CoonsPatch surface.
    Vector3 q = tfi_surface->eval(r, s);
    // Use TFI close to the boundaries.  
    double tol = 1.0e-12;
    if ( r < tol || s < tol || r > 1.0-tol || s > 1.0-tol ) {
	return q;
    }
    // Set up for projection.
    Vector3 point = q;
    Vector3 t1 = tfi_surface->dpdr(r, s);
    Vector3 t2 = tfi_surface->dpds(r, s);
    Vector3 qr = unit(cross(t1, t2));
    // Project onto planes defined by the set of triangular facets 
    // and stop when we land inside one of the triangles.
    int project_flag, inside_flag;
    for ( int i = 0; i < ntri; ++i ) {
	Vector3 qtest = q;
	Vector3 a = p[itri[i*3]];
	Vector3 b = p[itri[i*3+1]];
	Vector3 c = p[itri[i*3+2]];
	project_flag = project_onto_plane(qtest, qr, a, b, c);
	if ( project_flag != 0 && project_flag != 1 ) continue;
	inside_flag = inside_triangle(qtest, a, b, c);
	if ( inside_flag ) {
	    point = qtest;
	    break;
	}
    }
    return point;  // projected point on a facet
}
string TrianglePatch::str() const
{
    ostringstream ost;
    size_t i;
    ost << "TrianglePatch(";
    ost << "p=[";
    for ( i = 0; i < p.size()-1; ++i ) ost << p[i] << ", ";
    ost << p[p.size()-1] << "], ";
    ost << "itri=[";
    for ( i = 0; i < itri.size()-1; ++i ) ost << itri[i] << ", ";
    ost << itri[itri.size()-1] << "], ";
    ost << "iA=[";
    for ( i = 0; i < iA.size()-1; ++i ) ost << iA[i] << ", ";
    ost << iA[iA.size()-1] << "], ";
    ost << "iB=[";
    for ( i = 0; i < iB.size()-1; ++i ) ost << iB[i] << ", ";
    ost << iB[iB.size()-1] << "], ";
    ost << "iC=[";
    for ( i = 0; i < iC.size()-1; ++i ) ost << iC[i] << ", ";
    ost << iC[iC.size()-1] << "], ";
    ost << "iD=[";
    for ( i = 0; i < iD.size()-1; ++i ) ost << iD[i] << ", ";
    ost << iD[iD.size()-1] << "], ";
    ost << " ntri=" << ntri << ", ";
    ost << "\"" << label << "\"" << ", ";
    ost << r0 << ", " << r1 << ", " << s0 << ", " << s1 << ")";
    return ost.str();
}
TrianglePatch* TrianglePatch::translate( const Vector3 &v )
{
    for ( size_t i = 0; i < p.size(); ++i ) p[i] += v;
    tfi_surface->translate(v);
    return this;
}
TrianglePatch* TrianglePatch::translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    return this;
}
TrianglePatch* TrianglePatch::mirror_image( const Vector3 &point, const Vector3 &normal )
{
    for ( size_t i = 0; i < p.size(); ++i ) p[i].mirror_image(point, normal);
    tfi_surface->mirror_image(point, normal);
    return this;
}
TrianglePatch* TrianglePatch::rotate_about_zaxis( double dtheta )
{
    for ( size_t i = 0; i < p.size(); ++i ) p[i].rotate_about_zaxis(dtheta);
    tfi_surface->rotate_about_zaxis(dtheta);
    return this;
}
void TrianglePatch::setup_tfi_surface()
{
    size_t i;
    vector<Vector3*> pA, pB, pC, pD;
    for ( i = 0; i < iA.size(); ++i ) pA.push_back(&p[iA[i]]);
    Polyline cA = Polyline(pA, 0);
    for ( i = 0; i < iB.size(); ++i ) pB.push_back(&p[iB[i]]);
    Polyline cB = Polyline(pB, 0);
    for ( i = 0; i < iC.size(); ++i ) pC.push_back(&p[iC[i]]);
    Polyline cC = Polyline(pC, 0);
    for ( i = 0; i < iD.size(); ++i ) pD.push_back(&p[iD[i]]);
    Polyline cD = Polyline(pD, 0);
    tfi_surface = new CoonsPatch(cA, cB, cC, cD);
    return;
}
int TrianglePatch::add( const TrianglePatch &surf )
{
    int i; // Have changed to using int so that we can count down and
           // use the comparison >= 0
    int np0 = p.size();  // number of points originally in this patch
    if ( edges_are_matched( *this, iD, surf, surf.iC, 0 ) ) {
	//          B                  B
	//    p01--->---p11      p01--->---p11
	//     |         |        |         |
	//   C ^   [0]   ^ D    C ^   [1]   ^ D
	//     |         |        |         |
	//    p00--->---p10      p00--->---p10
	//          A                  A
	// Replace this edge D with the added-surface edge D
	iD.clear();
	for ( i = 0; i < (int)surf.iD.size(); ++i ) iD.push_back(np0 + surf.iD[i]);
	// Extend edges A and B.
	for ( i = 0; i < (int)surf.iA.size(); ++i ) iA.push_back(np0 + surf.iA[i]);
	for ( i = 0; i < (int)surf.iB.size(); ++i ) iB.push_back(np0 + surf.iB[i]);
    } else if ( edges_are_matched( *this, iD, surf, surf.iD, 1 ) ) {
	iD.clear();
	for ( i = (int)surf.iC.size()-1; i >= 0; --i ) iD.push_back(np0 + surf.iC[i]);
	for ( i = (int)surf.iB.size()-1; i >= 0; --i ) iA.push_back(np0 + surf.iB[i]);
	for ( i = (int)surf.iA.size()-1; i >= 0; --i ) iB.push_back(np0 + surf.iA[i]);
    } else if ( edges_are_matched( *this, iD, surf, surf.iA, 1 ) ) {
	iD.clear();
	for ( i = (int)surf.iB.size()-1; i >= 0; --i ) iD.push_back(np0 + surf.iB[i]);
	for ( i = 0; i < (int)surf.iD.size(); ++i ) iA.push_back(np0 + surf.iD[i]);
	for ( i = 0; i < (int)surf.iC.size(); ++i ) iB.push_back(np0 + surf.iC[i]);
    } else if ( edges_are_matched( *this, iD, surf, surf.iB, 0 ) ) {
	iD.clear();
	for ( i = 0; i < (int)surf.iA.size(); ++i ) iD.push_back(np0 + surf.iA[i]);
	for ( i = (int)surf.iD.size()-1; i >= 0; --i ) iA.push_back(np0 + surf.iD[i]);
	for ( i = (int)surf.iC.size()-1; i >= 0; --i ) iB.push_back(np0 + surf.iC[i]);
    } else if ( edges_are_matched( *this, iD, surf, surf.iC, 1 ) ) {
	//          B                  A
	//    p01--->---p11      p00--->---p10
	//     |         |        |         |
	//   C ^   [0]   ^ D    C v   [1]   v D
	//     |         |        |         |
	//    p00--->---p10      p01--->---p11
	//          A                  B
	// Replace this edge D with the added-surface edge D
	iD.clear();
	for ( i = (int)surf.iD.size()-1; i >= 0; --i ) iD.push_back(np0 + surf.iD[i]);
	// Extend edges A and B.
	for ( i = 0; i < (int)surf.iB.size(); ++i ) iA.push_back(np0 + surf.iB[i]);
	for ( i = 0; i < (int)surf.iA.size(); ++i ) iB.push_back(np0 + surf.iA[i]);
    } else if ( edges_are_matched( *this, iD, surf, surf.iD, 0 ) ) {
	iD.clear();
	for ( i = 0; i < (int)surf.iC.size(); ++i ) iD.push_back(np0 + surf.iC[i]);
	for ( i = (int)surf.iA.size()-1; i >= 0; --i ) iA.push_back(np0 + surf.iA[i]);
	for ( i = (int)surf.iB.size()-1; i >= 0; --i ) iB.push_back(np0 + surf.iB[i]);
    } else if ( edges_are_matched( *this, iD, surf, surf.iA, 0 ) ) {
	iD.clear();
	for ( i = 0; i < (int)surf.iB.size(); ++i ) iD.push_back(np0 + surf.iB[i]);
	for ( i = 0; i < (int)surf.iC.size(); ++i ) iA.push_back(np0 + surf.iC[i]);
	for ( i = 0; i < (int)surf.iD.size(); ++i ) iB.push_back(np0 + surf.iD[i]);
    } else if ( edges_are_matched( *this, iD, surf, surf.iB, 1 ) ) {
	iD.clear();
	for ( i = (int)surf.iA.size()-1; i >= 0; --i ) iD.push_back(np0 + surf.iA[i]);
	for ( i = (int)surf.iD.size()-1; i >= 0; --i ) iA.push_back(np0 + surf.iD[i]);
	for ( i = (int)surf.iC.size()-1; i >= 0; --i ) iB.push_back(np0 + surf.iC[i]);
    } else if ( edges_are_matched( *this, iB, surf, surf.iC, 1 ) ) {
	//          C                  B
	//    p00--->---p01      p01--->---p11
	//     |         |        |         |
	//   A v   [0]   v B    C ^   [1]   ^ D
	//     |         |        |         |
	//    p10--->---p11      p00--->---p10
	//          D                  A
	// Replace this edge B with the added-surface edge D
	iB.clear();
	for ( i = (int)surf.iD.size()-1; i >= 0; --i ) iB.push_back(np0 + surf.iD[i]);
	// Extend edges D and C.
	for ( i = 0; i < (int)surf.iA.size(); ++i ) iD.push_back(np0 + surf.iA[i]);
	for ( i = 0; i < (int)surf.iB.size(); ++i ) iC.push_back(np0 + surf.iB[i]);
    } else if ( edges_are_matched( *this, iB, surf, surf.iD, 0 ) ) {
	iB.clear();
	for ( i = 0; i < (int)surf.iC.size(); ++i ) iB.push_back(np0 + surf.iC[i]);
	for ( i = (int)surf.iB.size()-1; i >= 0; --i ) iD.push_back(np0 + surf.iB[i]);
	for ( i = (int)surf.iA.size()-1; i >= 0; --i ) iC.push_back(np0 + surf.iA[i]);
    } else if ( edges_are_matched( *this, iB, surf, surf.iA, 0 ) ) {
	iB.clear();
	for ( i = 0; i < (int)surf.iB.size(); ++i ) iB.push_back(np0 + surf.iB[i]);
	for ( i = 0; i < (int)surf.iD.size(); ++i ) iD.push_back(np0 + surf.iD[i]);
	for ( i = 0; i < (int)surf.iC.size(); ++i ) iC.push_back(np0 + surf.iC[i]);
    } else if ( edges_are_matched( *this, iB, surf, surf.iB, 1 ) ) {
	iB.clear();
	for ( i = (int)surf.iA.size()-1; i >= 0; --i ) iB.push_back(np0 + surf.iA[i]);
	for ( i = (int)surf.iC.size()-1; i >= 0; --i ) iD.push_back(np0 + surf.iC[i]);
	for ( i = (int)surf.iD.size()-1; i >= 0; --i ) iC.push_back(np0 + surf.iD[i]);
    } else if ( edges_are_matched( *this, iB, surf, surf.iC, 0 ) ) {
	//          C                  A
	//    p00--->---p01      p00--->---p10
	//     |         |        |         |
	//   A v   [0]   v B    C v   [1]   v D
	//     |         |        |         |
	//    p10--->---p11      p01--->---p11
	//          D                  B
	// Replace this edge B with the added-surface edge D
	iB.clear();
	for ( i = 0; i < (int)surf.iD.size(); ++i ) iB.push_back(np0 + surf.iD[i]);
	// Extend edges D and C.
	for ( i = 0; i < (int)surf.iB.size(); ++i ) iD.push_back(np0 + surf.iB[i]);
	for ( i = 0; i < (int)surf.iA.size(); ++i ) iC.push_back(np0 + surf.iA[i]);
    } else if ( edges_are_matched( *this, iB, surf, surf.iD, 1 ) ) {
	iB.clear();
	for ( i = (int)surf.iC.size()-1; i >= 0; --i ) iB.push_back(np0 + surf.iC[i]);
	for ( i = (int)surf.iA.size()-1; i >= 0; --i ) iD.push_back(np0 + surf.iA[i]);
	for ( i = (int)surf.iB.size()-1; i >= 0; --i ) iC.push_back(np0 + surf.iB[i]);
    } else if ( edges_are_matched( *this, iB, surf, surf.iA, 1 ) ) {
	iB.clear();
	for ( i = (int)surf.iB.size()-1; i >= 0; --i ) iB.push_back(np0 + surf.iB[i]);
	for ( i = 0; i < (int)surf.iC.size(); ++i ) iD.push_back(np0 + surf.iC[i]);
	for ( i = 0; i < (int)surf.iD.size(); ++i ) iC.push_back(np0 + surf.iD[i]);
    } else if ( edges_are_matched( *this, iB, surf, surf.iB, 0 ) ) {
	iB.clear();
	for ( i = 0; i < (int)surf.iA.size(); ++i ) iB.push_back(np0 + surf.iA[i]);
	for ( i = (int)surf.iD.size()-1; i >= 0; --i ) iD.push_back(np0 + surf.iD[i]);
	for ( i = (int)surf.iC.size()-1; i >= 0; --i ) iC.push_back(np0 + surf.iC[i]);
    } else {
	cout << "TrianglePatch::add() : no matching edges located." << endl;
	return 1;
    }
    // Append the new points and vertex indices to this TrianglePatch.
    // We will include redundant points from the common edge because
    // it makes the translation of indices easier.
    for ( i = 0; i < (int)surf.p.size(); ++i ) p.push_back(surf.p[i]);
    for ( i = 0; i < (int)surf.itri.size(); ++i ) itri.push_back(np0 + surf.itri[i]);
    ntri += surf.ntri;
    // Finally, rebuild the query surface.
    delete tfi_surface;
    setup_tfi_surface();
    return 0;
}

// helper function
bool edges_are_matched( const TrianglePatch &surf0, const vector<int> &edge0, 
			const TrianglePatch &surf1, const vector<int> &edge1,
			int opposite_direction )
{
    Vector3 p0, p1;
    int i1;
    double tolerance = 1.0e-9;
    // Test for point-by-point matching of the edges, 
    // returning a fail at the first difference.
    if ( edge0.size() != edge1.size() ) return 0;
    for ( size_t i = 0; i < edge0.size(); ++i ) {
	p0 = surf0.p[edge0[i]];
	if ( opposite_direction ) {
	    i1 = edge1.size() - 1 - i;
	} else {
	    i1 = i;
	}
	p1 = surf1.p[edge1[i1]];
	if ( vabs(p1 - p0) > tolerance ) return 0;  
    }
    // We have managed to get past all of the tests, so assume it is a match.
    return 1;
}

//----------------------------------------------------------------------------------

BezierPatch::BezierPatch( const vector<Vector3*> &_Q, int _n, int _m, string label,
			  double r0, double r1, double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), Q(vector<Vector3>()), n(_n), m(_m) 
{
    for ( size_t i = 0; i < _Q.size(); ++i ) Q.push_back(*_Q[i]);
    if ( (n+1)*(m+1) != (int)Q.size() ) {
	cout << "BezierPatch error: control net is not correct size." << endl;
	cout << "                   n=" << n << ", m=" << m;
	cout << ", Q.size()=" << Q.size() << endl;
    }
}
BezierPatch::BezierPatch( const BezierPatch &surf )
    : ParametricSurface(surf.label, surf.r0, surf.r1, surf.s0, surf.s1), 
      Q(vector<Vector3>()), n(surf.n), m(surf.m) 
{
    for ( size_t i = 0; i < surf.Q.size(); ++i ) Q.push_back(surf.Q[i]);
}
BezierPatch::~BezierPatch() {}
BezierPatch* BezierPatch::clone() const
{
    return new BezierPatch(*this);
} 
BezierPatch* BezierPatch::copy() const
{
    return new BezierPatch(*this);
} 
Vector3 BezierPatch::eval( double r, double s ) const 
{
    r = r0 + (r1 - r0) * r;
    s = s0 + (s1 - s0) * s;

    vector<Vector3> Qi;
    for ( int i = 0; i <= n; ++i ) {
	vector<Vector3> Qij;
	for ( int j = 0; j <= m; ++j ) {
	    Qij.push_back(Q[i*(m+1) + j]);
	}
	Bezier b = Bezier(Qij);
	Qi.push_back(b.eval(s));
    }
    Bezier b2 = Bezier(Qi);
    return b2.eval(r);
}
string BezierPatch::str() const
{
    ostringstream ost;
    ost << "BezierPatch(" << "[";
    for ( size_t i = 0; i < Q.size(); ++i )
	ost << Q[i].str() << (( i < Q.size()-1 ) ? ", " : "]");
    ost << ", n=" << n << ", m=" << m << ", ";
    ost << "\"" << label << "\"" << ", ";
    ost << r0 << ", " << r1 << ", " << s0 << ", " << s1 << ")";
    return ost.str();
}
BezierPatch* BezierPatch::translate( const Vector3 &v )
{
    for ( size_t i = 0; i < Q.size(); ++i ) Q[i] += v;
    return this;
}
BezierPatch* BezierPatch::translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    return this;
}
BezierPatch* BezierPatch::mirror_image( const Vector3 &point, const Vector3 &normal )
{
    for ( size_t i = 0; i < Q.size(); ++i ) Q[i].mirror_image(point, normal);
    return this;
}
BezierPatch* BezierPatch::rotate_about_zaxis( double dtheta )
{
    for ( size_t i = 0; i < Q.size(); ++i ) Q[i].rotate_about_zaxis(dtheta);
    return this;
}

//----------------------------------------------------------------------------------

RevolvedSurface::RevolvedSurface( const Path* pathptr, string label,
				  double r0, double r1, double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), path(pathptr->clone()) {}
RevolvedSurface::RevolvedSurface( const RevolvedSurface &surf )
    : ParametricSurface(surf.label, surf.r0, surf.r1, surf.s0, surf.s1), 
      path(surf.path->clone()) {}
RevolvedSurface::~RevolvedSurface()
{
    delete path;
}
RevolvedSurface* RevolvedSurface::clone() const
{
    return new RevolvedSurface(*this);
} 
RevolvedSurface* RevolvedSurface::copy() const
{
    return new RevolvedSurface(*this);
} 
Vector3 RevolvedSurface::eval( double r, double s ) const
{
    r = r0 + (r1 - r0) * r;
    s = s0 + (s1 - s0) * s;
    // r is the parametric distance along the surface
    Vector3 p = path->eval(r);
    // coordinates in the y,z-plane will be rotated
    double radius = sqrt(p.y * p.y + p.z * p.z);
    double theta = atan2(p.z, p.y);
    // 0 <= s <= 1.0 is mapped to angle 0 <= theta <= 2 pi, in radians
    theta += (s * 2.0 * M_PI);
    p.y = radius * cos(theta);
    p.z = radius * sin(theta);
    return p;
}
string RevolvedSurface::str() const
{
    ostringstream ost;
    ost << "RevolvedSurface(" << path->str() << ", ";
    ost << "\"" << label << "\"" << ", ";
    ost << r0 << ", " << r1 << ", " << s0 << ", " << s1 << ")";
    return ost.str();
}

//----------------------------------------------------------------------------------

class ErrorFunctionForMappedSurface: public MultivariateFunction {
public:
    Vector3 query_point;
    Vector3 query_vector;
    ParametricSurface *surf;
    ErrorFunctionForMappedSurface( Vector3 &qp, Vector3 &qr, ParametricSurface *true_surf )
	: query_point(qp), query_vector(qr), surf(true_surf) {}
    double eval(vector<double> &x)
    {
	double r = x[0];
	double s = x[1];
	// Keep within the nominal parameter range by applying a big penalty.
	if ( r < 0.0 || r > 1.0 || s < 0.0 || s > 1.0 ) return 1.0e6;
	Vector3 p = surf->eval(r, s);
	// The true surface point will be aligned with the query vector.
	double err = vabs(cross(query_vector, p - query_point));
	return err;
    }
};

MappedSurface::MappedSurface( const ParametricSurface* qsptr, 
			      const ParametricSurface* tsptr,
			      string label, double r0, double r1,
			      double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), 
      query_surf(qsptr->clone()),
      true_surf(tsptr->clone()) {}
MappedSurface::MappedSurface( const MappedSurface &surf )
    : ParametricSurface(surf.label, surf.r0, surf.r1, surf.s0, surf.s1), 
      query_surf(surf.query_surf->clone()),
      true_surf(surf.true_surf->clone()) {}
MappedSurface::~MappedSurface()
{
    delete query_surf;
    delete true_surf;
}
MappedSurface* MappedSurface::clone() const
{
    return new MappedSurface(*this);
} 
MappedSurface* MappedSurface::copy() const
{
    return new MappedSurface(*this);
} 
Vector3 MappedSurface::eval( double r, double s ) const
{
    r = r0 + (r1 - r0) * r;
    s = s0 + (s1 - s0) * s;
    Vector3 q = query_surf->eval(r,s);
    Vector3 t1 = unit(query_surf->dpdr(r,s));
    Vector3 t2 = unit(query_surf->dpds(r,s));
    Vector3 qr = cross(t1, t2);
    Vector3 p;  // point on true surface
    vector<double> x = vector<double>(2);
    double f_min = 0.0;
    int n_fe = 0;
    int n_restart = 0;
    int flag = 0;
    ErrorFunctionForMappedSurface error_f = ErrorFunctionForMappedSurface(q, qr, true_surf);
    // Some functions can have local minima nowhere near where we want to end up.
    // We should have a mesh of start points and only quit when the minimizer 
    // has found a good solution.
    int nr = 5; int ns = 5;
    double dr = 1.0 / nr;
    double ds = 1.0 / ns;
    for ( int ir = 0; ir < nr; ++ir ) {
	for ( int is = 0; is < ns; ++is ) {
	    x[0] = dr/2 + ir * dr;
	    x[1] = ds/2 + is * ds;
	    flag = minimize(&error_f, x, &f_min, &n_fe, &n_restart);
	    p = true_surf->eval(x[0], x[1]);
	    if ( fabs(f_min) < 1.0e-6 ) break;
	}
	if ( fabs(f_min) < 1.0e-6 ) break;
    }
    if ( fabs(f_min) > 1.0e-6 ) {
	cout << "MappedSurface::eval() after minimizer: f_min=" << f_min << endl;
	cout << "    flag=" << flag << ", n_fe=" << n_fe << ", n_restart=" << n_restart << endl;
	cout << "    q=" << q << ", qr=" << qr << endl;
	cout << "    x[0]=" << x[0] << ", x[1]=" << x[1] << ", p=" << p << endl;
    }
    return p;
}
string MappedSurface::str() const
{
    ostringstream ost;
    ost << "MappedSurface(";
    ost << query_surf->str() << ", ";
    ost << true_surf->str() << ", ";
    ost << "\"" << label << "\"" << ", ";
    ost << r0 << ", " << r1 << ", " << s0 << ", " << s1 << ")";
    return ost.str();
}

// ---------------------------------------------------------------------------------------

PolarSurface::PolarSurface( const ParametricSurface* surf, 
			    double H, string label,
			    double r0, double r1, double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), 
      original_surf(surf->clone()), H(H)
{}
PolarSurface::PolarSurface( const PolarSurface &surf )
    : ParametricSurface(surf.label, surf.r0, surf.r1, surf.s0, surf.s1), 
      original_surf(surf.original_surf->clone()), H(surf.H)
{}
PolarSurface::~PolarSurface()
{
    delete original_surf;
}
PolarSurface* PolarSurface::clone() const
{
    return new PolarSurface(*this);
}
PolarSurface* PolarSurface::copy() const
{
    return new PolarSurface(*this);
}
Vector3 PolarSurface::eval( double r, double s ) const
{
    Vector3 p = original_surf->eval(r, s);
    map_neutral_plane_to_cylinder(p, H);
    return p;
}
string PolarSurface::str() const
{
    ostringstream ost;
    ost << "PolarSurface(";
    ost << original_surf->str() << ", ";
    ost << "\"" << label << "\"" << ", ";
    ost << r0 << ", " << r1 << ", " << s0 << ", " << s1 << ")";
    return ost.str();
}

// ---------------------------------------------------------------------------------------

SurfaceThruVolume::SurfaceThruVolume( const ParametricVolume& _pvol,
				      const BivariateFunction& _fr,
				      const BivariateFunction& _fs,
				      const BivariateFunction& _ft,
				      string label,
				      double r0, double r1, double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), 
      pvol(_pvol.clone()), fr(_fr.clone()), fs(_fs.clone()), ft(_ft.clone())
{}
SurfaceThruVolume::SurfaceThruVolume( const SurfaceThruVolume &surf )
    : ParametricSurface(surf.label, surf.r0, surf.r1, surf.s0, surf.s1), 
      pvol(surf.pvol->clone()), 
      fr(surf.fr->clone()), fs(surf.fs->clone()), ft(surf.ft->clone())
{}
SurfaceThruVolume::~SurfaceThruVolume()
{
    delete pvol;
    delete fr;
    delete fs;
    delete ft;
}
SurfaceThruVolume* SurfaceThruVolume::clone() const
{
    return new SurfaceThruVolume(*this);
}
SurfaceThruVolume* SurfaceThruVolume::copy() const
{
    return new SurfaceThruVolume(*this);
}
Vector3 SurfaceThruVolume::eval( double r, double s ) const
{
    Vector3 p = pvol->eval(fr->eval(r,s), fs->eval(r,s), ft->eval(r,s));
    return p;
}
string SurfaceThruVolume::str() const
{
    ostringstream ost;
    ost << "SurfaceThruVolume(";
    ost << pvol->str() << ", ";
    ost << fr->str() << ", ";
    ost << fs->str() << ", ";
    ost << ft->str() << ", ";
    ost << "\"" << label << "\"" << ", ";
    ost << r0 << ", " << r1 << ", " << s0 << ", " << s1 << ")";
    return ost.str();
}

NurbsSurface::
NurbsSurface( const vector<vector<Vector3> > &Q, const vector<vector<double> > &w,
	      int p, const vector<double> &U, int q, const vector<double> &V,
	      string label,
	      double r0, double r1,
	      double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), Q(Q), w(w), p(p), U(U),
      q(q), V(V)
{
    umin_ = U.front();
    umax_ = U.back();
    vmin_ = V.front();
    vmax_ = V.back();

    Qw_.resize(Q.size());

    for ( size_t i = 0; i < Q.size(); ++i ) {

	Qw_[i].resize(Q[0].size());

	for ( size_t j = 0; j < Q[0].size(); ++j ) {
	    Qw_[i][j].wx = w[i][j] * Q[i][j].x;
	    Qw_[i][j].wy = w[i][j] * Q[i][j].y;
	    Qw_[i][j].wz = w[i][j] * Q[i][j].z;
	    Qw_[i][j].w = w[i][j];
	}
    }

}

NurbsSurface::
NurbsSurface( const vector<vector<Mapped_point> > &Qw,
	      int p, const vector<double> &U, int q, const vector<double> &V,
	      string label,
	      double r0, double r1,
	      double s0, double s1 )
    : ParametricSurface(label, r0, r1, s0, s1), p(p), U(U),
      q(q), V(V), Qw_(Qw)
{
    umin_ = U.front();
    umax_ = U.back();
    vmin_ = V.front();
    vmax_ = V.back();

    Q.resize(Qw_.size());
    w.resize(Qw_.size());

    for ( size_t i = 0; i < Q.size(); ++i ) {

	Q[i].resize(Qw_[0].size());
	w[i].resize(Qw_[0].size());

	for ( size_t j = 0; j < Q[0].size(); ++j ) {
	    if ( Qw_[i][j].w == 0.0 ) {
		Q[i][j].x = Qw_[i][j].wx;
		Q[i][j].y = Qw_[i][j].wy;
		Q[i][j].z = Qw_[i][j].wz;
		w[i][j] = Qw_[i][j].w;
	    }
	    else {
		Q[i][j].x = Qw_[i][j].wx / Qw_[i][j].w;
		Q[i][j].y = Qw_[i][j].wy / Qw_[i][j].w;
		Q[i][j].z = Qw_[i][j].wz / Qw_[i][j].w;
		w[i][j] = Qw_[i][j].w;
	    }
	}
    }
}

NurbsSurface::
NurbsSurface( const NurbsSurface &n)
    : ParametricSurface(n.label, n.r0, n.r1, n.s0, n.s1),
      Q(n.Q), w(n.w), p(n.p), U(n.U), q(n.q), V(n.V),
      Qw_(n.Qw_), umin_(n.umin_), umax_(n.umax_),
      vmin_(n.vmin_), vmax_(n.vmax_) {}


NurbsSurface::
~NurbsSurface() {}

NurbsSurface*
NurbsSurface::
clone() const
{
    return new NurbsSurface(*this);
}

NurbsSurface*
NurbsSurface::
copy() const
{
    return new NurbsSurface(*this);
}

string
NurbsSurface::
str() const
{
    ostringstream ost;
    ost << "NurbsSurface::str(): NOT IMPLEMENTED PROPERLY.\n";
    return ost.str();
}
NurbsSurface*
NurbsSurface::
translate( const Vector3 &v )
{
    for ( size_t i = 0; i < Q.size(); ++i ) {
	for ( size_t j = 0; j < Q[0].size(); ++j ) {
	    Q[i][j] += v;
	    Qw_[i][j].wx = w[i][j] * Q[i][j].x;
	    Qw_[i][j].wy = w[i][j] * Q[i][j].y;
	    Qw_[i][j].wz = w[i][j] * Q[i][j].z;
	    Qw_[i][j].w = w[i][j];
	}
    }
    return this;
}

NurbsSurface*
NurbsSurface::
translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    return this;
}

NurbsSurface*
NurbsSurface::
mirror_image( const Vector3 &point, const Vector3 &normal )
{
    for ( size_t i = 0; i < Q.size(); ++i ) {
	for ( size_t j = 0; j < Q[0].size(); ++j ) {
	    Q[i][j].mirror_image(point, normal);
	    Qw_[i][j].wx = w[i][j] * Q[i][j].x;
	    Qw_[i][j].wy = w[i][j] * Q[i][j].y;
	    Qw_[i][j].wz = w[i][j] * Q[i][j].z;
	    Qw_[i][j].w = w[i][j];
	}
    }
    return this;
}

double
NurbsSurface::
map_r2u(double r) const
{
    return umin_ + r*(umax_ - umin_);
}

double
NurbsSurface::
map_s2v(double s) const
{
    return vmin_ + s*(vmax_ - vmin_);
}


Vector3
NurbsSurface::
eval( double r, double s ) const
{
    double u = map_r2u(r);
    double v = map_s2v(s);
    return nurbs_surface_point(u, p, U, v, q, V, Qw_);
}


int write_STL(const ParametricSurface &surf, int nr, int ns, string fname)
{
    // s: surface
    // nr: number of 'cells' in r direction
    // ns: number of 'cells' in s direction
    // fname: output name for STL surface
    //
    // This function first divides the surface into quad-shaped
    // cells. Each of these is subdivided to give the triangles
    // required for STL output.

    ofstream f(fname.c_str());
    if ( f.fail() ) {
	cout << "Problem opening file: " << fname << " for writing in write_STL().\n";
	cout << "Bailing out!\n";
	exit(1);
    }

    f << setprecision(12) << scientific;
    f << "solid dummy-name\n";

    double r = 0.0;
    double dr = 1.0/nr;
    double s = 0.0;
    double ds = 1.0/ns;
    Vector3 A, B, C, D, n1, n2;
    
    for ( int ir = 0; ir < nr; ++ir ) {
	r = ir*dr;
	for ( int is = 0; is < ns; ++is ) {
	    s = is*ds;
	    A = surf.eval(r, s);
	    B = surf.eval(r+dr, s);
	    C = surf.eval(r+dr, s+ds);
	    D = surf.eval(r, s+ds);
	    
	    // Tri: ADC
	    n1 = unit(cross(D-A, C-D));
	    f << "facet normal ";
	    f << setw(20) << n1.x << " ";
	    f << setw(20) << n1.y << " ";
	    f << setw(20) << n1.z << endl;
	    f << "outer loop\n";
	    f << "vertex ";
	    f << setw(20) << A.x << " ";
	    f << setw(20) << A.y << " ";
	    f << setw(20) << A.z << endl;
	    f << "vertex ";
	    f << setw(20) << D.x << " ";
	    f << setw(20) << D.y << " ";
	    f << setw(20) << D.z << endl;
	    f << "vertex ";
	    f << setw(20) << C.x << " ";
	    f << setw(20) << C.y << " ";
	    f << setw(20) << C.z << endl;
	    f << "endloop\n";
	    f << "endfacet\n";

	    // Tri: BAC
	    n2 = unit(cross(A-B, C-A));
	    f << "facet normal ";
	    f << setw(20) << n2.x << " ";
	    f << setw(20) << n2.y << " ";
	    f << setw(20) << n2.z << endl;
	    f << "outer loop\n";
	    f << "vertex ";
	    f << setw(20) << B.x << " ";
	    f << setw(20) << B.y << " ";
	    f << setw(20) << B.z << endl;
	    f << "vertex ";
	    f << setw(20) << A.x << " ";
	    f << setw(20) << A.y << " ";
	    f << setw(20) << A.z << endl;
	    f << "vertex ";
	    f << setw(20) << C.x << " ";
	    f << setw(20) << C.y << " ";
	    f << setw(20) << C.z << endl;
	    f << "endloop\n";
	    f << "endfacet\n";
	}
    }

    f << "endsolid dummy-name\n";
    f.close();

    return 0;
}
