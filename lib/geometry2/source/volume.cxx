/** \file volume.cxx
 *  \ingroup libgeom2
 *  \brief Implementation of the hexahedral volume class. 
 *  \author PJ
 *  \version 12-Jan-2006 initial coding
 *
 * This module is a C++ replacement for the combined C+Python geometry code that
 * was built for mb_cns and Elmer.
 */
#include <cstdio>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include "geom.hh"
#include "gpath.hh"
#include "../../nm/source/secant.hh"
#include "surface.hh"
#include "volume.hh"
using namespace std;

// Base class for geometric surfaces.
ParametricVolume::ParametricVolume( const ParametricSurface* faceN, 
				    const ParametricSurface* faceE,
				    const ParametricSurface* faceS,
				    const ParametricSurface* faceW,
				    const ParametricSurface* faceT,
				    const ParametricSurface* faceB,
				    string label,
				    double r0, double r1, double s0, double s1,
				    double t0, double t1 )
    : label(label), r0(r0), r1(r1), s0(s0), s1(s1), t0(t0), t1(t1)
{
    // We clone the boundary surfaces to ensure their persistence.
    // Beware that the pointer values may be NULL.
    if ( faceS ) south = faceS->clone();
    if ( faceB ) bottom = faceB->clone();
    if ( faceW ) west = faceW->clone();
    if ( faceE ) east = faceE->clone(); 
    if ( faceN ) north = faceN->clone();
    if ( faceT ) top = faceT->clone();
    set_and_check_corners();
}
// Faces supplied as a vector of pointers with old (confusing) order.
ParametricVolume::ParametricVolume( const vector<ParametricSurface*> &face, string label,
				    double r0, double r1, double s0, double s1,
				    double t0, double t1 )
    : label(label), r0(r0), r1(r1), s0(s0), s1(s1), t0(t0), t1(t1)
{
    // We clone the boundary surfaces to ensure their persistence.
    // Beware that the pointer values may be NULL.
    if ( face[0] ) south = face[0]->clone();
    if ( face[1] ) bottom = face[1]->clone();
    if ( face[2] ) west = face[2]->clone();
    if ( face[3] ) east = face[3]->clone(); 
    if ( face[4] ) north = face[4]->clone();
    if ( face[5] ) top = face[5]->clone();
    set_and_check_corners();
}
// Derived classes will set up the surfaces themselves.
ParametricVolume::ParametricVolume( string label, double r0, double r1,
				    double s0, double s1, double t0, double t1 )
    : label(label), 
      south(0), bottom(0), west(0), east(0), north(0), top(0), 
      r0(r0), r1(r1), s0(s0), s1(s1), t0(t0), t1(t1) {}
ParametricVolume::ParametricVolume( const ParametricVolume &vol )
    : label(vol.label),
      south(vol.south->clone()), bottom(vol.bottom->clone()),
      west(vol.west->clone()), east(vol.east->clone()), 
      north(vol.north->clone()), top(vol.top->clone()),
      r0(vol.r0), r1(vol.r1), s0(vol.s0), s1(vol.s1), t0(vol.t0), t1(vol.t1) 
{
    set_and_check_corners();
}
ParametricVolume::~ParametricVolume() 
{
    delete south; delete bottom; delete west;
    delete east; delete north; delete top;
}
ParametricVolume* ParametricVolume::clone() const
{
    return new ParametricVolume(*this);
}
ParametricVolume* ParametricVolume::copy() const
{
    return new ParametricVolume(*this);
}
Vector3 ParametricVolume::eval( double r, double s, double t ) const 
{ 
    // Transform to parameter range of subsection of the volume.
    r = r0 + (r1 - r0) * r;
    s = s0 + (s1 - s0) * s;
    t = t0 + (t1 - t0) * t;
    // Locate a point within the volume by blended linear interpolation.
    // @param r: interpolation parameter i-direction west-->east, 0.0<=r<=1.0
    // @param s: interpolation parameter j-direction south-->north, 0.0<=s<=1.0 
    // @param t: interpolation parameter k-direction bottom-->top, 0.0<=t<=1.0
    // @returns: a L{Vector3} value for the point.
    Vector3 pW = west->eval(s,t); 
    Vector3 pE = east->eval(s,t);
    Vector3 pS = south->eval(r,t); 
    Vector3 pN = north->eval(r,t);
    Vector3 pB = bottom->eval(r,s);
    Vector3 pT = top->eval(r,s);
    double omr = 1.0 - r; double oms = 1.0 - s; double omt = 1.0 - t;
    Vector3 BigC = omr * oms * omt * p000 + omr * oms * t * p001 +
	omr * s * omt * p010   + omr * s * t * p011 +
	r * oms * omt * p100   + r * oms * t * p101 +
	r * s * omt * p110     + r * s * t * p111;
    Vector3 p = 0.5 * ( omr * pW + r * pE +
                        oms * pS + s * pN +
                        omt * pB + t * pT ) - 0.5 * BigC;
    return p; 
}
string ParametricVolume::str() const
{
    ostringstream ost;
    ost << "ParametricVolume("
	<< "[ " << south->str()
	<< ", " << bottom->str()
	<< ", " << west->str()
	<< ", " << east->str()
	<< ", " << north->str()
	<< ", " << top->str()
	<< "], \"" << label << ", ";
    ost << r0 << ", " << r1 << ", " 
	<< s0 << ", " << s1 << ", " 
	<< t0 << ", " << t1 << "\")";
    return ost.str();
}
ParametricVolume* ParametricVolume::translate( const Vector3 &v ) 
{
    south->translate(v);
    bottom->translate(v);
    west->translate(v);
    east->translate(v);
    north->translate(v);
    top->translate(v);
    set_and_check_corners();
    return this;
}
ParametricVolume* ParametricVolume::translate( double vx, double vy, double vz ) 
{
    translate(Vector3(vx, vy, vz));
    return this;
}
ParametricVolume* ParametricVolume::mirror_image( const Vector3 &point, const Vector3 &normal ) 
{
    south->mirror_image(point, normal);
    bottom->mirror_image(point, normal);
    west->mirror_image(point, normal);
    east->mirror_image(point, normal);
    north->mirror_image(point, normal);
    top->mirror_image(point, normal);
    set_and_check_corners();
    return this;
}
ParametricVolume* ParametricVolume::rotate_about_zaxis( double dtheta ) 
{
    south->rotate_about_zaxis(dtheta);
    bottom->rotate_about_zaxis(dtheta);
    west->rotate_about_zaxis(dtheta);
    east->rotate_about_zaxis(dtheta);
    north->rotate_about_zaxis(dtheta);
    top->rotate_about_zaxis(dtheta);
    set_and_check_corners();
    return this;
}
int ParametricVolume::set_and_check_corners()
{
    // Compute the corners (for later use in TFI).
    p000 = bottom->eval(0.0, 0.0);
    p100 = bottom->eval(1.0, 0.0);
    p110 = bottom->eval(1.0, 1.0);
    p010 = bottom->eval(0.0, 1.0);
    p001 = top->eval(0.0, 0.0);
    p101 = top->eval(1.0, 0.0);
    p111 = top->eval(1.0, 1.0);
    p011 = top->eval(0.0, 1.0);
    // Check that the corners close by evaluating the same points
    // with east and west surfaces.
    Vector3 p000a = west->eval(0.0, 0.0);
    Vector3 p100a = east->eval(0.0, 0.0);
    Vector3 p110a = east->eval(1.0, 0.0);
    Vector3 p010a = west->eval(1.0, 0.0);
    Vector3 p001a = west->eval(0.0, 1.0);
    Vector3 p101a = east->eval(0.0, 1.0);
    Vector3 p111a = east->eval(1.0, 1.0);
    Vector3 p011a = west->eval(1.0, 1.0);
    // ... and with north, south surfaces.
    Vector3 p000b = south->eval(0.0, 0.0);
    Vector3 p100b = south->eval(1.0, 0.0);
    Vector3 p110b = north->eval(1.0, 0.0);
    Vector3 p010b = north->eval(0.0, 0.0);
    Vector3 p001b = south->eval(0.0, 1.0);
    Vector3 p101b = south->eval(1.0, 1.0);
    Vector3 p111b = north->eval(1.0, 1.0);
    Vector3 p011b = north->eval(0.0, 1.0);
    double average_length = (vabs(p100-p000) + vabs(p010-p000) +
			     vabs(p110-p100) + vabs(p110-p010) +
			     vabs(p101-p001) + vabs(p011-p001) +
			     vabs(p111-p101) + vabs(p111-p011) +
			     vabs(p001-p000) + vabs(p101-p100) +
			     vabs(p111-p110) + vabs(p011-p010)) / 12.0;
    double tol = 1.0e-6 * (average_length + 1.0);
    if ( vabs(p000-p000a) > tol || vabs(p000-p000b) > tol ) {
	cout << "p000 corner appears open: " 
	     << p000 << ", " << p000a << ", " << p000b << endl;
    }
    if ( vabs(p100-p100a) > tol || vabs(p100-p100b) > tol ) {
	cout << "p100 corner appears open: " 
	     << p100 << ", " << p100a << ", " << p100b << endl;
    }
    if ( vabs(p010-p010a) > tol || vabs(p010-p010b) > tol ) {
	cout << "p010 corner appears open: " 
	     << p010 << ", " << p010a << ", " << p010b << endl;
    }
    if ( vabs(p110-p110a) > tol || vabs(p110-p110a) > tol ) {
	cout << "p110 corner appears open: " 
	     << p110 << ", " << p110a << ", " << p110b << endl;
    }
    if ( vabs(p001-p001a) > tol || vabs(p001-p001b) > tol ) {
	cout << "p001 corner appears open: " 
	     << p001 << ", " << p001a << ", " << p001b << endl;
    }
    if ( vabs(p101-p101a) > tol || vabs(p101-p101b) > tol ) {
	cout << "p101 corner appears open: " 
	     << p101 << ", " << p101a << ", " << p101b << endl;
    }
    if ( vabs(p011-p011a) > tol || vabs(p011-p011b) > tol ) {
	cout << "p011 corner appears open: " 
	     << p011 << ", " << p011a << ", " << p011b << endl;
    }
    if ( vabs(p111-p111a) > tol || vabs(p111-p111a) > tol ) {
	cout << "p111 corner appears open: " 
	     << p111 << ", " << p111a << ", " << p111b << endl;
    }
    return 0;
}

// Helper functions...

// Overload stream output for Surface objects
ostream& operator<<( ostream &os, const ParametricVolume &vol )
{
    os << vol.str();
    return os;
}

// Special cases...

WireFrameVolume::WireFrameVolume( const Path &c01, const Path &c12, const Path &c32, const Path &c03, 
				  const Path &c45, const Path &c56, const Path &c76, const Path &c47,
				  const Path &c04, const Path &c15, const Path &c26, const Path &c37,
				  const string label, double r0, double r1, double s0, double s1,
				  double t0, double t1 )
    : ParametricVolume(label, r0, r1, s0, s1, t0, t1)
{
    // Make the boundary surfaces from the 12 edges.
    south  = new CoonsPatch(c01, c45, c04, c15);
    bottom = new CoonsPatch(c01, c32, c03, c12);
    west   = new CoonsPatch(c03, c47, c04, c37);
    east   = new CoonsPatch(c12, c56, c15, c26);
    north  = new CoonsPatch(c32, c76, c37, c26);
    top    = new CoonsPatch(c45, c76, c47, c56);
    set_and_check_corners();
}
WireFrameVolume::WireFrameVolume( const CoonsPatch &base_surf, 
				  Path &extrude_path, 
				  const string direction,
				  const string label, 
				  double r0, double r1, 
				  double s0, double s1,
				  double t0, double t1 )
    : ParametricVolume(label, r0, r1, s0, s1, t0, t1)
{
    Path *c01, *c32, *c03, *c12;
    Path *c45, *c76, *c47, *c56;
    Path *c04, *c15, *c26, *c37;
    Vector3 original_p0, delta;
    Vector3 new_p0, new_p1, new_p2, new_p3;
    Vector3 new_p4, new_p5, new_p6, new_p7;
    // Make the boundary surfaces with the supplied base surface and
    // an extrusion Path.
    original_p0 = base_surf.eval(0.0, 0.0);
    new_p0 = extrude_path.eval(0.0);
    delta = new_p0 - original_p0;
    if ( direction == "k" || direction == "+k" ) {
	// The base surface becomes the BOTTOM surface and the
	// rest of the block is extruded in the positive-k direction.
	c01 = base_surf.cA->clone(); c01->translate(delta);
	c32 = base_surf.cB->clone(); c32->translate(delta);
	c03 = base_surf.cC->clone(); c03->translate(delta);
	c12 = base_surf.cD->clone(); c12->translate(delta);
	//
	new_p4 = extrude_path.eval(1.0);
	delta = new_p4 - new_p0;
	c45 = c01->clone(); c45->translate(delta);
	c76 = c32->clone(); c76->translate(delta);
	c47 = c03->clone(); c47->translate(delta);
	c56 = c12->clone(); c56->translate(delta);
	// connecting lines
	new_p1 = c01->eval(1.0);
	new_p2 = c32->eval(1.0);
	new_p3 = c32->eval(0.0);
	c04 = extrude_path.clone();
	c15 = extrude_path.clone(); c15->translate(new_p1 - new_p0);
	c26 = extrude_path.clone(); c26->translate(new_p2 - new_p0);
	c37 = extrude_path.clone(); c37->translate(new_p3 - new_p0);
    } else if ( direction == "j" || direction == "+j" ) {
	// The base surface becomes the SOUTH surface and the
	// rest of the block is extruded in the positive-j direction.
	c01 = base_surf.cA->clone(); c01->translate(delta);
	c45 = base_surf.cB->clone(); c45->translate(delta);
	c04 = base_surf.cC->clone(); c04->translate(delta);
	c15 = base_surf.cD->clone(); c15->translate(delta);
	//
	new_p3 = extrude_path.eval(1.0);
	delta = new_p3 - new_p0;
	c32 = c01->clone(); c32->translate(delta);
	c76 = c45->clone(); c76->translate(delta);
	c37 = c04->clone(); c37->translate(delta);
	c26 = c15->clone(); c26->translate(delta);
	// connecting lines
	new_p1 = c01->eval(1.0);
	new_p5 = c15->eval(1.0);
	new_p4 = c04->eval(1.0);
	c03 = extrude_path.clone();
	c12 = extrude_path.clone(); c12->translate(new_p1 - new_p0);
	c56 = extrude_path.clone(); c56->translate(new_p5 - new_p0);
	c47 = extrude_path.clone(); c47->translate(new_p4 - new_p0);
    } else if ( direction == "i" || direction == "+i" ) {
	c03 = base_surf.cA->clone(); c03->translate(delta);
	c47 = base_surf.cB->clone(); c47->translate(delta);
	c04 = base_surf.cC->clone(); c04->translate(delta);
	c37 = base_surf.cD->clone(); c37->translate(delta);
	//
	new_p1 = extrude_path.eval(1.0);
	delta = new_p1 - new_p0;
	c12 = c03->clone(); c12->translate(delta);
	c56 = c47->clone(); c56->translate(delta);
	c15 = c04->clone(); c15->translate(delta);
	c26 = c37->clone(); c26->translate(delta);
	// connecting lines
	new_p3 = c03->eval(1.0);
	new_p4 = c04->eval(1.0);
	new_p7 = c37->eval(1.0);
	c01 = extrude_path.clone();
	c32 = extrude_path.clone(); c32->translate(new_p3 - new_p0);
	c76 = extrude_path.clone(); c76->translate(new_p7 - new_p0);
	c45 = extrude_path.clone(); c45->translate(new_p4 - new_p0);
    } else {
	cout << "Requested extrusion is not implemented." << endl;
	c01 = 0; c32 = 0; c03 = 0; c12 = 0;
	c45 = 0; c76 = 0; c47 = 0; c56 = 0;
	c04 = 0; c15 = 0; c26 = 0; c37 = 0;
    }
    south  = new CoonsPatch(*c01, *c45, *c04, *c15);
    bottom = new CoonsPatch(*c01, *c32, *c03, *c12);
    west   = new CoonsPatch(*c03, *c47, *c04, *c37);
    east   = new CoonsPatch(*c12, *c56, *c15, *c26);
    north  = new CoonsPatch(*c32, *c76, *c37, *c26);
    top    = new CoonsPatch(*c45, *c76, *c47, *c56);
    set_and_check_corners();
}


SimpleBoxVolume::SimpleBoxVolume( const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3,
				  const Vector3 &p4, const Vector3 &p5, const Vector3 &p6, const Vector3 &p7, 
				  const string label, double r0, double r1, double s0, double s1,
				  double t0, double t1 )
    : ParametricVolume(label, r0, r1, s0, s1, t0, t1)
{
    // cout << "Construct a Volume from the corner points." << endl;
    // Make the boundary surfaces from the 8 vertices.
    south  = new CoonsPatch(p0, p1, p5, p4);
    bottom = new CoonsPatch(p0, p1, p2, p3);
    west   = new CoonsPatch(p0, p3, p7, p4);
    east   = new CoonsPatch(p1, p2, p6, p5);
    north  = new CoonsPatch(p3, p2, p6, p7);
    top    = new CoonsPatch(p4, p5, p6, p7);
    set_and_check_corners();
    // cout << "Finished volume construction." << endl;
}


MeshVolume::MeshVolume( const vector<Vector3*> _p, int _ni, int _nj, int _nk,
			const string label, double r0, double r1, double s0, double s1, 
			double t0, double t1 )
    : ParametricVolume(label, r0, r1, s0, s1, t0, t1),
      p(vector<Vector3>()), ni(_ni), nj(_nj), nk(_nk)
{
    cout << "MeshVolume constructor: build from vertices." << endl;
    if ( (size_t)(ni*nj*nk) != _p.size() ) {
	cout << "MeshVolume Warning: ni=" << ni << " nj=" << nj 
	     << " p.size()=" << _p.size() << endl;
    }
    // Copy vertices from the supplied list.
    for ( size_t i = 0; i < _p.size(); ++i ) p.push_back(*(_p[i]));
    set_surfaces();
    set_and_check_corners();
    cout << "Finished MeshVolume construction." << endl;
}
MeshVolume::MeshVolume( const string file_name, int vtk_format,
			const string label, double r0, double r1, double s0, double s1, 
			double t0, double t1 )
    : ParametricVolume(label, r0, r1, s0, s1, t0, t1),
      p(vector<Vector3>()), ni(0), nj(0), nk(0)
{
    cout << "MeshVolume constructor: read from file " << file_name << endl;
    FILE* fp = fopen(file_name.c_str(), "r");
    if ( fp == NULL ) {
	cout << "MeshVolume constructor: could not open file " << file_name << endl;
    }

    char buffer[256];
    if ( vtk_format == 1 ) {
	// Discard the first 5 lines of the VTK-format file.
	// expect  # vtk DataFile Version 2.0
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cerr << "MeshPatch constructor: fgets failure." << endl;
	}
	// expect  some title string
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cerr << "MeshPatch constructor: fgets failure." << endl;
	}
	// expect  ASCII
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cerr << "MeshPatch constructor: fgets failure." << endl;
	}
	// expect  blank line
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cerr << "MeshPatch constructor: fgets failure." << endl;
	}
	// expect  DATASET STRUCTURED_GRID
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cerr << "MeshPatch constructor: fgets failure." << endl;
	}
	// Now we get to some interesting data.
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cerr << "MeshPatch constructor: fgets failure." << endl;
	}
	sscanf(buffer, "DIMENSIONS %d %d %d", &ni, &nj, &nk);
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cerr << "MeshPatch constructor: fgets failure." << endl;
	}
	int total_points;
	sscanf(buffer, "POINTS %d float\n", &total_points);
	if ( total_points != ni * nj * nk ) {
	    cout << "MeshVolume constructor: wrong number of total points=" 
		 << total_points << " ni*nj*nk=" << ni*nj*nk << endl;
	}
    } else {
	// Assume TECPLOT format.
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cerr << "MeshPatch constructor: fgets failure." << endl;
	}
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cerr << "MeshPatch constructor: fgets failure." << endl;
	}
	// Now, we get to interesting data.
	if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
	    cerr << "MeshPatch constructor: fgets failure." << endl;
	}
	sscanf(buffer, "ZONE I=%d, J=%d, K=%d, F=POINT", &ni, &nj, &nk);
    }
    // Read the vertex coordinates.
    // Note the order of the loops: i index runs fastest.
    double x, y, z;
    for ( int k = 0; k < nk; ++k ) {
	for ( int j = 0; j < nj; ++j ) {
	    for ( int i = 0; i < ni; ++i) {
		if ( fgets(buffer, sizeof(buffer), fp) == NULL ) {
		    cerr << "MeshPatch constructor: fgets failure." << endl;
		}
		sscanf(buffer, "%lf %lf %lf", &x, &y, &z );
		p.push_back( Vector3(x, y, z) );
	    }
	}
    }
    fclose(fp);
    set_surfaces();
    set_and_check_corners();
    cout << "Finished MeshVolume construction." << endl;
}
MeshVolume::~MeshVolume()
{
    // base class deletes the surfaces
}
int MeshVolume::set_surfaces()
{
    int i, j, k;
    cout << "MeshVolume::set_surfaces()" << endl;
    vector<Vector3*> south_plist = vector<Vector3*>();
    vector<Vector3*> north_plist = vector<Vector3*>();
    for ( k = 0; k < nk; ++k ) {
	for ( i = 0; i < ni; ++i ) {
	    j = 0;
	    south_plist.push_back(&p[k*ni*nj + j*ni + i]);
	    j = nj-1;
	    north_plist.push_back(&p[k*ni*nj + j*ni + i]);
	}
    }
    south = new MeshPatch(south_plist, ni, nk);
    north = new MeshPatch(north_plist, ni, nk);
    vector<Vector3*> west_plist = vector<Vector3*>();
    vector<Vector3*> east_plist = vector<Vector3*>();
    for ( k = 0; k < nk; ++k ) {
	for ( j = 0; j < nj; ++j ) {
	    i = 0;
	    west_plist.push_back(&p[k*ni*nj + j*ni + i]);
	    i = ni-1;
	    east_plist.push_back(&p[k*ni*nj + j*ni + i]);
	}
    }
    west = new MeshPatch(west_plist, nj, nk);
    east = new MeshPatch(east_plist, nj, nk);
    vector<Vector3*> bottom_plist = vector<Vector3*>();
    vector<Vector3*> top_plist = vector<Vector3*>();
    for ( j = 0; j < nj; ++j ) {
	for ( i = 0; i < ni; ++i ) {
	    k = 0;
	    bottom_plist.push_back(&p[k*ni*nj + j*ni + i]);
	    k = nk-1;
	    top_plist.push_back(&p[k*ni*nj + j*ni + i]);
	}
    }
    bottom = new MeshPatch(bottom_plist, ni, nj);
    top = new MeshPatch(top_plist, ni, nj);
    return 0;
}
Vector3 MeshVolume::eval( double r, double s, double t ) const 
{ 
    // Transform to parameter range of subsection of the volume.
    r = r0 + (r1 - r0) * r;
    s = s0 + (s1 - s0) * s;
    t = t0 + (t1 - t0) * t;

    // Interpolate within the background mesh.
    // This involves finding the relevant cell and
    // using a bilinear interpolation within that cell.
    double dr = 1.0 / (ni-1);
    double ds = 1.0 / (nj-1);
    double dt = 1.0 / (nk-1);
    int i = (int) (r / dr); i = (i < 0) ? 0 : i; i = (i > ni-2) ? ni-2 : i;
    int j = (int) (s / ds); j = (j < 0) ? 0 : j; j = (j > nj-2) ? nj-2 : j;
    int k = (int) (t / dt); k = (k < 0) ? 0 : k; k = (k > nk-2) ? nk-2 : k;
    // Identify the four corners of the cell.
    int indx000 = k*nj*ni + j*ni + i;
    int indx100 = k*nj*ni + j*ni + i+1;
    int indx110 = k*nj*ni + (j+1)*ni + i+1;
    int indx010 = k*nj*ni + (j+1)*ni + i;
    int indx001 = (k+1)*nj*ni + j*ni + i;
    int indx101 = (k+1)*nj*ni + j*ni + i+1;
    int indx111 = (k+1)*nj*ni + (j+1)*ni + i+1;
    int indx011 = (k+1)*nj*ni + (j+1)*ni + i;
    // Parametric coordinate within the coarse cell.
    double local_r = (r - dr * i) / dr;
    double local_s = (s - ds * j) / ds;
    double local_t = (t - dt * k) / dt;
    // Weight the corner points with co-volumes.
    Vector3 point = (1.0 - local_r) * (1.0 - local_s) * (1.0 - local_t) * p[indx000] +
	(1.0 - local_r) * local_s * (1.0 - local_t) * p[indx010] +
	local_r * (1.0 - local_s) * (1.0 - local_t) * p[indx100] +
	local_r * local_s * (1.0 - local_t) * p[indx110] +
	(1.0 - local_r) * (1.0 - local_s) * local_t * p[indx001] +
	(1.0 - local_r) * local_s * local_t * p[indx011] +
	local_r * (1.0 - local_s) * local_t * p[indx101] +
	local_r * local_s * local_t * p[indx111];
    return point;
}
string MeshVolume::str() const
{
    ostringstream ost;
    size_t i;
    ost << "MeshVolume(";
    ost << "p=[";
    if ( p.size() <= 10 ) {
	// Print full complement of points.
	for ( i = 0; i < p.size()-1; ++i ) ost << p[i] << ", ";
    } else {
	// There are too many points to print.
	for ( i = 0; i < 5; ++i ) ost << p[i] << ", ";
	ost << " ...<snipped>... , ";
	for ( i = p.size()-4; i < p.size()-1; ++i ) ost << p[i] << ", ";
    }
    ost << p[p.size()-1] << "], ";
    ost << "ni=" << ni << ", nj=" << nj << ", nk=" << nk << ", ";
    ost << "\"" << label << "\"" << ", ";
    ost << r0 << ", " << r1 << ", " 
	<< s0 << ", " << s1 << ", " 
	<< t0 << ", " << t1 << ")";
    return ost.str();
}
MeshVolume* MeshVolume::translate( const Vector3 &v )
{
    for ( size_t i = 0; i < p.size(); ++i ) p[i] += v;
    south->translate(v);
    bottom->translate(v);
    west->translate(v);
    east->translate(v);
    north->translate(v);
    top->translate(v);
    set_and_check_corners();
    return this;
}
MeshVolume* MeshVolume::translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    return this;
}
MeshVolume* MeshVolume::mirror_image( const Vector3 &point, const Vector3 &normal )
{
    for ( size_t i = 0; i < p.size(); ++i ) p[i].mirror_image(point, normal);
    south->mirror_image(point, normal);
    bottom->mirror_image(point, normal);
    west->mirror_image(point, normal);
    east->mirror_image(point, normal);
    north->mirror_image(point, normal);
    top->mirror_image(point, normal);
    set_and_check_corners();
    return this;
}
MeshVolume* MeshVolume::rotate_about_zaxis( double dtheta )
{
    for ( size_t i = 0; i < p.size(); ++i ) p[i].rotate_about_zaxis(dtheta);
    south->rotate_about_zaxis(dtheta);
    bottom->rotate_about_zaxis(dtheta);
    west->rotate_about_zaxis(dtheta);
    east->rotate_about_zaxis(dtheta);
    north->rotate_about_zaxis(dtheta);
    top->rotate_about_zaxis(dtheta);
    set_and_check_corners();
    return this;
}
