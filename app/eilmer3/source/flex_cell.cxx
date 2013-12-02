#include <string>
#include <stdio.h>
#include "../../../lib/util/source/useful.h"
#include "flex_cell.hh"
#include "cell.hh"
#include "kernel.hh"

// ----------------------------------------------------------
/*
 * for all intents and purposes, flex_cell behaves the same way as a cell.
 * probably all cell functions can be called using a flex_cell, 
 * since it is a derived struct.
 */

int read_solution_for_flex_cell(flex_cell_center* fc, FILE *infile)
{
#if 0
#   define NCHAR 4000
    char line[NCHAR];
    if (fgets( line, NCHAR, infile ) == NULL) {
	printf( "read_solution_for_flex_cell(): Empty flow field.\n" );
	exit( BAD_INPUT_ERROR );
    }
    return fc->scan_values_from_string(line);
#   undef NCHAR
#endif
    return 0;
}

int write_solution_for_flex_cell(flex_cell_center* fc, FILE *outfile)
{
#if 0
    string str = fc->write_values_to_string();
    fputs( str.c_str(), outfile );
    fputc( '\n', outfile );
#endif
    return SUCCESS;
}

int set_array_sizes_for_flex_cell(flex_cell_center &fc, size_t nsp, size_t nvib)
{
#if 0
    int nb = fc.set_array_sizes(nsp, nvib);

    fc.mf.resize(nsp);
    fc.mf_old.resize(nsp);

    nb += 3 * nsp * sizeof(double);

    for ( size_t i = 0; i < N_LEVEL; ++i ) {
	fc.DmfDt[i].resize(nsp);
	nb += nsp * sizeof(double);
    }
    fc.menergies.resize(nvib);
    fc.menergies_old.resize(nvib);
    nb += 3 * nvib * sizeof(double);
    for ( size_t i = 0; i < N_LEVEL; ++i ) {
	fc.DmenergiesDt[i].resize(nvib);
	nb += nvib * sizeof(double);
    }
    
    // set the flex-cell status to shadowed
    // so that when data is copied to a grid cell, this status is retained
    fc.status = SHADOWED_CELL;
    return nb;
#endif
    return 0;
}

int record_conserved_for_flex_cell(flex_cell_center* fc)
{
#if 0
    fc->record_conserved();

    fc->mass_old = fc->mass;
    fc->mv_old.x = fc->mv.x;
    fc->mv_old.y = fc->mv.y;
    fc->mv_old.z = fc->mv.z;
    fc->mE_old = fc->mE;
    for ( size_t isp = 0; isp < fc->mf.size(); ++isp ) {
	fc->mf_old[isp] = fc->mf[isp];
    }
    for ( size_t imode = 0; imode < fc->menergies.size(); ++imode ) {
	fc->menergies_old[imode] = fc->menergies[imode];
    }
#endif
    return SUCCESS;
}

int restore_conserved_for_flex_cell(flex_cell_center* fc) 
{
#if 0
    fc->mass = fc->mass_old;
    fc->mv.x = fc->mv_old.x;
    fc->mv.y = fc->mv_old.y;
    fc->mv.z = fc->mv_old.z;
    fc->mE = fc->mE_old;
    for ( size_t isp = 0; isp < fc->mf.size(); ++isp ) {
	fc->mf[isp] = fc->mf_old[isp];
    }
    for ( size_t imode = 0; imode < fc->menergies.size(); ++imode ) {
	fc->menergies[imode] = fc->menergies_old[imode];
    }
#endif
    return SUCCESS;
}

int encode_intensive_to_extensive(flex_cell_center* fc) 
{
#if 0
    calculate_flex_cell_volume(fc);
    double vol = fc->volume;

    // encode conserved
    fc->encode_conserved();

    // calculate conserved quantities
    // total mass
    fc->mass = vol * fc->U.mass;
    // total momentum
    fc->mv.x = vol * fc->U.momentum.x;
    fc->mv.y = vol * fc->U.momentum.y;
    fc->mv.z = vol * fc->U.momentum.z;
    // total energy
    fc->mE = vol*fc->U.total_energy;
    for ( size_t isp = 0; isp < fc->mf.size(); ++isp) {
	fc->mf[isp] = vol * fc->U.massf[isp];
    }
    for ( size_t isp = 0; isp < fc->U.energies.size(); ++isp) {
	fc->menergies[isp] = vol * fc->U.energies[isp];
    }
#endif
    return SUCCESS;
}

int decode_extensive_to_intensive(flex_cell_center* fc) 
{
#if 0
    calculate_flex_cell_volume(fc);
    double vinv = 1.0/fc->volume;

    // calculate conserved quantities
    // intensive mass
    fc->U.mass = fc->mass * vinv;
    // intensive momentum
    fc->U.momentum.x = fc->mv.x * vinv;
    fc->U.momentum.y = fc->mv.y * vinv;
    fc->U.momentum.z = fc->mv.z * vinv;
    // intensive total energy
    fc->U.total_energy = fc->mE*vinv;

    for (size_t isp = 0; isp < fc->mf.size(); ++isp) {
	fc->U.massf[isp] = fc->mf[isp] * vinv;
    }
    for (size_t isp = 0; isp < fc->menergies.size(); ++isp) {
	fc->U.energies[isp] = fc->menergies[isp] * vinv;
    }
    
    // decode conserved
    fc->decode_conserved();
#endif
    return SUCCESS;
}

int copy_shadowed_cell_data_to_flex_cell(FV_Cell* src, flex_cell_center* dest, int type_of_copy) 
{
#if 0
    if (type_of_copy != COPY_FLOW_STATE) {
	cout << "copy_shadowed_cell_data_to_flex_cell called" << endl;
	cout << "with type_of_copy = " << type_of_copy << endl;
	cout << "not yet implemented." << endl;
	exit(NOT_IMPLEMENTED_ERROR);
    }
    dest->copy_values_from(src, type_of_copy);
    encode_intensive_to_extensive(dest);
#endif
    return SUCCESS;
}

int copy_flex_cell_data_to_shadowed_cell(flex_cell_center* src, FV_Cell* dest, int type_of_copy) 
{
#if 0
    if (type_of_copy != COPY_FLOW_STATE) {
	cout << "copy_flex_cell_data_to_shadowed_cell called" << endl;
	cout << "with type_of_copy = " << type_of_copy << endl;
	cout << "not yet implemented." << endl;
	exit(NOT_IMPLEMENTED_ERROR);
    }
    decode_extensive_to_intensive(src);
    dest->copy_values_from(src, type_of_copy);
    dest->encode_conserved();
#endif
    return SUCCESS;
}

int calculate_flex_cell_volume(flex_cell_center* fc)
{
#if 0
    global_data &G = *get_global_data_ptr();
    double xA, yA, xB, yB, xC, yC, xD, yD;
    double xyarea;

    // These are the corners.
    xA = fc->vtx[1]->pos.x;
    yA = fc->vtx[1]->pos.y;
    xB = fc->vtx[2]->pos.x;
    yB = fc->vtx[2]->pos.y;
    xC = fc->vtx[3]->pos.x;
    yC = fc->vtx[3]->pos.y;
    xD = fc->vtx[0]->pos.x;
    yD = fc->vtx[0]->pos.y;

    // Cell area in the (x,y)-plane.
    xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
		    (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
    // Cell Centroid.
    fc->pos.x = 1.0 / (xyarea * 6.0) * 
	   ((yB - yA) * (xA * xA + xA * xB + xB * xB) + 
	    (yC - yB) * (xB * xB + xB * xC + xC * xC) +
	    (yD - yC) * (xC * xC + xC * xD + xD * xD) + 
	    (yA - yD) * (xD * xD + xD * xA + xA * xA));
    fc->pos.y = -1.0 / (xyarea * 6.0) * 
           ((xB - xA) * (yA * yA + yA * yB + yB * yB) + 
	    (xC - xB) * (yB * yB + yB * yC + yC * yC) +
	    (xD - xC) * (yC * yC + yC * yD + yD * yD) + 
	    (xA - xD) * (yD * yD + yD * yA + yA * yA));
    
    fc->area = xyarea;

    if ( G.axisymmetric ) {
	fc->volume = xyarea * fc->pos.y; // volume per unit radian
    } else {
	fc->volume = xyarea * 1.0; // volume per unit depth
    }

    if (xyarea < 0.0) {
	cout << "Negative flex_cell volume = " << xyarea << endl;
	cout << "pos = " 
	     << "(" << fc->pos.x << ", " << fc->pos.y << ")"
	     << endl;
	cout << "vertices: \n"
	     << "(" << xD << ", " << yD << ")\n"
	     << "(" << xA << ", " << yA << ")\n"
	     << "(" << xB << ", " << yB << ")\n"
	     << "(" << xC << ", " << yC << ")"
	     << endl;
	return 1;
    }

    fc->pos = (fc->vtx[0]->pos + 
	       fc->vtx[1]->pos + 
	       fc->vtx[2]->pos + 
	       fc->vtx[3]->pos) / 4.0;
#endif
    return SUCCESS;
}

int calculate_flex_cell_length(flex_cell_center* fc)
{
#if 0
    fc->iLength = ( fc->vtx[1]->pos.x - fc->vtx[0]->pos.x );
    fc->jLength = ( fc->vtx[3]->pos.y - fc->vtx[0]->pos.y );
    
    if ((fc->jLength) < 0.0) {
	cout << "Negative flex_cell jLength = " << fc->jLength << endl;
	return 1;
    }	

    if ((fc->iLength) < 0.0) {
	cout << "Negative flex_cell iLength = " << fc->iLength << endl;
	return 1;
    }	
#endif
    return SUCCESS;
} // end of calculate_flex_cell_length

int add_shadowed_cell_and_flex_cell_data(FV_Cell* src, flex_cell_center* dest)
{
#if 0
//     cout << "dest->mass = " << dest->mass;
//     cout << ", src->mass = " << src->fs.gas.rho*src->volume << endl;
//     cout << "dest->mv.x = " << dest->mv.x;
//     cout << ", src->mv.x = " << src->rv.x*src->volume << endl;
//     cout << "dest->mE = " << dest->mE;
//     cout << ", src->mE = " << src->volume*(src->rE) << endl;

    // add flow states
    dest->mass += src->U.mass*src->volume;

    dest->mv.x += src->U.momentum.x*src->volume;
    dest->mv.y += src->U.momentum.y*src->volume;
    dest->mv.z += src->U.momentum.z*src->volume;
    dest->mE += src->volume*(src->U.total_energy);
    
    /* Species densities: mass of species isp per unit volume. */
    for (size_t isp = 0; isp < dest->mf.size(); ++isp) {
	dest->mf[isp] += src->U.massf[isp]*src->volume;
    }
    for (size_t isp = 0; isp < dest->menergies.size(); ++isp) {
	dest->menergies[isp] += src->U.energies[isp]*src->volume;
    }
#endif
    return SUCCESS;
}

int time_derivatives_for_flex_cell(flex_cell_center *fc, int ftl, size_t dimensions) 
{
#if 0
    if (dimensions != 2) {
	cout << "time_derivatives_for_flex_cell called" << endl;
	cout << "with dimensions = " << dimensions << endl;
	cout << "not yet implemented." << endl;
	exit(NOT_IMPLEMENTED_ERROR);
    }

    calculate_flex_cell_volume(fc);
    double vol = fc->volume;

    // masquerade as a normal cell
    fc->status = NORMAL_CELL;
    fc->time_derivatives(ftl, dimensions);
    fc->status = MASKED_CELL;
    
    // calculate extensive time derivatives
    fc->DmDt[ftl] = vol * fc->dUdt[ftl].mass;
    
    fc->DmvDt[ftl].x = vol * fc->dUdt[ftl].momentum.x;
    fc->DmvDt[ftl].y = vol * fc->dUdt[ftl].momentum.y;
    fc->DmvDt[ftl].z = 0.0;
    
    fc->DmEDt[ftl] = vol * fc->dUdt[ftl].total_energy;

    for ( size_t isp = 0; isp < fc->DmfDt[ftl].size(); ++isp ) {
	fc->DmfDt[ftl][isp] = vol * fc->dUdt[ftl].massf[isp];
    }
    for ( size_t imode = 0; imode < fc->DmenergiesDt[ftl].size(); ++imode ) {
	fc->DmenergiesDt[ftl][imode] = vol * fc->dUdt[ftl].energies[imode];
    }
#endif
    return SUCCESS;
}

int calculate_wall_flux(flex_cell_center *fc, double u)
{
#if 0
/**
\brief
Here we use flex_cell values of pressure and piston velocity to set the flux
vector. These are dirichlet boundary conditions for a moving wall. 

**/ 
    FV_Interface* ifp;

    // if x-direction is negative, this is a west flex_cell, point to the east interface
    if (fc->n->x < 0) ifp = fc->iface[EAST];
    // if x-direction is positive, this is an east flex_cell, point to the west interface
    else ifp = fc->iface[WEST];

    // This is the flux vector for a moving wall boundary condition
    ifp->F.mass = 0.0;
    ifp->F.momentum.x = fc->fs.gas.p;
    ifp->F.total_energy = fc->fs.gas.p*u;
#endif
    return SUCCESS;
}

double signal_frequency_for_flex_cell(flex_cell_center *fc, size_t dimensions)
{
#if 0
    if (dimensions != 2) {
	cout << "signal_frequency_for_flex_cell called" << endl;
	cout << "with dimensions = " << dimensions << endl;
	cout << "not yet implemented." << endl;
	exit(NOT_IMPLEMENTED_ERROR);
    }
    return fc->signal_frequency(dimensions);
#endif
    return 0.0;
}

int update_for_flex_cell(flex_cell_center *fc, double u, double dt)
{
#if 0
    double A_av = 0.5*(fc->iface[EAST]->area + fc->iface[WEST]->area);

    fc->mass = fc->mass_old + dt * fc->DmDt[0];
    fc->mv.x = fc->mv_old.x + dt * fc->DmvDt[0].x;
    fc->mv.y = fc->mv_old.y + dt * fc->DmvDt[0].y;
    fc->mv.z = fc->mv_old.z + dt * fc->DmvDt[0].z;
    fc->mE = fc->mE_old + dt * fc->DmEDt[0];

    for ( size_t isp = 0; isp < fc->mf.size(); ++isp ) {
	fc->mf[isp] = fc->mf_old[isp] + dt * fc->DmfDt[0][isp];
    }
    for ( size_t isp = 0; isp < fc->menergies.size(); ++isp ) {
	fc->menergies[isp] = fc->menergies_old[isp] + dt * fc->DmenergiesDt[0][isp];
    }
    fc->iLength = fc->iLength - fc->n->x*u*dt;
    fc->volume = fc->volume - fc->n->x*u*dt*A_av;
#endif
    return SUCCESS;
}

int inviscid_source_vector_for_flex_cell(flex_cell_center *fc)
{
#if 0
    return fc->inviscid_source_vector();
#endif
    return 0;
}

int update_flex_cell_geometry_based_on_vertex_movement(flex_cell_center *src) {
#if 0
    global_data &G = *get_global_data_ptr();
    // Vertex numeral convention.
    //
    // 3,C .--. 2,B
    //     |  |
    // 0,D .--. 1,A

    double xA, xB, yA, yB, xC, yC, xD, yD, LAB;
    xA = src->vtx[1]->pos.x;
    yA = src->vtx[1]->pos.y;
    xB = src->vtx[2]->pos.x; 
    yB = src->vtx[2]->pos.y;
    xC = src->vtx[3]->pos.x;
    yC = src->vtx[3]->pos.y;
    xD = src->vtx[0]->pos.x;
    yD = src->vtx[0]->pos.y;

    // -----------------------------------------------------------------
    // EAST FACE
    LAB = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA));
    if (LAB < 1.0e-9) {
	cout << "Zero ifi length " << LAB << endl;
    }
    // Direction cosines for the unit normal.
    src->iface[EAST]->n.x = (yB - yA) / LAB;
    src->iface[EAST]->n.y = -(xB - xA) / LAB;
    src->iface[EAST]->n.z = 0.0;  // 2D plane
    // Length in the XY-plane.
    src->iface[EAST]->length = LAB;
    // Mid-point and area.
    src->iface[EAST]->Ybar = 0.5 * (yA + yB);
    if ( G.axisymmetric ) {
	// Interface area per radian.
	src->iface[EAST]->area = LAB * src->iface[EAST]->Ybar;
    } else {
	// Assume unit depth in the Z-direction.
	src->iface[EAST]->area = LAB;
    }
    src->iface[EAST]->pos =  (src->vtx[1]->pos + src->vtx[2]->pos)/2.0;
       
    // -----------------------------------------------------------------
    // WEST FACE
    LAB = sqrt((xC - xD) * (xC - xD) + (yC - yD) * (yC - yD));
    if (LAB < 1.0e-9) {
	cout << "Zero ifi length " << LAB << endl;
    }
    // Direction cosines for the unit normal.
    src->iface[WEST]->n.x = (yC - yD) / LAB;
    src->iface[WEST]->n.y = -(xC - xD) / LAB;
    src->iface[WEST]->n.z = 0.0;  // 2D plane
    // Length in the XY-plane.
    src->iface[WEST]->length = LAB;
    // Mid-point and area.
    src->iface[WEST]->Ybar = 0.5 * (yC + yD);
    if ( G.axisymmetric ) {
	// Interface area per radian.
	src->iface[WEST]->area = LAB * src->iface[WEST]->Ybar;
    } else {
	// Assume unit depth in the Z-direction.
	src->iface[WEST]->area = LAB;
    }
    src->iface[WEST]->pos =  (src->vtx[0]->pos + src->vtx[3]->pos)/2.0;
       
    // -----------------------------------------------------------------
    // NORTH FACE
    LAB = sqrt((xB - xC) * (xB - xC) + (yB - yC) * (yB - yC));
    if (LAB < 1.0e-9) {
	cout << "Zero ifi length " << LAB << endl;
    }
    // Direction cosines for the unit normal.
    src->iface[NORTH]->n.x = (yB - yC) / LAB;
    src->iface[NORTH]->n.y = -(xB - xC) / LAB;
    src->iface[NORTH]->n.z = 0.0;  // 2D plane
    // Length in the XY-plane.
    src->iface[NORTH]->length = LAB;
    // Mid-point and area.
    src->iface[NORTH]->Ybar = 0.5 * (yB + yC);
    if ( G.axisymmetric ) {
	// Interface area per radian.
	src->iface[NORTH]->area = LAB * src->iface[NORTH]->Ybar;
    } else {
	// Assume unit depth in the Z-direction.
	src->iface[NORTH]->area = LAB;
    }
    src->iface[NORTH]->pos =  (src->vtx[2]->pos + src->vtx[3]->pos)/2.0;

    // -----------------------------------------------------------------
    // SOUTH FACE
    LAB = sqrt((xA - xD) * (xA - xD) + (yA - yD) * (yA - yD));
    if (LAB < 1.0e-9) {
	cout << "Zero ifi length " << LAB << endl;
    }
    // Direction cosines for the unit normal.
    src->iface[SOUTH]->n.x = (yA - yD) / LAB;
    src->iface[SOUTH]->n.y = -(xA - xD) / LAB;
    src->iface[SOUTH]->n.z = 0.0;  // 2D plane
    // Length in the XY-plane.
    src->iface[SOUTH]->length = LAB;
    // Mid-point and area.
    src->iface[SOUTH]->Ybar = 0.5 * (yA + yD);
    if ( G.axisymmetric ) {
	// Interface area per radian.
	src->iface[SOUTH]->area = LAB * src->iface[SOUTH]->Ybar;
    } else {
	// Assume unit depth in the Z-direction.
	src->iface[SOUTH]->area = LAB;
    }
    src->iface[SOUTH]->pos = (src->vtx[0]->pos + src->vtx[1]->pos)/2.0;

    // -----------------------------------------------------
    // update some cell geometry

    src->pos = 0.25 * (src->iface[WEST]->pos + src->iface[EAST]->pos + 
		       src->iface[SOUTH]->pos + src->iface[NORTH]->pos);

    calculate_flex_cell_length(src);
    calculate_flex_cell_volume(src);
#endif
    return SUCCESS;
}

int print_data_for_flex_cell(flex_cell_center *fc)
{
#if 0
    cout << "--------------------------------"<< endl;
    cout << "--- begin data for flex-cell ---" << endl;

    cout << "pos = (" << fc->pos.x << ", "
	 << fc->pos.y << ", "
	 << fc->pos.z << ")" << endl;
    cout << "base_qdot = " << fc->base_qdot << endl;
    cout << endl;
    cout << "rho = " << fc->fs.gas.rho << "\n"
         << "p   = " << fc->fs.gas.p << "\n"
         << "T[0]= " << fc->fs.gas.T[0] << "\n"
	 << "u   = " << fc->fs.vel.x << "\n"
	 << "v   = " << fc->fs.vel.y << "\n"
	 << "w   = " << fc->fs.vel.z << "\n"
	 << "e[0]= " << fc->fs.gas.e[0] << endl;
    cout << endl;
    double f_sum = 0.0;
    size_t nsp = fc->fs.gas.massf.size();
    for ( size_t isp = 0; isp < nsp; ++isp) {
	cout << "f[" << isp << "]  = " << fc->fs.gas.massf[isp] << endl;
	f_sum += fc->fs.gas.massf[isp];
    }
    cout << "f_sum = " << f_sum << endl;

    cout << endl;
    cout << "Conserved quantities:\n"
	 << "rv.x = " << fc->U.momentum.x << "\n"
	 << "rv.y = " << fc->U.momentum.y << "\n"
	 << "rv.z = " << fc->U.momentum.z << "\n"
	 << "rE   = " << fc->U.total_energy << endl;
    cout << endl;
    cout << "flex_cell status = " << fc->status << endl;
    cout << "---- end data for flex-cell ----" << endl;
    cout << "--------------------------------"<< endl;
#endif
    return SUCCESS;
}
