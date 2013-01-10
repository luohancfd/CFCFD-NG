/** file cell.cxx
 * \ingroup eilmer3
 * \brief Classes and functions to deal with cell-related data.
 *
 * \author PJ
 * \version 05-Aug-04 reconstituted from cns_bc.c, cns_tstp.c, etc.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <sstream>

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/geometry2/source/geom.hh"
#include "one_d_interp_scalar.hh"
#include "cell.hh"
#include "kernel.hh"
#include "flux_calc.hh"
#include "diffusion.hh"
#include "bgk.hh"

#define VISCOUS_TIME_LIMIT_MODEL 0 // (0) original Swanson model, (1) Ramshaw model

/*----------------------------------------------------------------*/

static string noname = "none";
static string face_name[] = { "north", "east", "south", "west", "top", "bottom" };

/// \brief Translates a face name to a numeric index.
int get_face_index(const string name)
{
    std::string  name_copy = name;
    for ( size_t i = 0; i < name_copy.length(); ++i )
	name_copy[i] = tolower(name_copy[i]);
    if ( name_copy == face_name[NORTH] ) return NORTH;
    if ( name_copy == face_name[SOUTH] ) return SOUTH;
    if ( name_copy == face_name[EAST] ) return EAST;
    if ( name_copy == face_name[WEST] ) return WEST;
    if ( name_copy == face_name[TOP] ) return TOP;
    if ( name_copy == face_name[BOTTOM] ) return BOTTOM;
    return -1; // not a valid face index
}

/// \brief Translates an index to a face name.
std::string get_face_name(int index)
{
    if ( index >= 0 && index <= 5 )
	return face_name[index];
    else
	return noname;
}

//---------------------------------------------------------------------------

FlowState::FlowState(Gas_model *gm)
{
    gas = new Gas_data(gm);
    vel.x = 0.0; vel.y = 0.0; vel.z = 0.0;
    B.x = 0.0; B.y = 0.0; B.z = 0.0;
    shock_vel.x = 0.0; shock_vel.y = 0.0; shock_vel.z = 0.0;
    S = 0;
    tke = 0.0;
    omega = 0.0;
    mu_t = 0.0;
    k_t = 0.0;
    G.resize(get_velocity_buckets(), 0.0);
    H.resize(get_velocity_buckets(), 0.0);
}

FlowState::~FlowState()
{
    delete gas;
    G.resize(0, 0.0);
    H.resize(0, 0.0);
}

int FlowState::print()
{
    printf("----------- Data for a flow state... ------------\n");
    gas->print_values();
    printf("v.x= %e, v.y= %e, v.z= %e \n", vel.x, vel.y, vel.z);
    if ( get_mhd_flag() == 1 ) {
	printf("B.x= %e, B.y=%e, B.z=%e \n", B.x, B.y, B.z);
    }
    if ( get_velocity_buckets() > 0) {
	printf("Velocity Distribution Partial Densities:\n");
	for ( size_t ipd = 0; ipd < G.size(); ++ipd ) {
	    printf("i=%d: G=%e, H=%e\n",(int) ipd, G[ipd], H[ipd]);
	}
    }
    printf("tke= %e, omega=%e\n", tke, omega);
    printf("S= %d, mu_t=%e, k_t=%e\n", S, mu_t, k_t);
    return SUCCESS;
}

int FlowState::copy_values_from(FlowState &src)
{
    gas->copy_values_from(*(src.gas));
    vel.x = src.vel.x; vel.y = src.vel.y; vel.z = src.vel.z;
    B.x = src.B.x; B.y = src.B.y; B.z = src.B.z;
    S = src.S;
    tke = src.tke;
    omega = src.omega;
    mu_t = src.mu_t;
    k_t = src.k_t;
    for ( size_t ipd = 0; ipd < src.G.size(); ++ipd ) {
	G[ipd] = src.G[ipd];
	H[ipd] = src.H[ipd];
    }
    return SUCCESS;
}

int FlowState::copy_values_from(CFlowCondition &src)
{
    gas->copy_values_from(*(src.gas));
    vel.x = src.u; vel.y = src.v; vel.z = src.w;
    if ( get_mhd_flag() == 1 ) {
	B.x = src.Bx; B.y = src.By; B.z = src.Bz;
    }
    S = src.S;
    tke = src.tke;
    omega = src.omega;
    mu_t = src.mu_t;
    k_t = src.k_t;
    return SUCCESS;
}

int FlowState::average_values_from(FlowState &src0, FlowState &src1, bool with_diff_coeff)
{
    gas->average_values_from(*(src0.gas), 0.5, *(src1.gas), 0.5, with_diff_coeff);
    vel.x = 0.5 * (src0.vel.x + src1.vel.x);
    vel.y = 0.5 * (src0.vel.y + src1.vel.y);
    vel.z = 0.5 * (src0.vel.z + src1.vel.z);
    if ( get_mhd_flag() == 1 ) {
	B.x = 0.5 * (src0.B.x + src1.B.x);
	B.y = 0.5 * (src0.B.y + src1.B.y);
	B.z = 0.5 * (src0.B.z + src1.B.z);
    }
    S = src0.S || src1.S; // err on detecting a shock
    tke = 0.5 * (src0.tke + src1.tke);
    omega = 0.5 * (src0.omega + src1.omega);
    mu_t = 0.5 * (src0.mu_t + src1.mu_t);
    k_t = 0.5 * (src0.k_t + src1.k_t);
    for ( size_t ipd = 0; ipd < src0.G.size(); ++ipd ) {
	G[ipd] = 0.5 * (src0.G[ipd] + src1.G[ipd]);
	H[ipd] = 0.5 * (src0.H[ipd] + src1.H[ipd]);
    }
    return SUCCESS;
}

/// \brief Copy the FlowState data into a linear data buffer.
/// \param buf : pointer to the current element somewhere in buffer
/// \returns a pointer to the next available location in the data buffer.
double * FlowState::copy_values_to_buffer(double *buf)
{
    buf = gas->copy_values_to_buffer(buf);
    *buf++ = vel.x;
    *buf++ = vel.y;
    *buf++ = vel.z;
    if ( get_mhd_flag() == 1 ) {
	*buf++ = B.x;
	*buf++ = B.y;
	*buf++ = B.z;
    }
    *buf++ = (double) S;
    *buf++ = tke;
    *buf++ = omega;
    *buf++ = mu_t;
    *buf++ = k_t;
    for ( size_t ipd = 0; ipd < G.size(); ++ipd ) {
	*buf++ = G[ipd];
	*buf++ = H[ipd];
    }
	
    return buf;
}

/// \brief Copy the data from a linear data buffer into the FlowState structure.
/// \param buf : pointer to the current element somewhere in buffer
/// \returns a pointer to the next available location in the data buffer.
double * FlowState::copy_values_from_buffer(double *buf)
{
    buf = gas->copy_values_from_buffer(buf);
    vel.x = *buf++;
    vel.y = *buf++;
    vel.z = *buf++;
    if ( get_mhd_flag() == 1 ) {
	B.x = *buf++;
	B.y = *buf++;
	B.z = *buf++;
    }
    S = (int)(*buf++);
    tke = *buf++;
    omega = *buf++;
    mu_t = *buf++;
    k_t = *buf++;
    for ( size_t ipd = 0; ipd < G.size(); ++ipd ) {
	G[ipd] = *buf++;
	H[ipd] = *buf++;
    }
    return buf;
}

/// \brief From the macroscopic quantities set the values of 
/// G and H to those for equilibrium.
int FlowState::BGK_equilibrium(void)
{

    //DARYL - FIX-ME: check the definition of temperature, heat flux coeff, what is "status" for

    Vector3 gh, uvw;
    int status;

    Gas_model *gmodel = get_gas_model_ptr();
    
    double R = gmodel->R(*gas, status);
    double Cp = gmodel->Cp(*gas, status);
    double Pr = gmodel->Prandtl(gas->mu, Cp, gas->k[0]);

    for (size_t iq = 0; iq < get_velocity_buckets(); ++iq) {
	uvw = get_vcoord(iq);	
	gh = Shakhov(gas->rho, vel.x, vel.y, gas->T[0], 0.0, 0.0, R, Pr, uvw.x, uvw.y);
	G[iq] = gh.x;
	H[iq] = gh.y;
    }
    return SUCCESS;
}

//----------------------------------------------------------------------------

ConservedQuantities::ConservedQuantities(Gas_model *gm)
{
    mass = 0.0;
    momentum.x = 0.0; momentum.y = 0.0; momentum.z = 0.0;
    B.x = 0.0; B.y = 0.0; B.z = 0.0;
    total_energy = 0.0;
    massf.resize(gm->get_number_of_species(), 0.0);
    energies.resize(gm->get_number_of_modes(), 0.0);
    tke = 0.0;
    omega = 0.0;
    G.resize(get_velocity_buckets(), 0.0);
    H.resize(get_velocity_buckets(), 0.0);
}

ConservedQuantities::~ConservedQuantities()
{
    massf.clear();
    energies.clear();
    G.clear();
    H.clear();
}

int ConservedQuantities::print()
{
    cout << "mass= " << mass << endl;
    cout << "momentum.x= " << momentum.x << " .y= " << momentum.y
	 << " .z= " << momentum.z << endl;
    if ( get_mhd_flag() == 1 ) {
	cout << "B.x= " << B.x << " .y= " << B.y << " .z= " << B.z << endl;
    }
    cout << "total_energy= " << total_energy << endl;
    cout << "massf= ";
    for ( size_t isp = 0; isp < massf.size(); ++isp )
	cout << isp << ":" << massf[isp] << " ";
    cout << endl;
    cout << "energies= ";
    for ( size_t imode = 0; imode < energies.size(); ++imode )
	cout << imode << ":" << energies[imode] << " ";
    cout << endl;
    cout << "tke= " << tke << " omega=" << omega << endl;
    if ( get_velocity_buckets() > 0) {
	printf("Velocity Distribution Partial Densities:\n");
	for ( size_t ipd = 0; ipd < G.size(); ++ipd ) {
	    cout << ipd << ": (" << G[ipd] << ", " << H[ipd] << ")";
	}
	cout << endl;
    }
    return SUCCESS;
}

int ConservedQuantities::copy_values_from(ConservedQuantities &src)
{
    mass = src.mass;
    momentum.x = src.momentum.x;
    momentum.y = src.momentum.y;
    momentum.z = src.momentum.z;
    B.x = src.B.x; B.y = src.B.y; B.z = src.B.z;
    total_energy = src.total_energy;
    for ( size_t isp = 0; isp < src.massf.size(); ++isp )
	massf[isp] = src.massf[isp];
    for ( size_t imode = 0; imode < src.energies.size(); ++imode )
	energies[imode] = src.energies[imode];
    tke = src.tke;
    omega = src.omega;
    for (size_t ipd = 0; ipd < src.G.size(); ++ipd) {
	G[ipd] = src.G[ipd];
	H[ipd] = src.H[ipd];
    }
    return SUCCESS;
}

int ConservedQuantities::clear_values()
{
    mass = 0.0;
    momentum.x = 0.0;
    momentum.y = 0.0;
    momentum.z = 0.0;
    B.x = 0.0; B.y = 0.0; B.z = 0.0;
    total_energy = 0.0;
    for ( size_t isp = 0; isp < massf.size(); ++isp )
	massf[isp] = 0.0;
    for ( size_t imode = 0; imode < energies.size(); ++imode )
	energies[imode] = 0.0;
    tke = 0.0;
    omega = 0.0;
    for (size_t ipd = 0; ipd < G.size(); ++ipd) {
	G[ipd] = 0.0;
	H[ipd] = 0.0;
    }
    return SUCCESS;
}

//----------------------------------------------------------------------------

FV_Interface::FV_Interface(Gas_model *gm)
{
    id = 0;
    status = 0;
    pos.x = 0.0; pos.y = 0.0; pos.z = 0.0;
    vel.x = 0.0; vel.y = 0.0; vel.z = 0.0;
    Ybar = 0.0;
    length = 0.0;
    for ( size_t i = 0; i < NL; ++i ) {
	ar[i] = 0.0;	
    }
    area = 0.0;
    n.x = 0.0; n.y = 0.0; n.z = 0.0;
    t1.x = 0.0; t1.y = 0.0; t1.z = 0.0;
    t2.x = 0.0; t2.y = 0.0; t2.z = 0.0;
    fs = new FlowState(gm);
    F = new ConservedQuantities(gm);
#   if WITH_IMPLICIT == 1
    // Ojas, may be good to initialize point-implicit variables
    // double Lambda[6][2], R[6][6], R_inv[6][6], sl[6][2], sl_min[6][2];
    // double A[6][6], J[6][6], hl[6][2], gl[6][2];
#   endif
}

FV_Interface::~FV_Interface()
{
    delete fs;
    delete F;
}

int FV_Interface::print(int to_stdout)
{
    printf( "----------- Begin data for interface -----------\n");
    // printf( "id = %i\n", iface->id);
    printf("x=%e, y=%e, z=%e\n", pos.x, pos.y, pos.z);
    printf("area=%e, Ybar=%e, length=%e\n", area, Ybar, length);
    printf("n.x=%e, n.y=%e, n.z=%e\n", n.x, n.y, n.z);
    printf("t1.x=%e, t1.y=%e, t1.z=%e\n", t1.x, t1.y, t1.z);
    printf("t2.x=%e, t2.y=%e, t2.z=%e\n", t2.x, t2.y, t2.z);
    fs->print();
    printf("Fluxes of conserved quantities: \n");
    F->print();
    printf("status= %d \n", status);
    printf("----------- End data for interface -----------\n");
    return SUCCESS;
}

/// \brief Copies data between two interface structures
int FV_Interface::copy_values_from(FV_Interface &src, int type_of_copy)
{
    if ( type_of_copy == COPY_ALL_CELL_DATA ||
	 type_of_copy == COPY_FLOW_STATE ) {
        fs->copy_values_from(*(src.fs));
	F->copy_values_from(*(src.F));
    }
    if ( type_of_copy == COPY_ALL_CELL_DATA ||
	 type_of_copy == COPY_CELL_LENGTHS ) {
	pos.x = src.pos.x; pos.y = src.pos.y; pos.z = src.pos.z;
	vel.x = src.vel.x; vel.y = src.vel.y; vel.z = src.vel.z;
	n.x = src.n.x; n.y = src.n.y; n.z = src.n.z;
	t1.x = src.t1.x; t1.y = src.t1.y; t1.z = src.t1.z;
	t2.x = src.t2.x; t2.y = src.t2.y; t2.z = src.t2.z;
	Ybar = src.Ybar; length = src.length; area = src.area;
    }
    return SUCCESS;
}

//----------------------------------------------------------------------------

FV_Vertex::FV_Vertex(Gas_model *gm)
{
    int nsp = gm->get_number_of_species();
    int nmodes = gm->get_number_of_modes();
    id = 0;
    pos.x = 0.0; pos.y = 0.0; pos.z = 0.0;
    for ( size_t i = 0; i < NL; ++i ) {
	    position[i].x = 0.0; position[i].y = 0.0; position[i].z = 0.0;
    }
    vel.x = 0.0; vel.y = 0.0; vel.z = 0.0;
    for ( size_t i = 0; i < NL; ++i ) {
	    velocity[i].x = 0.0; velocity[i].y = 0.0; velocity[i].z = 0.0;
    }
    area = 0.0;
    volume = 0.0;
    dudx = 0.0; dudy = 0.0; dudz = 0.0;
    dvdx = 0.0; dvdy = 0.0; dvdz = 0.0;
    dwdx = 0.0; dwdy = 0.0; dwdz = 0.0;
    dTdx.resize(nmodes, 0.0); dTdy.resize(nmodes, 0.0); dTdz.resize(nmodes, 0.0);
    dtkedx = 0.0; dtkedy = 0.0; dtkedz = 0.0;
    domegadx = 0.0; domegady = 0.0; domegadz = 0.0;
    dfdx.resize(nsp, 0.0); dfdy.resize(nsp, 0.0); dfdz.resize(nsp, 0.0);
}

FV_Vertex::~FV_Vertex()
{
    dTdx.clear(); dTdy.clear(); dTdz.clear();
    dfdx.clear(); dfdy.clear(); dfdz.clear();
}

int FV_Vertex::copy_values_from(FV_Vertex &src)
{
    // don't copy id
    pos = src.pos; area = src.area; volume = src.volume;
    dudx = src.dudx; dudy = src.dudy; dudz = src.dudz;
    dvdx = src.dvdx; dvdy = src.dvdy; dvdz = src.dvdz;
    dwdx = src.dwdx; dwdy = src.dwdy; dwdz = src.dwdz;
    for ( size_t imode = 0; imode < dTdx.size(); ++imode ) {
	dTdx[imode] = src.dTdx[imode];
	dTdy[imode] = src.dTdy[imode];
	dTdz[imode] = src.dTdz[imode];
    }
    for ( size_t isp = 0; isp < dfdx.size(); ++isp ) {
	dfdx[isp] = src.dfdx[isp];
	dfdy[isp] = src.dfdy[isp];
	dfdz[isp] = src.dfdz[isp];
    }
    return SUCCESS;
}

//----------------------------------------------------------------------------

FV_Cell::FV_Cell(Gas_model *gm)
{
    id = 0;
    status = NORMAL_CELL;
    fr_reactions_allowed = 0;
    dt_chem = 0.0;
    dt_therm = 0.0;
    in_turbulent_zone = 0;
    base_qdot = 0.0;
    pos.x = 0.0; pos.y = 0.0; pos.z = 0.0;
    for ( size_t i = 0; i < NL; ++i ) {
	position[i].x = 0.0; position[i].y = 0.0; position[i].z = 0.0;
    }
    for ( size_t i = 0; i < NL; ++i ) {
	vol[i] = 0.0;
    }
    volume = 0.0;
    area = 0.0;
    for ( size_t i = 0; i < NL; ++i ) {
	ar[i] = 0.0;	
    }
    uf = 0.0;
    iLength = 0.0; jLength = 0.0; kLength = 0.0;
    L_min = 0.0;
    distance_to_nearest_wall = 0.0;
    half_cell_width_at_wall = 0.0;
    cell_at_nearest_wall = (FV_Cell *) NULL;
    for ( size_t i = 0; i < 6; ++i ) iface[i] = (FV_Interface *) NULL;
    for ( size_t i = 0; i < 8; ++i ) vtx[i] = (FV_Vertex *) NULL;
    fs = new FlowState(gm);
    U = new ConservedQuantities(gm);
    U_old = new ConservedQuantities(gm);
    for ( size_t i = 0; i < NL; ++i ) {
	dUdt[i] = new ConservedQuantities(gm);
    }
    Q = new ConservedQuantities(gm);
    Q_rad_org = 0.0;
    f_rad_org = 0.0;
    Q_rE_rad = 0.0;
    rho_at_start_of_step = 0.0;
    rE_at_start_of_step = 0.0;
#   if WITH_IMPLICIT == 1
    // Ojas, it may be good to initialize point-implicit variables.
    // double piM[6][6], pir[6][2], qL[6][2], g_Ll[6][2], gjvec[6][2], gjmtx[6][6];
#   endif
}

FV_Cell::~FV_Cell()
{
    delete fs;
    delete U;
    delete U_old;
    for ( size_t i = 0; i < NL; ++i ) {
	delete dUdt[i];
    }
    delete Q;
}

int FV_Cell::print()
{
    Gas_model *gmodel = get_gas_model_ptr();
    printf("----------- Begin data for cell -----------\n");
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();
    printf("nsp=%d, nmodes=%d\n", nsp, nmodes);
    printf("x=%e, y=%e, z=%e, base_qdot=%e\n", pos.x, pos.y, pos.z, base_qdot);
    fs->print();
    if ( get_radiation_flag() == 1 ) {
	printf("radiation source Q_rE_rad=%e\n", Q_rE_rad);
    }
    printf("Conserved quantities: \n");
    U->print();
    printf("status= %d \n", status);
    printf( "----------- End data for cell -----------\n");
    return SUCCESS;
}

int FV_Cell::point_is_inside(Vector3 &p, int dimensions)
/// \brief Returns 1 if the point p is inside or on the cell surface.
{
    if ( dimensions == 2 ) {
	// In 2 dimensions,
	// we split the x,y-plane into half-planes and check which side p is on.
	double xA = vtx[1]->pos.x; double yA = vtx[1]->pos.y;
	double xB = vtx[1]->pos.x; double yB = vtx[2]->pos.y;
	double xC = vtx[3]->pos.x; double yC = vtx[3]->pos.y;
	double xD = vtx[0]->pos.x; double yD = vtx[0]->pos.y;
	// Now, check to see if the specified point is on the
	// left of (or on) each bloundary line AB, BC, CD and DA.
	if ((p.x - xB) * (yA - yB) >= (p.y - yB) * (xA - xB) &&
	    (p.x - xC) * (yB - yC) >= (p.y - yC) * (xB - xC) &&
	    (p.x - xD) * (yC - yD) >= (p.y - yD) * (xC - xD) &&
	    (p.x - xA) * (yD - yA) >= (p.y - yA) * (xD - xA)) {
	    return 1;
	} else {
	    return 0;
	}
    } else {
	// In 3 dimensions,
	// the test consists of dividing the 6 cell faces into triangular facets
	// with outwardly-facing normals and then computing the volumes of the
	// tetrahedra formed by these facets and the sample point p.
	// If any of the tetrahedra volumes are positive
	// (i.e. p is on the positive side of a facet) and we assume a convex cell,
	// it means that the point is outside the cell and we may say so
	// without further testing.

	// North
	if ( tetrahedron_volume(vtx[2]->pos, vtx[3]->pos, vtx[7]->pos, p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[7]->pos, vtx[6]->pos, vtx[2]->pos, p) > 0.0 ) return 0;
	// East
	if ( tetrahedron_volume(vtx[1]->pos, vtx[2]->pos, vtx[6]->pos, p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[6]->pos, vtx[5]->pos, vtx[1]->pos, p) > 0.0 ) return 0;
	// South
	if ( tetrahedron_volume(vtx[0]->pos, vtx[1]->pos, vtx[5]->pos, p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[5]->pos, vtx[4]->pos, vtx[0]->pos, p) > 0.0 ) return 0;
	// West
	if ( tetrahedron_volume(vtx[3]->pos, vtx[0]->pos, vtx[4]->pos, p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[4]->pos, vtx[7]->pos, vtx[3]->pos, p) > 0.0 ) return 0;
	// Bottom
	if ( tetrahedron_volume(vtx[1]->pos, vtx[0]->pos, vtx[3]->pos, p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[3]->pos, vtx[2]->pos, vtx[1]->pos, p) > 0.0 ) return 0;
	// Top
	if ( tetrahedron_volume(vtx[4]->pos, vtx[5]->pos, vtx[6]->pos, p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[6]->pos, vtx[7]->pos, vtx[4]->pos, p) > 0.0 ) return 0;
	// If we arrive here, we haven't determined that the point is outside...
	return 1;
    } // end dimensions != 2
} // end point_is_inside()

int FV_Cell::copy_values_from(CFlowCondition &src)
{
    fs->gas->copy_values_from(*(src.gas));
    fs->vel.x = src.u; fs->vel.y = src.v; fs->vel.z = src.w;
    if ( get_mhd_flag() == 1 ) {
	fs->B.x = src.Bx; fs->B.y = src.By; fs->B.z = src.Bz;
    }
    fs->S = src.S;
    fs->mu_t = src.mu_t;
    fs->k_t = src.k_t;
    fs->tke = src.tke;
    fs->omega = src.omega;
    return SUCCESS;
}


/// \brief Copy the data from one FV_Cell structure to another.
///
/// \param src: pointer to the source cell structure
///\param type_of_copy indicates whether we copy the cell geometry along with the FlowState data
int FV_Cell::copy_values_from(FV_Cell &src, int type_of_copy)
{
    if ( type_of_copy == COPY_ALL_CELL_DATA ||
	 type_of_copy == COPY_FLOW_STATE ) {
	fs->copy_values_from(*(src.fs));
	status = src.status;
	Q_rE_rad = src.Q_rE_rad;
    }
    if ( type_of_copy == COPY_ALL_CELL_DATA ||
	 type_of_copy == COPY_CELL_LENGTHS ) {
       	iLength = src.iLength; jLength = src.jLength; kLength = src.kLength;
	pos.x = src.pos.x; pos.y = src.pos.y; pos.z = src.pos.z;
	if ( get_shock_fitting_flag() ) {
	    for ( size_t j = 0; j < 6; ++j ) {
		if ( src.iface[j] == 0 || iface[j] == 0 ) { // When copying from ghost cell which may
		    continue;                               // not have initialised interfaces
		}
		iface[j]->copy_values_from(*(src.iface[j]), COPY_ALL_CELL_DATA);
	    }
	}
    }
    return SUCCESS;
}


/// \brief Copy the cell data into a linear data buffer.
///
/// \param buf : pointer to the current element somewhere in buffer
/// \returns a pointer to the next available location in the data buffer.
double * FV_Cell::copy_values_to_buffer(double *buf, int type_of_copy)
{
    if (type_of_copy == COPY_ALL_CELL_DATA ||
	type_of_copy == COPY_FLOW_STATE) {
        buf = fs->copy_values_to_buffer(buf);
	*buf++ = (double)status;
	*buf++ = Q_rE_rad;
    }
    if (type_of_copy == COPY_ALL_CELL_DATA ||
        type_of_copy == COPY_CELL_LENGTHS) {
        *buf++ = iLength; *buf++ = jLength; *buf++ = kLength;
        *buf++ = pos.x; *buf++ = pos.y; *buf++ = pos.z;
	if ( get_shock_fitting_flag() ) {
	    for ( int j = 0; j < 4; ++j ) {
		if ( iface[j] == 0 ) { // When copying from ghost cell which may
		    continue;          // not have initialised interfaces
		}
		*buf++ = iface[j]->pos.x; *buf++ = iface[j]->pos.y; *buf++ = iface[j]->pos.z;
		*buf++ = iface[j]->vel.x; *buf++ = iface[j]->vel.y; *buf++ = iface[j]->vel.z;
		*buf++ = iface[j]->length;
		iface[j]->fs->copy_values_to_buffer(buf);
	    }
	}
    }
    return buf;
}

/// \brief Copy the data from a linear data buffer into the cell.
/// \param buf : pointer to the current element somewhere in buffer
/// \returns a pointer to the next available location in the data buffer.
double * FV_Cell::copy_values_from_buffer(double *buf, int type_of_copy)
{
    if (type_of_copy == COPY_ALL_CELL_DATA ||
	type_of_copy == COPY_FLOW_STATE) {
	buf = fs->copy_values_from_buffer(buf);
	status = (int)(*buf++);
	Q_rE_rad = *buf++;
    }
    if (type_of_copy == COPY_ALL_CELL_DATA ||
        type_of_copy == COPY_CELL_LENGTHS) {
        iLength = *buf++; jLength = *buf++; kLength = *buf++;
        pos.x = *buf++; pos.y = *buf++; pos.z = *buf++;
	if ( get_shock_fitting_flag() ) {
	    for ( int j = 0; j < 4; ++j ) {
		if ( iface[j] == 0 ) { // When copying from ghost cell which may
		    continue;          // not have initialised interfaces
		}
		iface[j]->pos.x = *buf++; iface[j]->pos.y = *buf++; iface[j]->pos.z = *buf++;
		iface[j]->vel.x = *buf++; iface[j]->vel.y = *buf++; iface[j]->vel.z = *buf++;
		iface[j]->length = *buf++;
		iface[j]->fs->copy_values_from_buffer(buf);
	    }
	}
    }
    return buf;
}


/// \brief Replace the flow data in a cell with the average from neighbour cells.
int FV_Cell::replace_flow_data_with_average(FV_Cell *src[], int ncell)
{
    int ii;
    if ( ncell < 1 ) return 0;  /* nothing to work with. */

    /* First, replace the crappy data with that from the first cell. */
    ii = 0;
    fs->gas->accumulate_values_from(*(src[ii]->fs->gas), 0.0);
    fs->vel.x = src[ii]->fs->vel.x;
    fs->vel.y = src[ii]->fs->vel.y;
    fs->vel.z = src[ii]->fs->vel.z;
    if ( get_mhd_flag() == 1 ) {
	fs->B.x = src[ii]->fs->B.x;
	fs->B.y = src[ii]->fs->B.y;
	fs->B.z = src[ii]->fs->B.z;
    }
    fs->mu_t = src[ii]->fs->mu_t;
    fs->k_t = src[ii]->fs->k_t;
    fs->tke = src[ii]->fs->tke;
    fs->omega = src[ii]->fs->omega;
    Q_rE_rad = src[ii]->Q_rE_rad;
    if ( ncell > 1 ) {
	/* Now, accumulate the remaining data with equal weight. */
	for ( ii = 1; ii < ncell; ++ii ) {
	    fs->gas->accumulate_values_from(*(src[ii]->fs->gas), 1.0 );
	    fs->vel.x += src[ii]->fs->vel.x;
	    fs->vel.y += src[ii]->fs->vel.y;
	    fs->vel.z += src[ii]->fs->vel.z;
	    if ( get_mhd_flag() == 1 ) {
		fs->B.x += src[ii]->fs->B.x;
		fs->B.y += src[ii]->fs->B.y;
		fs->B.z += src[ii]->fs->B.z;
	    }
	    fs->mu_t += src[ii]->fs->mu_t;
	    fs->k_t += src[ii]->fs->k_t;
	    fs->tke += src[ii]->fs->tke;
	    fs->omega += src[ii]->fs->omega;
	    Q_rE_rad += src[ii]->Q_rE_rad;
	}   /* end for */
	/* Effectively divides the result by ncell to get the average. */
	fs->gas->accumulate_values_from(*(fs->gas), (1.0-ncell)/((double)ncell) );
	fs->vel.x /= ((double)ncell);
	fs->vel.y /= ((double)ncell);
	fs->vel.z /= ((double)ncell);
	if ( get_mhd_flag() == 1 ) {
	    fs->B.x /= ((double)ncell);
	    fs->B.y /= ((double)ncell);
	    fs->B.z /= ((double)ncell);
	}
	fs->mu_t /= ((double)ncell);
	fs->k_t /= ((double)ncell);
	fs->tke /= ((double)ncell);
	fs->omega /= ((double)ncell);
	Q_rE_rad /= ((double)ncell);
    }
    // The following calls are expensive but getting to this point should be very rare.
    // If it is common, we have debugging to do...
    Gas_model *gmodel = get_gas_model_ptr();
    gmodel->eval_thermo_state_pT(*(fs->gas));
    if ( get_viscous_flag() ) gmodel->eval_transport_coefficients(*(fs->gas));
    if ( get_diffusion_flag() ) gmodel->eval_diffusion_coefficients(*(fs->gas));

    return SUCCESS;
} // end function replace_cell_data_with_average()


/// \brief Scan a string, extracting the flow solution (i.e. the primary variables).
int FV_Cell::scan_values_from_string(char *bufptr)
// There isn't any checking of the file content.
// If anything gets out of place, the result is wrong data.
{
    // Look for a new-line character and truncate the string there.
    char *cptr = strchr(bufptr, '\n');
    if ( cptr != NULL ) cptr = '\0';
    // Now, we should have a string with only numbers separated by spaces.
    pos.x = atof(strtok( bufptr, " " )); // tokenize on space characters
    pos.y = atof(strtok( NULL, " " ));
    pos.z = atof(strtok( NULL, " " ));
    volume = atof(strtok( NULL, " " ));
    fs->gas->rho = atof(strtok( NULL, " " ));
    fs->vel.x = atof(strtok( NULL, " " ));
    fs->vel.y = atof(strtok( NULL, " " ));
    fs->vel.z = atof(strtok( NULL, " " ));
    if ( get_mhd_flag() == 1 ) {
	fs->B.x = atof(strtok( NULL, " " ));
	fs->B.y = atof(strtok( NULL, " " ));
	fs->B.z = atof(strtok( NULL, " " ));
    }
    fs->gas->p = atof(strtok( NULL, " " ));
    fs->gas->a = atof(strtok( NULL, " " ));
    fs->gas->mu = atof(strtok( NULL, " " ));
    fs->gas->k[0] = atof(strtok( NULL, " " ));
    fs->mu_t = atof(strtok( NULL, " " ));
    fs->k_t = atof(strtok( NULL, " " ));
    fs->S = atoi(strtok( NULL, " " ));
    if ( get_radiation_flag() == 1 ) {
    	Q_rad_org = atof(strtok( NULL, " " ));
    	f_rad_org = atof(strtok( NULL, " " ));
	Q_rE_rad = atof(strtok( NULL, " " ));
    } else {
    	Q_rad_org = 0.0;
    	f_rad_org = 0.0;
	Q_rE_rad = 0.0;
    }
    fs->tke = atof(strtok( NULL, " " ));
    fs->omega = atof(strtok( NULL, " " ));
    size_t nsp = fs->gas->massf.size();
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	fs->gas->massf[isp] = atof(strtok( NULL, " " ));
    }
    if ( nsp > 1 ) dt_chem = atof(strtok( NULL, " " ));
    size_t nmodes = fs->gas->T.size();
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	fs->gas->e[imode] = atof(strtok( NULL, " " ));
	fs->gas->T[imode] = atof(strtok( NULL, " " ));
    }
    if ( nmodes > 1 ) dt_therm = atof(strtok( NULL, " " ));
    return SUCCESS;
} // end scan_values_from_string()

/// \brief Write the flow solution (i.e. the primary variables) to a string.
std::string FV_Cell::write_values_to_string()
{
    // The new format for Elmer3 puts everything onto one line.
    ostringstream ost;
    ost.setf(ios_base::scientific);
    ost.precision(12);
    ost << pos.x << " " << pos.y << " " << pos.z
	<< " " << volume << " " <<  fs->gas->rho
	<< " " << fs->vel.x << " " << fs->vel.y << " " << fs->vel.z;
    if ( get_mhd_flag() == 1 ) {
	ost << " " << fs->B.x << " " << fs->B.y << " " << fs->B.z;
    }
    ost << " " << fs->gas->p << " " << fs->gas->a << " " << fs->gas->mu
	<< " " << fs->gas->k[0] << " " << fs->mu_t << " " << fs->k_t
	<< " " << fs->S;
    if ( get_radiation_flag() == 1 ) {
	ost << " " << Q_rad_org << " " << f_rad_org << " " << Q_rE_rad;
    }
    ost << " " << fs->tke << " " << fs->omega;
    // Species mass fractions.
    size_t nsp = fs->gas->massf.size();
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	ost << " " << fs->gas->massf[isp];
    }
    if ( nsp > 1 ) ost << " " << dt_chem;
    // Individual energies (in e, T pairs)
    size_t nmodes = fs->gas->T.size();
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	ost << " " << fs->gas->e[imode] << " " << fs->gas->T[imode];
    }
    if ( nmodes > 1 ) ost << " " << dt_therm;
    // Don't put the newline char on the end.
    return ost.str();
} // end of write_values_to_string()

/// \brief Scan a string, extracting the discrete samples of the BGK velocity distribution function
int FV_Cell::scan_BGK_from_string(char *bufptr)
// There isn't any checking of the file content.
// If anything gets out of place, the result is wrong data.
{
    // Look for a new-line character and truncate the string there.
    char *cptr = strchr(bufptr, '\n');
    if ( cptr != NULL ) cptr = '\0';
    // Now, we should have a string with only numbers separated by spaces.
    // include the position data, duplicate of "flow", as insurance against 
    // finding this file in isolation
    pos.x = atof(strtok( bufptr, " " )); // tokenize on space characters
    pos.y = atof(strtok( NULL, " " ));
    pos.z = atof(strtok( NULL, " " ));
    volume = atof(strtok( NULL, " " ));
    
    
    // values of G and H are interleaved
    for (size_t iGH = 0; iGH < get_velocity_buckets(); ++iGH) {
	fs->G[iGH] = atof(strtok( NULL, " " ));
	fs->H[iGH] = atof(strtok( NULL, " " ));
    }

    return SUCCESS;
} // end scan_BGK_from_string()

/// \brief Write the discrete samples of the BGK velocity distribution function to a string.
std::string FV_Cell::write_BGK_to_string()
{
    // The new format for Elmer3 puts everything onto one line.
    ostringstream ost;
    ost.setf(ios_base::scientific);
    ost.precision(12);
    ost << pos.x << " " << pos.y << " " << pos.z << " " << volume;
    // BGK discrete samples of velocity distribution function
    // interleave G and H
    for ( size_t iGH = 0; iGH < get_velocity_buckets(); ++iGH ) {
	ost << " " << fs->G[iGH];
	ost << " " << fs->H[iGH];
    }
    // Don't put the newline char on the end.
    return ost.str();
} // end of write_BGK_to_string()

int FV_Cell::impose_chemistry_timestep(double dt)
{
    dt_chem = dt;
    return SUCCESS;
}


int FV_Cell::impose_thermal_timestep(double dt)
{
    dt_therm = dt;
    return SUCCESS;
}


int FV_Cell::set_fr_reactions_allowed(int flag)
{
    fr_reactions_allowed = flag;
    return SUCCESS;
}


int FV_Cell::record_conserved(void)
// Just in case they need to be reinstated later in the time step.
{
    U_old->copy_values_from(*U);
    return SUCCESS;
}


int FV_Cell::restore_conserved(void)
{
    U->copy_values_from(*U_old);
    fs->gas->rho = U->mass; // restore this copy of density also
    return SUCCESS;
}


int FV_Cell::encode_conserved(double omegaz)
{
    U->mass = fs->gas->rho;
    // X-, Y- and Z-momentum per unit volume.
    U->momentum.x = fs->gas->rho * fs->vel.x;
    U->momentum.y = fs->gas->rho * fs->vel.y;
    U->momentum.z = fs->gas->rho * fs->vel.z;
    // Magnetic field
    U->B.x = fs->B.x;
    U->B.y = fs->B.y;
    U->B.z = fs->B.z;
    // Total Energy / unit volume = density
    // (specific internal energy + kinetic energy/unit mass).
    double ke = 0.5 * (fs->vel.x * fs->vel.x
		       + fs->vel.y * fs->vel.y
		       + fs->vel.z * fs->vel.z);
    if ( get_k_omega_flag() ) {
	U->tke = fs->gas->rho * fs->tke;
	U->omega = fs->gas->rho * fs->omega;
	U->total_energy = fs->gas->rho * (fs->gas->e[0] + ke + fs->tke);
    } else {
	U->tke = 0.0;
	U->omega = fs->gas->rho * 1.0;
	U->total_energy = fs->gas->rho * (fs->gas->e[0] + ke);
    }
    if ( get_mhd_flag() == 1) {
	double me = 0.5 * (fs->B.x * fs->B.x
			   + fs->B.y * fs->B.y
			   + fs->B.z * fs->B.z);
	U->total_energy += me;
    }
    // Species densities: mass of species is per unit volume.
    for ( size_t isp = 0; isp < U->massf.size(); ++isp ) {
	U->massf[isp] = fs->gas->rho * fs->gas->massf[isp];
    }
    // Individual energies
    Gas_model *gmodel = get_gas_model_ptr();
    gmodel->encode_conserved_energy(*(fs->gas), U->energies);
    
    if ( omegaz != 0.0 ) {
	// Rotating frame.
	// Finally, we adjust the total energy to make rothalpy.
	// We do this last because the gas models don't know anything
	// about rotating frames and we don't want to mess their
	// energy calculations around.
	double rho = fs->gas->rho;
	double x = pos.x;
	double y = pos.y;
	double rsq = x*x + y*y;
	// The conserved quantity is rothalpy. I = E - (u**2)/2
	// where rotating frame velocity  u = omegaz * r.
	U->total_energy -= rho * 0.5 * omegaz * omegaz * rsq;
    }
    return SUCCESS;
} // end of encode_conserved()


int FV_Cell::decode_conserved(double omegaz)
{
    Gas_model *gmodel = get_gas_model_ptr();
    double ke, dinv, rE, me;

    // Mass / unit volume = Density
    double rho = U->mass;
    fs->gas->rho = rho;
    // This is limited to nonnegative and finite values.
    if ( get_bad_cell_complain_flag() && (rho <= 0.0) ) {
	printf("FV_Cell::decode_conserved(): Density is below minimum rho=%e\n", rho);
	printf("x=%g, y=%g, z=%g\n", pos.x, pos.y, pos.z);
	fs->gas->print_values();
    }
    dinv = 1.0 / rho;
    if ( omegaz != 0.0 ) {
	// Rotating frame.
	// The conserved quantity is rothalpy so we need to convert
	// back to enthalpy to do the rest of the decode.
	double x = pos.x;
	double y = pos.y;
	double rsq = x*x + y*y;
	rE = U->total_energy + rho * 0.5 * omegaz * omegaz * rsq;
    } else {
	// Non-rotating frame.
	rE = U->total_energy;
    }
    // Velocities from momenta.
    fs->vel.x = U->momentum.x * dinv;
    fs->vel.y = U->momentum.y * dinv;
    fs->vel.z = U->momentum.z * dinv;
    // Magnetic field
    fs->B.x = U->B.x;
    fs->B.y = U->B.y;
    fs->B.z = U->B.z;
    // Specific internal energy from total energy per unit volume.
    ke = 0.5 * (fs->vel.x * fs->vel.x + fs->vel.y * fs->vel.y + fs->vel.z * fs->vel.z);
    if (get_mhd_flag() == 1) {
        me = 0.5*(fs->B.x*fs->B.x + fs->B.y*fs->B.y + fs->B.z*fs->B.z);
    } else {
        me = 0.0;
    }
    if ( get_k_omega_flag() ) {
        fs->tke = U->tke * dinv;
        fs->omega = U->omega * dinv;
        fs->gas->e[0] = (rE - U->tke - me) * dinv - ke;
    } else {
        fs->tke = 0.0;
        fs->omega = 1.0;
        fs->gas->e[0] = (rE - me) * dinv - ke;
    }
    size_t nsp = U->massf.size();
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	fs->gas->massf[isp] = U->massf[isp] * dinv;
    }
    if ( nsp > 1 ) scale_mass_fractions( fs->gas->massf );
    
    // NOTE: - gas->e[0] is the total internal energy (sum of all modes)
    //         and has already been calculated above
    //       - renergies[0] is being calculated but never used.
    //         We've decided to leave it that way.
    double e0_save = fs->gas->e[0];
    gmodel->decode_conserved_energy(*(fs->gas), U->energies);
    fs->gas->e[0] = e0_save;
    // Fill out the other variables; P, T, a and
    // check the species mass fractions.
    // Update the viscous transport coefficients.
    gmodel->eval_thermo_state_rhoe(*(fs->gas));
    if ( get_viscous_flag() ) gmodel->eval_transport_coefficients(*(fs->gas));
    if ( get_diffusion_flag() ) gmodel->eval_diffusion_coefficients(*(fs->gas));

    return SUCCESS;
} // end of decode_conserved()

int FV_Cell::decode_conserved( int time_level, double omegaz)
{
    Gas_model *gmodel = get_gas_model_ptr();
    double ke, dinv, rE, me;
    
    // Mass / unit volume = Density
    double rho = U->mass;
    fs->gas->rho = rho;
    // This is limited to nonnegative and finite values.
    if ( get_bad_cell_complain_flag() && (rho <= 0.0) ) {
	printf("FV_Cell::decode_conserved(): Density is below minimum rho=%e\n", rho);
	printf("x=%g, y=%g, z=%g\n", pos.x, pos.y, pos.z);
	fs->gas->print_values();
    }
    dinv = 1.0 / rho;
    if ( omegaz != 0.0 ) {
	// Rotating frame.
	// The conserved quantity is rothalpy so we need to convert
	// back to enthalpy to do the rest of the decode.
	double x = position[time_level].x;
	double y = position[time_level].y;
	double rsq = x*x + y*y;
	rE = U->total_energy + rho * 0.5 * omegaz * omegaz * rsq;
    } else {
	// Non-rotating frame.
	rE = U->total_energy;
    }
    // Velocities from momenta.
    fs->vel.x = U->momentum.x * dinv;
    fs->vel.y = U->momentum.y * dinv;
    fs->vel.z = U->momentum.z * dinv;
    // Magnetic field
    fs->B.x = U->B.x;
    fs->B.y = U->B.y;
    fs->B.z = U->B.z;
    // Specific internal energy from total energy per unit volume.
    ke = 0.5 * (fs->vel.x * fs->vel.x + fs->vel.y * fs->vel.y + fs->vel.z * fs->vel.z);
    if (get_mhd_flag() == 1) {
        me = 0.5*(fs->B.x*fs->B.x + fs->B.y*fs->B.y + fs->B.z*fs->B.z);
    } else {
        me = 0.0;
    }
    if ( get_k_omega_flag() ) {
        fs->tke = U->tke * dinv;
        fs->omega = U->omega * dinv;
        fs->gas->e[0] = (rE - U->tke - me) * dinv - ke;
    } else {
        fs->tke = 0.0;
        fs->omega = 1.0;
        fs->gas->e[0] = (rE - me) * dinv - ke;
    }
    size_t nsp = U->massf.size();
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	fs->gas->massf[isp] = U->massf[isp] * dinv;
    }
    if ( nsp > 1 ) scale_mass_fractions( fs->gas->massf );

    // NOTE: - gas->e[0] is the total internal energy (sum of all modes)
    //         and has already been calculated above
    //       - renergies[0] is being calculated but never used.
    //         We've decided to leave it that way.
    double e0_save = fs->gas->e[0];
    gmodel->decode_conserved_energy(*(fs->gas), U->energies);
    fs->gas->e[0] = e0_save;
    // Fill out the other variables; P, T, a and
    // check the species mass fractions.
    // Update the viscous transport coefficients.
    gmodel->eval_thermo_state_rhoe(*(fs->gas));
    if ( get_viscous_flag() ) gmodel->eval_transport_coefficients(*(fs->gas));
    if ( get_diffusion_flag() ) gmodel->eval_diffusion_coefficients(*(fs->gas));
    
    return SUCCESS;
} // end of decode_conserved()


/// \brief Check the primary flow data for a specified cell.
///  \returns 1 for valid data, 0 for bad data.
int FV_Cell::check_flow_data(void)
{
    int data_valid = 1;
    bool print_message = get_bad_cell_complain_flag() == 1;
    
    data_valid = fs->gas->check_values(print_message);
#   define MAXVEL 30000.0
    if (fabs(fs->vel.x) > MAXVEL || fabs(fs->vel.y) > MAXVEL || fabs(fs->vel.z) > MAXVEL) {
	if ( print_message )
	    cout << "Velocity bad " << fs->vel.x << " " << fs->vel.y << " " << fs->vel.z << endl;
	data_valid = 0;
    }
    if ( !(finite(fs->tke)) ) {
	if ( print_message ) cout << "Turbulence KE invalid number " << fs->tke << endl;
	data_valid = 0;
    }
    if ( fs->tke < 0.0 ) {
	if ( print_message ) cout << "Turbulence KE negative " << fs->tke << endl;
	data_valid = 0;
    }
    if ( !(finite(fs->omega)) ) {
	if ( print_message ) cout << "Turbulence frequency invalid number " << fs->omega << endl;
	data_valid = 0;
    }
    if ( fs->omega <= 0.0 ) {
	if ( print_message ) cout << "Turbulence frequency nonpositive " << fs->omega << endl;
	data_valid = 0;
    }
    if ( !data_valid && print_message ) {
	cout << "cell pos=(" << pos.x << "," << pos.y << "," << pos.z << ")" << endl;
	fs->print();
	cout << "----------------------------------------------------------" << endl;
    }
    return data_valid;
} // end of check_flow_data()

int FV_Cell::set_geometry_to_time_level( void )
{
    for ( size_t j = 0; j < NL; ++j ) {
	ar[j] = area;
	vol[j] = volume;
	position[j] = pos;
	for ( int i = 0; i <= 3; ++i ) {
	    vtx[i]->position[j] = vtx[i]->pos;
	    iface[i]->ar[j] = iface[i]->area;
	}
    }
    return SUCCESS;
}

/// \brief Set the cell geometry to the values calculated at the specified time level.
///
int FV_Cell::set_geometry_from_time_level(int time_level)
{
    area = ar[time_level];
    volume = vol[time_level];
    pos = position[time_level];
    for ( int i = 0; i <= 3; ++i ) {
	vtx[i]->pos = vtx[i]->position[time_level];
	iface[i]->area = iface[i]->ar[time_level];
    }
    return SUCCESS;
}

/// \brief Compute the time derivatives for the conserved quantities.
///
/// These are the spatial (RHS) terms in the semi-discrete governing equations.
/// \param time_level : specifies where the computed derivatives are to be stored.
/// \param dimensions : number of space dimensions: 2 for mbcns, 3 for eilmer
int FV_Cell::time_derivatives( int time_level, int dimensions )
{
    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t nmodes = gmodel->get_number_of_modes();
    FV_Interface *IFn = iface[NORTH];
    FV_Interface *IFe = iface[EAST];
    FV_Interface *IFs = iface[SOUTH];
    FV_Interface *IFw = iface[WEST];
    FV_Interface *IFt = iface[TOP];
    FV_Interface *IFb = iface[BOTTOM];
    // Cell volume (inverted).
    double vol_inv = 1.0 / vol[time_level];
    double integral;

    if (status == MASKED_CELL || status == SHADOWED_CELL) {
	dUdt[time_level]->clear_values();
    	return SUCCESS; // do not update cell if it is covered with a piston
    }
    
    // Time-derivative for Mass/unit volume.
    // Note that the unit normals for the interfaces are oriented
    // such that the unit normals for the east, north and top faces
    // are outward and the unit normals for the south, west and
    // bottom faces are inward.
    integral = -IFe->F->mass * IFe->ar[time_level] - IFn->F->mass * IFn->ar[time_level]
	+ IFw->F->mass * IFw->ar[time_level] + IFs->F->mass * IFs->ar[time_level];
    if ( dimensions == 3 )
	integral += IFb->F->mass * IFb->ar[time_level] - IFt->F->mass * IFt->ar[time_level];
    dUdt[time_level]->mass = vol_inv * integral + Q->mass;

    // Time-derivative for X-Momentum/unit volume.
    integral = -IFe->F->momentum.x * IFe->ar[time_level] - IFn->F->momentum.x * IFn->ar[time_level]
	+ IFw->F->momentum.x * IFw->ar[time_level] + IFs->F->momentum.x * IFs->ar[time_level];
    if ( dimensions == 3 )
	integral += IFb->F->momentum.x * IFb->ar[time_level] - IFt->F->momentum.x * IFt->ar[time_level];
    dUdt[time_level]->momentum.x = vol_inv * integral + Q->momentum.x;
    // Time-derivative for Y-Momentum/unit volume.
    integral = -IFe->F->momentum.y * IFe->ar[time_level] - IFn->F->momentum.y * IFn->ar[time_level]
	+ IFw->F->momentum.y * IFw->ar[time_level] + IFs->F->momentum.y * IFs->ar[time_level];
    if ( dimensions == 3 )
	integral += IFb->F->momentum.y * IFb->ar[time_level] - IFt->F->momentum.y * IFt->ar[time_level];
    dUdt[time_level]->momentum.y = vol_inv * integral + Q->momentum.y;
    
    // we require the z-momentum for MHD even in 2D
    if ((dimensions == 3) || (get_mhd_flag() == 1)) {
	// Time-derivative for Z-Momentum/unit volume.
	integral = -IFe->F->momentum.z * IFe->ar[time_level] - IFn->F->momentum.z * IFn->ar[time_level]
	    + IFw->F->momentum.z * IFw->ar[time_level] + IFs->F->momentum.z * IFs->ar[time_level];
    }
    if ( dimensions == 3) {
	integral += IFb->F->momentum.z * IFb->ar[time_level] - IFt->F->momentum.z * IFt->ar[time_level];
    }
    if ((dimensions == 3) || (get_mhd_flag() == 1)) {
	dUdt[time_level]->momentum.z = vol_inv * integral + Q->momentum.z;
    } else {
	dUdt[time_level]->momentum.z = 0.0;
    }
    
    if (get_mhd_flag() == 1) {
	// Time-derivative for X-Magnetic Field/unit volume.
	integral = -IFe->F->B.x * IFe->ar[time_level] - IFn->F->B.x * IFn->ar[time_level]
	    + IFw->F->B.x * IFw->ar[time_level] + IFs->F->B.x * IFs->ar[time_level];
	if ( dimensions == 3 )
	    integral += IFb->F->B.x * IFb->ar[time_level] - IFt->F->B.x * IFt->ar[time_level];
	dUdt[time_level]->B.x = vol_inv * integral + Q->B.x;
	// Time-derivative for Y-Magnetic Field/unit volume.
	integral = -IFe->F->B.y * IFe->ar[time_level] - IFn->F->B.y * IFn->ar[time_level]
	    + IFw->F->B.y * IFw->ar[time_level] + IFs->F->B.y * IFs->ar[time_level];
	if ( dimensions == 3 )
	    integral += IFb->F->B.y * IFb->ar[time_level] - IFt->F->B.y * IFt->ar[time_level];
	dUdt[time_level]->B.y = vol_inv * integral + Q->B.y;
	// Time-derivative for Z-Magnetic Field/unit volume.
	integral = -IFe->F->B.z * IFe->ar[time_level] - IFn->F->B.z * IFn->ar[time_level]
	    + IFw->F->B.z * IFw->ar[time_level] + IFs->F->B.z * IFs->ar[time_level];
	if ( dimensions == 3 ) {
	    integral += IFb->F->B.z * IFb->ar[time_level] - IFt->F->B.z * IFt->ar[time_level];
	}
	dUdt[time_level]->B.z = vol_inv * integral + Q->B.z;
    }
    else {
	dUdt[time_level]->B.x = 0.0;
	dUdt[time_level]->B.y = 0.0;
	dUdt[time_level]->B.z = 0.0;
    }

    // Time-derivative for Total Energy/unit volume.
    integral = -IFe->F->total_energy * IFe->ar[time_level] - IFn->F->total_energy * IFn->ar[time_level]
	+ IFw->F->total_energy * IFw->ar[time_level] + IFs->F->total_energy * IFs->ar[time_level];
    if ( dimensions == 3 )
	integral += IFb->F->total_energy * IFb->ar[time_level] - IFt->F->total_energy * IFt->ar[time_level];
    dUdt[time_level]->total_energy = vol_inv * integral + Q->total_energy;
    
    if ( get_k_omega_flag() ) {
	integral = -IFe->F->tke * IFe->ar[time_level] - IFn->F->tke * IFn->ar[time_level]
	    + IFw->F->tke * IFw->ar[time_level] + IFs->F->tke * IFs->ar[time_level];
	if ( dimensions == 3 )
	    integral += IFb->F->tke * IFb->ar[time_level] - IFt->F->tke * IFt->ar[time_level];
	dUdt[time_level]->tke = vol_inv * integral + Q->tke;
	
	integral = -IFe->F->omega * IFe->ar[time_level] - IFn->F->omega * IFn->ar[time_level]
	    + IFw->F->omega * IFw->ar[time_level] + IFs->F->omega * IFs->ar[time_level];
	if ( dimensions == 3 )
	    integral += IFb->F->omega * IFb->ar[time_level] - IFt->F->omega * IFt->ar[time_level];
	dUdt[time_level]->omega = vol_inv * integral + Q->omega;
    } else {
	dUdt[time_level]->tke = 0.0;
	dUdt[time_level]->omega = 0.0;
    }
    // Time-derivative for individual species.
    // The conserved quantity is the mass per unit
    // volume of species isp and
    // the fluxes are mass/unit-time/unit-area.
    // Units of DmassfDt are 1/sec.
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	integral =
	    -IFe->F->massf[isp] * IFe->ar[time_level]
	    - IFn->F->massf[isp] * IFn->ar[time_level]
	    + IFw->F->massf[isp] * IFw->ar[time_level]
	    + IFs->F->massf[isp] * IFs->ar[time_level];
	if ( dimensions == 3 )
	    integral += IFb->F->massf[isp] * IFb->ar[time_level] - IFt->F->massf[isp] * IFt->ar[time_level];
	dUdt[time_level]->massf[isp] = vol_inv * integral + Q->massf[isp];
    }
    // Individual energies.
    // NOTE: energies[0] is never used so skipping (DFP 10/12/09)
    for ( size_t imode = 1; imode < nmodes; ++imode ) {
	integral =
	    -IFe->F->energies[imode] * IFe->ar[time_level]
	    - IFn->F->energies[imode] * IFn->ar[time_level]
	    + IFw->F->energies[imode] * IFw->ar[time_level]
	    + IFs->F->energies[imode] * IFs->ar[time_level];
	if ( dimensions == 3 )
	    integral += IFb->F->energies[imode] * IFb->ar[time_level] - IFt->F->energies[imode] * IFt->ar[time_level];
	dUdt[time_level]->energies[imode] = vol_inv * integral + Q->energies[imode];
    }
    return SUCCESS;
} // end of time_derivatives()


/// \brief Apply the predictor-stage of the update for a specified cell.
/// \param dt   : size of the time-step
int FV_Cell::predictor_update(double dt)
{
    ConservedQuantities &dUdt0 = *(dUdt[0]);
    double gamma_1;
    double vr = 1.0;
    if (get_Torder_flag() == 3) {
	/* 3rd order Runge-Kutta */
	gamma_1 = 8.0 / 15.0;
    } else {
	/* Normal Predictor-Corrector or Euler */
	gamma_1 = 1.0;
    }
    if ( get_shock_fitting_flag() == 1 ) {
        vr = vol[0]/vol[1];
    }
    U->mass = vr * (U_old->mass + dt * gamma_1 * dUdt0.mass);
    U->momentum.x = vr * (U_old->momentum.x + dt * gamma_1 * dUdt0.momentum.x);
    U->momentum.y = vr * (U_old->momentum.y + dt * gamma_1 * dUdt0.momentum.y);
    U->momentum.z = vr * (U_old->momentum.z + dt * gamma_1 * dUdt0.momentum.z);

    // Magnetic field
    if (get_mhd_flag() == 1) {
	U->B.x = vr * (U_old->B.x + dt * gamma_1 * dUdt0.B.x);
	U->B.y = vr * (U_old->B.y + dt * gamma_1 * dUdt0.B.y);
	U->B.z = vr * (U_old->B.z + dt * gamma_1 * dUdt0.B.z);
    }
    
    U->total_energy = vr * (U_old->total_energy + dt * gamma_1 * dUdt0.total_energy);
    if ( get_k_omega_flag() ) {
	U->tke = vr * (U_old->tke + dt * gamma_1 * dUdt0.tke);
	U->tke = MAXIMUM(U->tke, 0.0);
	U->omega = vr * (U_old->omega + dt * gamma_1 * dUdt0.omega);
	U->omega = MAXIMUM(U->omega, U_old->mass);
	// ...assuming a minimum value of 1.0 for omega
	// It may occur (near steps in the wall) that a large flux of romega
	// through one of the cell interfaces causes romega within the cell
	// to drop rapidly.
	// The large values of omega come from Menter's near-wall correction that may be
	// applied outside the control of this finite-volume core code.
	// These large values of omega will be convected along the wall and,
	// if they are convected past a corner with a strong expansion,
	// there will be an unreasonably-large flux out of the cell.
    } else {
	U->tke = U_old->tke;
	U->omega = U_old->omega;
    }
    for ( size_t isp = 0; isp < U->massf.size(); ++isp ) {
	U->massf[isp] = vr * (U_old->massf[isp] + dt * gamma_1 * dUdt0.massf[isp]);
    }
    // NOTE: energies[0] is never used so skipping (DFP 10/12/09)
    for ( size_t imode = 1; imode < U->energies.size(); ++imode ) {
	U->energies[imode] = vr * (U_old->energies[imode] + dt * gamma_1 * dUdt0.energies[imode]);
    }
    return SUCCESS;
} // end of predictor_update()


/// \brief Apply the corrector-stage of the update for a specified cell.
/// \param dt   : size of the time-step
int FV_Cell::corrector_update(double dt)
{
    ConservedQuantities &dUdt0 = *(dUdt[0]);
    ConservedQuantities &dUdt1 = *(dUdt[1]);
    double th, th_inv;
    double v0 = 1.0;
    double v1 = 1.0;
    double vol_inv = 1.0;
    /*
     * Set the type of time-stepping
     * th = 0.0 : Euler
     * th = 0.5 : 2nd order
     * th = 1.0 : sort-of-implicit?
     */
    if ( get_shock_fitting_flag() == 1 ) {
        v0 = vol[0]; // The volume at different 
	v1 = vol[1];
	vol_inv = 1.0 / vol[2];
    }
    if (get_Torder_flag() == 3) {
	/* 3rd order Runge-Kutta */
	th = v1 * 5.0 / 12.0;
	th_inv = v0 * -17.0 / 60.0;
    } else {
	/* Normal Predictor-Corrector or Euler */
	th = v1 * 0.5;
	th_inv = v0 * 0.5;
    }
    
    U->mass = vol_inv * (v0 * U_old->mass + dt * (th_inv * dUdt0.mass + th * dUdt1.mass));
    U->momentum.x = vol_inv * (v0 * U_old->momentum.x + 
			       dt * (th_inv * dUdt0.momentum.x + th * dUdt1.momentum.x));
    U->momentum.y = vol_inv * (v0 * U_old->momentum.y + 
			       dt * (th_inv * dUdt0.momentum.y + th * dUdt1.momentum.y));
    U->momentum.z = vol_inv * (v0 * U_old->momentum.z + 
			       dt * (th_inv * dUdt0.momentum.z + th * dUdt1.momentum.z));

    // Magnetic field
    if (get_mhd_flag() == 1) {
	U->B.x = vol_inv * (v0 * U_old->B.x + dt * (th_inv * dUdt0.B.x + th * dUdt1.B.x));
	U->B.y = vol_inv * (v0 * U_old->B.y + dt * (th_inv * dUdt0.B.y + th * dUdt1.B.y));
	U->B.z = vol_inv * (v0 * U_old->B.z + dt * (th_inv * dUdt0.B.z + th * dUdt1.B.z));
    }
    
    U->total_energy = vol_inv * (v0 * U_old->total_energy + 
				 dt * (th_inv * dUdt0.total_energy + th * dUdt1.total_energy));
    if ( get_k_omega_flag() ) {
	U->tke = vol_inv * (v0 * U_old->tke + dt * (th_inv * dUdt0.tke + th * dUdt1.tke));
	U->tke = MAXIMUM(U->tke, 0.0);
	U->omega = vol_inv * (v0 * U_old->omega + dt * (th_inv * dUdt0.omega + th * dUdt1.omega));
	U->omega = MAXIMUM(U->omega, U_old->mass);
    } else {
	U->tke = vol_inv * (v0 * U_old->tke);
	U->omega = vol_inv * (v0 * U_old->omega);
    }
    for ( size_t isp = 0; isp < U->massf.size(); ++isp ) {
	U->massf[isp] = vol_inv * (v0 * U_old->massf[isp] +
				   dt * (th_inv * dUdt0.massf[isp] + th * dUdt1.massf[isp]));
    }
    for ( size_t imode = 0; imode < U->energies.size(); ++imode ) {
	U->energies[imode] = vol_inv * (v0 * U_old->energies[imode] +
					dt * (th_inv * dUdt0.energies[imode] + th * dUdt1.energies[imode]));
    }
    return SUCCESS;
} // end of corrector_update()


/// \brief Apply the final Runge-Kutta (3rd order) stage of the update for a specified cell.
/// \param dt : size of the time-step
int FV_Cell::rk3_update(double dt)
{
    ConservedQuantities &dUdt1 = *(dUdt[1]);
    ConservedQuantities &dUdt2 = *(dUdt[2]);
    double gamma_3 = 3.0 / 4.0;
    double psi_2 = -5.0 / 12.0;
    double vr = 1.0;
    if ( get_shock_fitting_flag() == 1 ) {
        vr = volume/vol[2];
    }
    U->mass = vr * (U_old->mass + dt * (psi_2 * dUdt1.mass + gamma_3 * dUdt2.mass));
    U->momentum.x = vr * (U_old->momentum.x + dt * (psi_2 * dUdt1.momentum.x + gamma_3 * dUdt2.momentum.x));
    U->momentum.y = vr * (U_old->momentum.y + dt * (psi_2 * dUdt1.momentum.y + gamma_3 * dUdt2.momentum.y));
    U->momentum.z = vr * (U_old->momentum.z + dt * (psi_2 * dUdt1.momentum.z + gamma_3 * dUdt2.momentum.z));
    
    // Magnetic field
    if (get_mhd_flag() == 1) {
	U->B.x = vr * (U_old->B.x + dt * (psi_2 * dUdt1.B.x + gamma_3 * dUdt2.B.x));
	U->B.y = vr * (U_old->B.y + dt * (psi_2 * dUdt1.B.y + gamma_3 * dUdt2.B.y));
	U->B.z = vr * (U_old->B.z + dt * (psi_2 * dUdt1.B.z + gamma_3 * dUdt2.B.z));
    }

    U->total_energy = vr * (U_old->total_energy + 
			    dt * (psi_2 * dUdt1.total_energy + gamma_3 * dUdt2.total_energy));
    if ( get_k_omega_flag() ) {
	U->tke = vr * (U_old->tke + dt * (psi_2 * dUdt1.tke + gamma_3 * dUdt2.tke));
	U->tke = MAXIMUM(U->tke, 0.0);
	U->omega = vr * (U_old->omega + dt * (psi_2 * dUdt1.omega + gamma_3 * dUdt2.omega));
	U->omega = MAXIMUM(U->omega, U_old->mass);
    } else {
	U->tke = vr * (U_old->tke);
	U->omega = vr * (U_old->omega);
    }
    for ( size_t isp = 0; isp < U->massf.size(); ++isp ) {
	U->massf[isp] = vr * (U_old->massf[isp] +
	    dt * (psi_2 * dUdt1.massf[isp] + gamma_3 * dUdt2.massf[isp]));
    }
    for ( size_t imode = 0; imode < U->energies.size(); ++imode ) {
	U->energies[imode] = vr * (U_old->energies[imode] +
	    dt * (psi_2 * dUdt1.energies[imode] + gamma_3 * dUdt2.energies[imode]));
    }
    return SUCCESS;
} // end of rk3_update()


/// \brief Apply the chemistry update for a specified cell.
///
/// Use the finite-rate chemistry module to update the
/// species fractions and the other thermochemical properties.
/// \param dt   : size of the time-step
int FV_Cell::chemical_increment(double dt)
{
    if ( !fr_reactions_allowed ) return SUCCESS;
    Gas_model *gmodel = get_gas_model_ptr();
    Reaction_update *rupdate = get_reaction_update_ptr();

    int flag = rupdate->update_state(*(fs->gas), dt, dt_chem, gmodel);

    // The update only changes mass fractions, we need to impose
    // a thermodynamic constraint based on a call to the equation
    // of state.
    gmodel->eval_thermo_state_rhoe(*(fs->gas));

    // If we are doing a viscous sim, we'll need to ensure
    // viscous properties are up-to-date
    if ( get_viscous_flag() ) gmodel->eval_transport_coefficients(*(fs->gas));
    if ( get_diffusion_flag() ) gmodel->eval_diffusion_coefficients(*(fs->gas));
    // ...but we have to manually update the conservation quantities
    // for the gas-dynamics time integration.
    // Species densities: mass of species isp per unit volume.
    for ( size_t isp = 0; isp < fs->gas->massf.size(); ++isp )
	U->massf[isp] = fs->gas->rho * fs->gas->massf[isp];
    return flag;
} // end of chemical_increment()


/// \brief Apply the thermal update for a specified cell.
///
/// Use the nonequilibrium multi-Temperature module to update the
/// energy values  and the other thermochemical properties.
/// \param dt   : size of the time-step
int FV_Cell::thermal_increment(double dt)
{
    if ( !fr_reactions_allowed ) return SUCCESS;
    Gas_model *gmodel = get_gas_model_ptr();
    Energy_exchange_update *eeupdate = get_energy_exchange_update_ptr();

    int flag = eeupdate->update_state(*(fs->gas), dt, dt_therm, gmodel);

    // The update only changes modal energies, we need to impose
    // a thermodynamic constraint based on a call to the equation
    // of state.
    gmodel->eval_thermo_state_rhoe(*(fs->gas));

    // If we are doing a viscous sim, we'll need to ensure
    // viscous properties are up-to-date
    if ( get_viscous_flag() ) gmodel->eval_transport_coefficients(*(fs->gas));
    if ( get_diffusion_flag() ) gmodel->eval_diffusion_coefficients(*(fs->gas));
    // ...but we have to manually update the conservation quantities
    // for the gas-dynamics time integration.
    // Independent energies energy: Joules per unit volume.
    gmodel->encode_conserved_energy(*(fs->gas), U->energies);
    return flag;
} // end of thermal_increment()


/// \brief Compute the signal frequency (1/seconds) for a cell.
///
/// \param dimensions: number of spatial dimensions: 2 for mbcns2, 3 for eilmer
/// The North and East faces are taken as the representative lengths the cells.
double FV_Cell::signal_frequency(int dimensions)
{
    double signal;
    double un_N, un_E, un_T, u_mag;
    double Bn_N, Bn_E, Bn_T, B_mag;
    double ca2, cfast;
    double gam_eff, viscous_factor;
    int statusf;
    Gas_model *gmodel = get_gas_model_ptr();
    FV_Interface *north = iface[NORTH];
    FV_Interface *east = iface[EAST];
    FV_Interface *top = iface[TOP];
    // Get the local normal velocities by rotating the
    // local frame of reference.
    // Also, compute the velocity magnitude and
    // recall the minimum length.
    un_N = fabs(dot(fs->vel, north->n));
    un_E = fabs(dot(fs->vel, east->n));
    if ( dimensions == 3 ) {
	un_T = fabs(dot(fs->vel, top->n));
	u_mag = sqrt(fs->vel.x * fs->vel.x + fs->vel.y * fs->vel.y + fs->vel.z * fs->vel.z);
    }  else {
	un_T = 0.0;
	u_mag = sqrt(fs->vel.x * fs->vel.x + fs->vel.y * fs->vel.y);
    }
    if (get_mhd_flag() == 1) {
	Bn_N = fabs(dot(fs->B, north->n));
	Bn_E = fabs(dot(fs->B, east->n));
	if ( dimensions == 3 ) {
	    Bn_T = fabs(dot(fs->B, top->n));
	}
	u_mag = sqrt(fs->vel.x * fs->vel.x + fs->vel.y * fs->vel.y + fs->vel.z * fs->vel.z);
	B_mag = sqrt(fs->B.x * fs->B.x + fs->B.y * fs->B.y + fs->B.z * fs->B.z);
    }
    // Check the INVISCID time step limit first,
    // then add a component to ensure viscous stability.
    if ( get_stringent_cfl_flag() ) {
	// Make the worst case.
	if (get_mhd_flag() == 1) {
	    // MHD
	    ca2 = B_mag*B_mag / fs->gas->rho;
	    cfast = sqrt( ca2 + fs->gas->a * fs->gas->a );
            signal = (u_mag + cfast) / L_min;
	}
	else {
	    // Hydrodynamics
	    signal = (u_mag + fs->gas->a) / L_min;
	}
    } else {
	// Standard signal speeds along each face.
	double signalN, signalE, signalT;
	if (get_mhd_flag() == 1) {
	    double ca2, catang2_N, catang2_E, cfast_N, cfast_E;
	    ca2 = B_mag * B_mag / fs->gas->rho;
	    ca2 = ca2 + fs->gas->a * fs->gas->a;
	    catang2_N = Bn_N * Bn_N / fs->gas->rho;
	    cfast_N = 0.5 * ( ca2 + sqrt( ca2*ca2 - 4.0 * (fs->gas->a * fs->gas->a * catang2_N) ) );
	    cfast_N = sqrt(cfast_N);
	    catang2_E = Bn_E * Bn_E / fs->gas->rho;
	    cfast_E = 0.5 * ( ca2 + sqrt( ca2*ca2 - 4.0 * (fs->gas->a * fs->gas->a * catang2_E) ) );
	    cfast_E = sqrt(cfast_E);
	    if ( dimensions == 3 ) {
		double catang2_T, cfast_T, signalT;
		catang2_T = Bn_T * Bn_T / fs->gas->rho;
		cfast_T = 0.5 * ( ca2 + sqrt( ca2*ca2 - 4.0 * (fs->gas->a * fs->gas->a * catang2_T) ) );
		cfast_T = sqrt(cfast_T);
		signalN = (un_N + cfast_N) / jLength;
		signal = signalN;
		signalE = (un_E + cfast_E) / iLength;
		if ( signalE > signal ) signal = signalE;
		signalT = (un_T + cfast_T) / kLength;
		if ( signalT > signal ) signal = signalT;
	    }
	    else {
		signalN = (un_N + cfast) / east->length;
		signalE = (un_E + cfast) / north->length;
		signal = MAXIMUM(signalN, signalE);
	    }
	}
	else if ( dimensions == 3 ) {
	    // eilmer -- 3D cells
	    signalN = (un_N + fs->gas->a) / jLength;
	    signal = signalN;
	    signalE = (un_E + fs->gas->a) / iLength;
	    if ( signalE > signal ) signal = signalE;
	    signalT = (un_T + fs->gas->a) / kLength;
	    if ( signalT > signal ) signal = signalT;
	} else {
	    // mbcns2 -- 2D cells
	    // The velocity normal to the north face is assumed to run
	    // along the length of the east face.
	    signalN = (un_N + fs->gas->a) / east->length;
	    signalE = (un_E + fs->gas->a) / north->length;
	    signal = MAXIMUM(signalN, signalE);
	}
    }
    if ( get_viscous_flag() == 1 && get_implicit_flag() == 0 && fs->gas->mu > 10e-23) {
	// Factor for the viscous time limit.
	// This factor is not included if viscosity is zero.
#	if VISCOUS_TIME_LIMIT_MODEL==0
	// See Swanson, Turkel and White (1991)
	gam_eff = gmodel->gamma(*(fs->gas), statusf);
	// Need to sum conductivities for TNE
	double k_total = 0.0;
	for ( size_t i=0; i<fs->gas->k.size(); ++i ) k_total += fs->gas->k[i];
	double Prandtl = fs->gas->mu * gmodel->Cp(*(fs->gas), statusf) / k_total;
	viscous_factor = get_viscous_factor();
	if ( dimensions == 3 ) {
	    signal += 4.0 * viscous_factor * (fs->gas->mu + fs->mu_t)
		* gam_eff / (Prandtl * fs->gas->rho)
		* (1.0 / (iLength * iLength) + 1.0 / (jLength * jLength) + 1.0 / (kLength * kLength));
	} else {
	    signal += 4.0 * viscous_factor * (fs->gas->mu + fs->mu_t) * gam_eff / (Prandtl * fs->gas->rho)
		* (1.0 / (east->length * east->length) + 1.0 / (north->length * north->length));
	}
#	elif VISCOUS_TIME_LIMIT_MODEL==1
	// A viscous time limit model incorporating diffusion effects
	// See Ramshaw and Chang PCPP V.12 n.3 1992 p314
	// 1. Viscosity signal frequency
	double D_a = 4.0 * viscous_factor * (gas->mu + -0.66667 * gas->mu +
					     mu_t - 0.66667 *mu_t ) / gas->rho;
	// 1. Conductivity signal frequency
	double D_b = gas->k[0];
	for ( int imode=1; imode<gmodel->get_number_of_modes(); ++imode )
	    D_b += gas->k[imode];
	int status;
	double c_v = gmodel->Cv( gas, status );
	D_b *= viscous_factor / ( gas->rho * c_v );
	// 3. calculate the largest mixture diffusivity
	double D_c = 0.0;
	if ( get_diffusion_flag() == 1 ) {
	    int nsp = gmodel->get_number_of_species();
	    vector<double> x(nsp);
	    vector<double> M(nsp);
	    vector<double> DAV_im(nsp);
	    for ( int isp=0; isp<nsp; ++isp )
	    	M[isp] = gmodel->molecular_weight(isp);
	    fill_in_x( gas->rho, gas->T, gas->massf, M, x );
	    fill_in_DAV_im( gas->D_AB, x, DAV_im);
	    for ( int isp=0; isp<nsp; ++isp )
	    	if ( DAV_im[isp] > D_c ) D_c = DAV_im[isp];
	}
	D_c *= viscous_factor;
	// 4. Find the maximum effective diffusion
	double D_max_i = MAXIMUM( D_a, D_b );
	double D_max_ii = MAXIMUM( D_b, D_c );
	double D_max = MAXIMUM( D_max_i, D_max_ii );
	// 5. Add signal frequency contribution
	if ( dimensions == 3 ) {
	    signal += D_max * (  1.0 / (iLength * iLength)
			   	+ 1.0 / (jLength * jLength)
				+ 1.0 / (kLength * kLength));
	} else {
	    signal += D_max * (  1.0 / (east->length * east->length)
		   		+ 1.0 / (north->length * north->length));
	}
#       endif
    }
    if ( get_k_omega_flag() == 1 && !SEPARATE_UPDATE_FOR_K_OMEGA_SOURCE ) {
	if ( fs->omega > signal ) signal = fs->omega;
    }
    return signal;
} // end of signal_frequency()


int FV_Cell::turbulence_viscosity_zero(void)
{
    fs->mu_t = 0.0;
    fs->k_t = 0.0;
    return SUCCESS;
}


int FV_Cell::turbulence_viscosity_zero_if_not_in_zone(void)
{
    if ( in_turbulent_zone ) {
	/* Do nothing, leaving the turbulence quantities as set. */ ;
    } else {
	/* Presume this part of the flow is laminar; clear turbulence quantities. */
	fs->mu_t = 0.0;
	fs->k_t = 0.0;
    }
    return SUCCESS;
}


// Limit the turbulent viscosity to reasonable values relative to
// the local molecular viscosity.
// In shock started flows, we seem to get crazy values on the
// starting shock structure and the simulations do not progress.
int FV_Cell::turbulence_viscosity_limit(double factor)
{
    fs->mu_t = MINIMUM(fs->mu_t, factor * fs->gas->mu);
    fs->k_t = MINIMUM(fs->k_t , factor * fs->gas->k[0]);
    return SUCCESS;
}

// Scale the turbulent viscosity to model effects
// such as not-fully-developed turbulence that might be expected
// in short-duration transient flows.
int FV_Cell::turbulence_viscosity_factor(double factor)
{
    fs->mu_t *= factor;
    fs->k_t *= factor;
    return SUCCESS;
}

/// \brief k-omega estimate of the turbulence viscosity in the cell.
///
/// Based on Wilcox' 2006 model.
/// Jan-Pieter Nap, PJ, January 2007.
///
/// Implementation of the 3D terms
/// Wilson Chan, December 2008
int FV_Cell::turbulence_viscosity_k_omega(void)
{
    if ( get_k_omega_flag() == 0 ) {
	fs->mu_t = 0.0;
	fs->k_t = 0.0;
	return SUCCESS;
    }
    global_data &G = *get_global_data_ptr();
    double dudx, dudy, dvdx, dvdy;
    double S_bar_squared;
    double C_lim = 0.875;
    double beta_star = 0.09;
    if ( G.dimensions == 2 ) {
        // 2D cartesian or 2D axisymmetric
        dudx = 0.25 * (vtx[0]->dudx + vtx[1]->dudx + vtx[2]->dudx + vtx[3]->dudx);
        dudy = 0.25 * (vtx[0]->dudy + vtx[1]->dudy + vtx[2]->dudy + vtx[3]->dudy);
        dvdx = 0.25 * (vtx[0]->dvdx + vtx[1]->dvdx + vtx[2]->dvdx + vtx[3]->dvdx);
        dvdy = 0.25 * (vtx[0]->dvdy + vtx[1]->dvdy + vtx[2]->dvdy + vtx[3]->dvdy);
        if ( get_axisymmetric_flag() == 1 ) {
            // 2D axisymmetric
            double v_over_y = fs->vel.y / pos.y;
            S_bar_squared = dudx*dudx + dvdy*dvdy + v_over_y*v_over_y
		- 1.0/3.0 * (dudx + dvdy + v_over_y)
		* (dudx + dvdy + v_over_y)
		+ 0.5 * (dudy + dvdx) * (dudy + dvdx) ;
        } else {
            // 2D cartesian
            S_bar_squared = dudx*dudx + dvdy*dvdy
		- 1.0/3.0 * (dudx + dvdy) * (dudx + dvdy)
		+ 0.5 * (dudy + dvdx) * (dudy + dvdx);
        }
    } else {
        // 3D cartesian
        double dudz, dvdz, dwdx, dwdy, dwdz;
        dudx = 0.125 * (vtx[0]->dudx + vtx[1]->dudx + vtx[2]->dudx + vtx[3]->dudx +
                        vtx[4]->dudx + vtx[5]->dudx + vtx[6]->dudx + vtx[7]->dudx);
        dudy = 0.125 * (vtx[0]->dudy + vtx[1]->dudy + vtx[2]->dudy + vtx[3]->dudy +
                        vtx[4]->dudy + vtx[5]->dudy + vtx[6]->dudy + vtx[7]->dudy);
        dudz = 0.125 * (vtx[0]->dudz + vtx[1]->dudz + vtx[2]->dudz + vtx[3]->dudz +
                        vtx[4]->dudz + vtx[5]->dudz + vtx[6]->dudz + vtx[7]->dudz);
        dvdx = 0.125 * (vtx[0]->dvdx + vtx[1]->dvdx + vtx[2]->dvdx + vtx[3]->dvdx +
                        vtx[4]->dvdx + vtx[5]->dvdx + vtx[6]->dvdx + vtx[7]->dvdx);
        dvdy = 0.125 * (vtx[0]->dvdy + vtx[1]->dvdy + vtx[2]->dvdy + vtx[3]->dvdy +
                        vtx[4]->dvdy + vtx[5]->dvdy + vtx[6]->dvdy + vtx[7]->dvdy);
        dvdz = 0.125 * (vtx[0]->dvdz + vtx[1]->dvdz + vtx[2]->dvdz + vtx[3]->dvdz +
                        vtx[4]->dvdz + vtx[5]->dvdz + vtx[6]->dvdz + vtx[7]->dvdz);
        dwdx = 0.125 * (vtx[0]->dwdx + vtx[1]->dwdx + vtx[2]->dwdx + vtx[3]->dwdx +
                        vtx[4]->dwdx + vtx[5]->dwdx + vtx[6]->dwdx + vtx[7]->dwdx);
        dwdy = 0.125 * (vtx[0]->dwdy + vtx[1]->dwdy + vtx[2]->dwdy + vtx[3]->dwdy +
                        vtx[4]->dwdy + vtx[5]->dwdy + vtx[6]->dwdy + vtx[7]->dwdy);
        dwdz = 0.125 * (vtx[0]->dwdz + vtx[1]->dwdz + vtx[2]->dwdz + vtx[3]->dwdz +
                        vtx[4]->dwdz + vtx[5]->dwdz + vtx[6]->dwdz + vtx[7]->dwdz);
        // 3D cartesian
        S_bar_squared =  dudx*dudx + dvdy*dvdy + dwdz*dwdz
                         - 1.0/3.0*(dudx + dvdy + dwdz)*(dudx + dvdy + dwdz)
                         + 0.5 * (dudy + dvdx) * (dudy + dvdx)
                         + 0.5 * (dudz + dwdx) * (dudz + dwdx)
                         + 0.5 * (dvdz + dwdy) * (dvdz + dwdy);
    }
    S_bar_squared = MAXIMUM(0.0, S_bar_squared);
    double omega_t = MAXIMUM(fs->omega, C_lim*sqrt(2.0*S_bar_squared/beta_star));
    fs->mu_t = fs->gas->rho * fs->tke / omega_t;
    double Pr_t = get_turbulence_prandtl_number();
    Gas_model *gmodel = get_gas_model_ptr();
    int status_flag;
    fs->k_t = gmodel->Cp(*(fs->gas), status_flag) * fs->mu_t / Pr_t;
    return SUCCESS;
} // end turbulence_viscosity_k_omega()


/// \brief Update the k-omega properties within the cell over a time step.
///
/// This routine is for the source terms only and may used instead of
/// the predictor-corrector update that is applied to the rest of
/// the convective terms. It is hoped that it will give a more stable
/// integration. Note that the spatial derivatives will not be updated.
///
/// Oct 2007: Added some brute-force limiting so that the integration
///           doesn't go so crazy.
/// Nov 2008: Added feature that stops update if cell is in a laminar
///           region.
/// Jan 2011: Improve existing update scheme by using an implicit update.
///           It should make the update "inherently stable". Have kept the
///           current explicit update scheme in the code, just in case the
///           implicit update scheme fails.
///
int FV_Cell::update_k_omega_properties(double dt)
{
    // Do not update k_omega properties if we are in laminar block
    if ( in_turbulent_zone == 0 ) return SUCCESS;

#   define KOMEGA_IMPLICIT_UPDATE_FLAG 1
#   if KOMEGA_IMPLICIT_UPDATE_FLAG == 1
    double DrtkeDt_perturbTke, DromegaDt_perturbTke;
    double DrtkeDt_perturbOmega, DromegaDt_perturbOmega;
    double DGkDzetak, DGkDzetaw, DGwDzetak, DGwDzetaw;
    double DfkDk, DfkDw, DfwDk, DfwDw;
    double Gk, Gw;
    double delta_rtke, delta_romega;
    double tke, omega;
    double tke_current, omega_current;
    double tke_updated, omega_updated;
    double DrtkeDt_current, DromegaDt_current;
    double DrtkeDt_updated, DromegaDt_updated;
    double perturbFactor = 1.01;  // Perturbation factor for perturbation
                                  // analysis to get derivatives
    double tol = 1.0e-6;          // Tolerance for the Newton-solve loop

    // Encode conserved quantities for cell.
    U->tke = fs->gas->rho * fs->tke;
    U->omega = fs->gas->rho * fs->omega;

    // Start of implicit updating scheme.
    tke_current = fs->tke; omega_current = fs->omega;  // Current values of tke and omega
    tke_updated = fs->tke; omega_updated = fs->omega;  // First guess of updated values

    // Work out values of Drtke_current and DromegaDt_current.
    this->k_omega_time_derivatives(&DrtkeDt_current, &DromegaDt_current, tke_current, omega_current);

    // Implicit updating scheme.
    // A for-loop is used to limit the Newton-solve to 20 steps
    // just in case convergence does not occur.
    for ( int i = 1; i <= 20; ++i ) {
        // Work out unperturbed values of Drtke_updated and DromegaDt_updated.
        this->k_omega_time_derivatives(&DrtkeDt_updated, &DromegaDt_updated, tke_updated, omega_updated);
        // Perturb tke and obtain perturbed values of DrtkeDt_updated and DromegaDt_updated.
        tke = perturbFactor * tke_updated; omega = omega_updated;
        this->k_omega_time_derivatives(&DrtkeDt_perturbTke, &DromegaDt_perturbTke, tke, omega);
        // Perturb omega and obtain perturbed values of DrtkeDt_updated and DromegaDt_updated.
        tke = tke_updated; omega = perturbFactor * omega_updated;
        this->k_omega_time_derivatives(&DrtkeDt_perturbOmega, &DromegaDt_perturbOmega, tke, omega);
        // Compute derivatives from perturb values.
        // FIX-ME : Dividing by tke and omega (instead of rtke and romega) seems to work (gives
        //          same results as explicit update scheme), but we will keep this note here for
        //          future reference..
        DfkDk = (DrtkeDt_perturbTke - DrtkeDt_updated) / ((perturbFactor - 1.0) * tke_updated);
        DfkDw = (DrtkeDt_perturbOmega - DrtkeDt_updated) / ((perturbFactor - 1.0) * omega_updated);
        DfwDk = (DromegaDt_perturbTke - DromegaDt_updated) / ((perturbFactor - 1.0) * tke_updated);
        DfwDw = (DromegaDt_perturbOmega - DromegaDt_updated) / ((perturbFactor - 1.0) * omega_updated);
        // Compute components in matrix A of Ax=B problem.
        DGkDzetak = -1.0 + 0.5 * dt * DfkDk;
        DGkDzetaw = 0.5 * dt * DfkDw;
        DGwDzetak = 0.5 * dt * DfwDk;
        DGwDzetaw = -1.0 + 0.5 * dt * DfwDw;
        // Compute vector B of Ax=B problem.
        Gk = fs->gas->rho * tke_updated - fs->gas->rho * tke_current -
              0.5 * dt * (DrtkeDt_updated + DrtkeDt_current);
        Gw = fs->gas->rho * omega_updated - fs->gas->rho * omega_current -
              0.5 * dt * (DromegaDt_updated + DromegaDt_current);
        // Solve Ax=B algebraically.
        delta_rtke = (DGkDzetaw * Gw - DGwDzetaw * Gk) /
                      (DGwDzetak * DGkDzetaw - DGwDzetaw * DGkDzetak);
        delta_romega = (Gk - DGkDzetak * delta_rtke) / DGkDzetaw;
        // Assign the updated tke and omega values if delta_rtke and
        // delta_romega are both smaller than the given tolerance, and
        // then break out from the Newton-solve loop.
        if (fabs(delta_rtke) <= tol && fabs(delta_romega) <= tol) {
            fs->tke = tke_updated;
            fs->omega = omega_updated;
            break;
        } else {
            // Compute next estimates for rtke and romega from
            // delta_rtke and delta_romega.
            if (delta_rtke + fs->gas->rho * tke_updated < 0.0) {
                // Don't let rtke go negative.
                U->tke = fs->gas->rho * tke_updated;
            } else {
                // Next estimate for rtke.
                U->tke = delta_rtke + fs->gas->rho * tke_updated;
            }
            if (delta_romega + fs->gas->rho * omega_updated < 0.0) {
                // Don't let romega go negative.
                U->omega = fs->gas->rho * omega_updated;
            } else {
                // Next estimate for romega.
                U->omega = delta_romega + fs->gas->rho * omega_updated;
            }
            // Decode for the next step of the Newton-solve loop
            tke_updated = U->tke / fs->gas->rho;
            omega_updated = U->omega / fs->gas->rho;
        }
    }  // End of Newton-solve loop for implicit update scheme

#   else
    // Explicit updating scheme.
    double DrtkeDt, DromegaDt, rtke_increment, romega_increment;
    int n_little_steps = 20;  // Make this as large as needed to get stability.
    double dt_little = dt / n_little_steps;
    
    // Encode conserved quantities for cell.
    U->tke = fs->gas->rho * fs->tke;
    U->omega = fs->gas->rho * fs->omega;
    for ( int i = 1; i <= n_little_steps; ++i ) {
        this->k_omega_time_derivatives(&DrtkeDt, &DromegaDt, fs->tke, fs->omega);
        rtke_increment = dt_little * DrtkeDt;
        romega_increment = dt_little * DromegaDt;
        if ( U->tke + rtke_increment < 0.0 ||
             (rtke_increment > 0.0 && U->tke + rtke_increment > 0.5 * U->total_energy) ) {
            // Don't let rtke go negative and don't let it grow too large.
            rtke_increment = 0.0;
	}
	if ( U->omega + romega_increment < 0.0 ) {
	    // Don't let romega go negative.
            romega_increment = 0.0;
            }
	U->tke += rtke_increment;
	U->omega += romega_increment;
	// Decode conserved quantities.
	fs->tke = U->tke / fs->gas->rho;
	fs->omega = U->omega / fs->gas->rho;
    }  // End of for-loop for explicit update scheme
#   endif

    return SUCCESS;
} // end update_k_omega_properties()


/// \brief Compute k-omega source terms.
///
/// Production and Dissipation expressions for turbulence kinetic energy
/// and turbulence frequency (or pseudo-vorticity). Based on Wilcox's 2006 model.
///
/// Jan 2007: Initial implementation (Jan-Pieter Nap, PJ)
/// Dec 2008: Implementation of the 3D terms (W Chan)
/// Jan 2011: Minor modification to allow for implicit updating of tke and omega (W Chan)
///           All "fs->tke" and "fs->omega" instances are replaced with tke and omega.
///
int FV_Cell::k_omega_time_derivatives(double *Q_rtke, double *Q_romega, double tke, double omega)
{
    if ( get_k_omega_flag() == 0 ) {
	*Q_rtke = 0.0;
	*Q_romega = 0.0;
	return SUCCESS;
    }
    global_data &G = *get_global_data_ptr();
    double dudx, dudy, dvdx, dvdy;
    double dtkedx, dtkedy, domegadx, domegady;
    double alpha = 0.52;
    double beta_0 = 0.0708;
    double beta;
    double beta_star = 0.09;
    double P_K, D_K, P_W, D_W;
    double cross_diff;
    double sigma_d = 0.0;
    double WWS, X_w, f_beta;
    if ( G.dimensions == 2 ) {
        // 2D cartesian or 2D axisymmetric
        dudx = 0.25 * (vtx[0]->dudx + vtx[1]->dudx + vtx[2]->dudx + vtx[3]->dudx);
        dudy = 0.25 * (vtx[0]->dudy + vtx[1]->dudy + vtx[2]->dudy + vtx[3]->dudy);
        dvdx = 0.25 * (vtx[0]->dvdx + vtx[1]->dvdx + vtx[2]->dvdx + vtx[3]->dvdx);
        dvdy = 0.25 * (vtx[0]->dvdy + vtx[1]->dvdy + vtx[2]->dvdy + vtx[3]->dvdy);
        dtkedx = 0.25 * (vtx[0]->dtkedx + vtx[1]->dtkedx + vtx[2]->dtkedx + vtx[3]->dtkedx);
        dtkedy = 0.25 * (vtx[0]->dtkedy + vtx[1]->dtkedy + vtx[2]->dtkedy + vtx[3]->dtkedy);
        domegadx = 0.25 * (vtx[0]->domegadx + vtx[1]->domegadx + vtx[2]->domegadx + vtx[3]->domegadx);
        domegady = 0.25 * (vtx[0]->domegady + vtx[1]->domegady + vtx[2]->domegady + vtx[3]->domegady);
        if ( get_axisymmetric_flag() == 1 ) {
            // 2D axisymmetric
            double v_over_y = fs->vel.y / pos.y;
            // JP.Nap correction from 03-May-2007 (-v_over_y in parentheses)
            // P_K -= 0.6667 * mu_t * v_over_y * (dudx+dvdy-v_over_y);
            // Wilson Chan correction to JP Nap's version (13 Dec 2008)
            P_K = 2.0 * fs->mu_t * (dudx*dudx + dvdy*dvdy)
                  + fs->mu_t * (dudy + dvdx) * (dudy + dvdx)
                  - 2.0/3.0 * fs->mu_t * (dudx + dvdy + v_over_y)
                  * (dudx + dvdy + v_over_y)
                  + 2.0 * fs->mu_t * (v_over_y) * (v_over_y)
                  - 2.0/3.0 * fs->gas->rho * tke * (dudx + dvdy + v_over_y);
            WWS = 0.25 * (dvdx - dudy) * (dvdx - dudy) * v_over_y ;
        } else {
            // 2D cartesian
            P_K = 1.3333 * fs->mu_t * (dudx*dudx - dudx*dvdy + dvdy*dvdy)
                  + fs->mu_t * (dudy + dvdx) * (dudy + dvdx)
                  - 0.66667 * fs->gas->rho * tke * (dudx + dvdy);
            WWS = 0.0 ;
        }
        cross_diff = dtkedx * domegadx + dtkedy * domegady ;
    } else {
        // 3D cartesian
        double dudz, dvdz, dwdx, dwdy, dwdz;
        double dtkedz, domegadz;
        dudx = 0.125 * (vtx[0]->dudx + vtx[1]->dudx + vtx[2]->dudx + vtx[3]->dudx +
                        vtx[4]->dudx + vtx[5]->dudx + vtx[6]->dudx + vtx[7]->dudx);
        dudy = 0.125 * (vtx[0]->dudy + vtx[1]->dudy + vtx[2]->dudy + vtx[3]->dudy +
                        vtx[4]->dudy + vtx[5]->dudy + vtx[6]->dudy + vtx[7]->dudy);
        dudz = 0.125 * (vtx[0]->dudz + vtx[1]->dudz + vtx[2]->dudz + vtx[3]->dudz +
                        vtx[4]->dudz + vtx[5]->dudz + vtx[6]->dudz + vtx[7]->dudz);
        dvdx = 0.125 * (vtx[0]->dvdx + vtx[1]->dvdx + vtx[2]->dvdx + vtx[3]->dvdx +
                        vtx[4]->dvdx + vtx[5]->dvdx + vtx[6]->dvdx + vtx[7]->dvdx);
        dvdy = 0.125 * (vtx[0]->dvdy + vtx[1]->dvdy + vtx[2]->dvdy + vtx[3]->dvdy +
                        vtx[4]->dvdy + vtx[5]->dvdy + vtx[6]->dvdy + vtx[7]->dvdy);
        dvdz = 0.125 * (vtx[0]->dvdz + vtx[1]->dvdz + vtx[2]->dvdz + vtx[3]->dvdz +
                        vtx[4]->dvdz + vtx[5]->dvdz + vtx[6]->dvdz + vtx[7]->dvdz);
        dwdx = 0.125 * (vtx[0]->dwdx + vtx[1]->dwdx + vtx[2]->dwdx + vtx[3]->dwdx +
                        vtx[4]->dwdx + vtx[5]->dwdx + vtx[6]->dwdx + vtx[7]->dwdx);
        dwdy = 0.125 * (vtx[0]->dwdy + vtx[1]->dwdy + vtx[2]->dwdy + vtx[3]->dwdy +
                        vtx[4]->dwdy + vtx[5]->dwdy + vtx[6]->dwdy + vtx[7]->dwdy);
        dwdz = 0.125 * (vtx[0]->dwdz + vtx[1]->dwdz + vtx[2]->dwdz + vtx[3]->dwdz +
                        vtx[4]->dwdz + vtx[5]->dwdz + vtx[6]->dwdz + vtx[7]->dwdz);
        dtkedx = 0.125 * (vtx[0]->dtkedx + vtx[1]->dtkedx + vtx[2]->dtkedx + vtx[3]->dtkedx +
                          vtx[4]->dtkedx + vtx[5]->dtkedx + vtx[6]->dtkedx + vtx[7]->dtkedx);
        dtkedy = 0.125 * (vtx[0]->dtkedy + vtx[1]->dtkedy + vtx[2]->dtkedy + vtx[3]->dtkedy +
                          vtx[4]->dtkedy + vtx[5]->dtkedy + vtx[6]->dtkedy + vtx[7]->dtkedy);
        dtkedz = 0.125 * (vtx[0]->dtkedz + vtx[1]->dtkedz + vtx[2]->dtkedz + vtx[3]->dtkedz +
                          vtx[4]->dtkedz + vtx[5]->dtkedz + vtx[6]->dtkedz + vtx[7]->dtkedz);
        domegadx = 0.125 * (vtx[0]->domegadx + vtx[1]->domegadx + vtx[2]->domegadx + vtx[3]->domegadx +
                            vtx[4]->domegadx + vtx[5]->domegadx + vtx[6]->domegadx + vtx[7]->domegadx);
        domegady = 0.125 * (vtx[0]->domegady + vtx[1]->domegady + vtx[2]->domegady + vtx[3]->domegady +
                            vtx[4]->domegady + vtx[5]->domegady + vtx[6]->domegady + vtx[7]->domegady);
        domegadz = 0.125 * (vtx[0]->domegadz + vtx[1]->domegadz + vtx[2]->domegadz + vtx[3]->domegadz +
                            vtx[4]->domegadz + vtx[5]->domegadz + vtx[6]->domegadz + vtx[7]->domegadz);
        P_K = 2.0 * fs->mu_t * (dudx*dudx + dvdy*dvdy + dwdz*dwdz)
              - 2.0/3.0 * fs->mu_t * (dudx + dvdy + dwdz) * (dudx + dvdy + dwdz)
              - 2.0/3.0 * fs->gas->rho * tke * (dudx + dvdy + dwdz)
              + fs->mu_t * (dudy + dvdx) * (dudy + dvdx)
              + fs->mu_t * (dudz + dwdx) * (dudz + dwdx)
              + fs->mu_t * (dvdz + dwdy) * (dvdz + dwdy) ;
        cross_diff = dtkedx * domegadx + dtkedy * domegady + dtkedz * domegadz ;
        WWS = 0.25 * (dudy - dvdx) * (dudy - dvdx) * dwdz
              + 0.25 * (dudz - dwdx) * (dudz - dwdx) * dvdy
              + 0.25 * (dvdz - dwdy) * (dvdz - dwdy) * dudx
              + 0.25 * (dudy - dvdx) * (dvdz - dwdy) * (dwdx + dudz)
              + 0.25 * (dudz - dwdx) * (dwdy - dvdz) * (dudy + dvdx)
              + 0.25 * (dvdx - dudy) * (dudz - dwdx) * (dwdy + dvdx) ;
    }

    D_K = beta_star * fs->gas->rho * tke * omega;
    
    // Apply a limit to the tke production as suggested by Jeff White, November 2007.
#   define P_OVER_D_LIMIT (25.0)
    P_K = MINIMUM(P_K, P_OVER_D_LIMIT*D_K);

    if ( cross_diff > 0 ) sigma_d = 0.125;
    P_W = alpha * omega / MAXIMUM(tke,SMALL_TKE) * P_K +
          sigma_d * fs->gas->rho / MAXIMUM(omega,SMALL_OMEGA) * cross_diff;

    X_w = FABS(WWS / pow(beta_star*omega, 3)) ;
    f_beta = (1.0 + 85.0 * X_w) / (1.0 + 100.0 * X_w) ;
    beta = beta_0 * f_beta;
    D_W = beta * fs->gas->rho * omega * omega;

    *Q_rtke = P_K - D_K;
    *Q_romega = P_W - D_W;
    return SUCCESS;
} // end k_omega_time_derivatives()


/// \brief Compute the components of the source vector, Q, for inviscid flow.
///
/// Currently, the axisymmetric equations include the
/// pressure contribution to the y-momentum equation
/// here rather than in the boundary fluxes.
int FV_Cell::inviscid_source_vector(int time_level, double omegaz)
{
    // By default, assume 2D-planar, or 3D-Cartesian flow.
    Q->mass = 0.0;
    Q->momentum.x = 0.0;
    Q->momentum.y = 0.0;
    Q->momentum.z = 0.0;
    // Magnetic field -- FIX-ME -- Daryl
    Q->B.x = 0.0;
    Q->B.y = 0.0;
    Q->B.z = 0.0;
    Q->total_energy = get_heat_factor() * base_qdot;
    Q->tke = 0.0;
    Q->omega = 0.0;
    if ( omegaz != 0.0 ) {
	// Rotating frame.
	double rho = fs->gas->rho;
	double x = pos.x;
	double y = pos.y;
        double wx = fs->vel.x;
	double wy = fs->vel.y;
	// Coriolis and centrifugal forces contribute to momenta.
	Q->momentum.x += rho * (omegaz*omegaz*x + 2.0*omegaz*wy);
	Q->momentum.y += rho * (omegaz*omegaz*y - 2.0*omegaz*wx);
	// There is no contribution to the energy equation in the rotating frame
        // because it is implicit in the use of rothalpy as the conserved quantity.
    }
    if ( get_axisymmetric_flag() == 1 ) {
	// For axisymmetric flow:
	// pressure contribution from the Front and Back (radial) interfaces.
	Q->momentum.y += fs->gas->p * ar[time_level] / vol[time_level];
    }
    // Species production (other than chemistry).
    // For the chemistry, see chemical_increment().
    for ( size_t isp = 0; isp < Q->massf.size(); ++isp ) Q->massf[isp] = 0.0;
    // Individual energies (other than energy exchange)
    // For the energy exchange, see thermal_increment()
    for ( size_t imode = 0; imode < Q->energies.size(); ++imode )
	Q->energies[imode] = 0.0;
    // Radiation can potentially be removed from both the electronic and
    // total energy source terms.
    if ( get_radiation_flag() == 1 ) {
	// Radiative source term should be already calculated
	// Add value to total energy
	// FIX-ME: - assuming electronic mode is the last in the vector of energies
	//         - what about Q_renergies[0]?
	Q->total_energy += Q_rE_rad;
	Q->energies.back() += Q_rE_rad;
    } else {
	// No radiation is being considered.
	Q_rE_rad = 0.0;
    }
    return SUCCESS;
} // end inviscid_source_vector()


/// \brief Compute the components of the source vector, Q, for viscous flow.
int FV_Cell::viscous_source_vector(void)
{
    double dudx, dvdy;
    double viscous_factor;
    double mu, lmbda, tau_00;
    double v_over_y;

    viscous_factor = get_viscous_factor();
    if ( get_axisymmetric_flag() == 1 ) {
	v_over_y = fs->vel.y / pos.y;
    } else {
	v_over_y = 0.0;
    }

    // By default, assume 2D-planar, or 3D-Cartesian flow.
    Q->mass = 0.0;
    Q->momentum.x = 0.0;
    Q->momentum.y = 0.0;
    Q->momentum.z = 0.0;
    // Magnetic field -- FIX-ME -- Daryl
    Q->B.x = 0.0;
    Q->B.y = 0.0;
    Q->B.z = 0.0;
    Q->total_energy = 0.0;

    for ( size_t isp = 0; isp < Q->massf.size(); ++isp ) Q->massf[isp] = 0.0;
    for ( size_t imode = 0; imode < Q->energies.size(); ++imode )
	Q->energies[imode] = 0.0;

    dudx = 0.25 * (vtx[0]->dudx + vtx[1]->dudx + vtx[2]->dudx + vtx[3]->dudx);
    dvdy = 0.25 * (vtx[0]->dvdy + vtx[1]->dvdy + vtx[2]->dvdy + vtx[3]->dvdy);

    if ( get_k_omega_flag() == 1 && !SEPARATE_UPDATE_FOR_K_OMEGA_SOURCE ) {
	this->k_omega_time_derivatives(&(Q->tke), &(Q->omega), fs->tke, fs->omega);
    } else {
	Q->tke = 0.0;
	Q->omega = 0.0;
    }

    if ( get_axisymmetric_flag() == 1 ) {
	// For viscous, axisymmetric flow:
	mu = 0.25 * (iface[EAST]->fs->gas->mu + iface[WEST]->fs->gas->mu +
		     iface[NORTH]->fs->gas->mu + iface[SOUTH]->fs->gas->mu) +
	    0.25 * (iface[EAST]->fs->mu_t + iface[WEST]->fs->mu_t +
		    iface[NORTH]->fs->mu_t + iface[SOUTH]->fs->mu_t);
	mu *= viscous_factor;
	lmbda = -2.0/3.0 * mu;
	tau_00 = 2.0 * mu * v_over_y + lmbda * (dudx + dvdy + v_over_y);
	// Y-Momentum; viscous stress contribution from the front and Back interfaces.
	// Note that these quantities are approximated at the
	// mid-point of the cell face and so should never be
	// singular -- at least I hope that this is so.
	Q->momentum.y -= tau_00 * area/ volume;
    }
    return SUCCESS;
} // end viscous_source_vector()


/// \brief Calculate the Reynolds number at a wall interface
double FV_Cell::calculate_wall_Reynolds_number(int which_boundary)
{
    FV_Interface * IFace = iface[which_boundary];
    double Re_wall, a_wall, cell_width = 0.0;
    Gas_model * gm = get_gas_model_ptr();
    gm->eval_thermo_state_rhoT(*(IFace->fs->gas));
    a_wall = IFace->fs->gas->a;
    if ( which_boundary==EAST || which_boundary==WEST )
	cell_width = iLength;
    else if ( which_boundary==NORTH || which_boundary==SOUTH )
	cell_width = jLength;
    else if ( which_boundary==TOP || which_boundary==BOTTOM )
	cell_width = kLength;
    Re_wall = IFace->fs->gas->rho * a_wall * cell_width / IFace->fs->gas->mu;

    return Re_wall;
}


/// \brief Store parameters for (re-)scaling of radiative source term
///
/// Simple rho x T**4 scaling seems to be adequate
int FV_Cell::store_rad_scaling_params(void)
{
    // 1. Store the freshly computed radiative flux as the 'original'
    Q_rad_org = Q_rE_rad;
    // 2. Compute the scaling factor based on local gas properties
    // NOTE: - The idea is that f_rad_org is proportional to actual value
    //       - Assuming that the last temperature is the electronic temperature
    double T = fs->gas->T.back();
    if ( Q_rad_org <= 0.0 ) {
	// This cell is a net emitter
        f_rad_org = fs->gas->rho * pow(T, 4);
    }
    else if ( Q_rad_org > 0.0 ) {
	// This cell is a net absorber
	f_rad_org = fs->gas->rho / pow(T, 4);
    }
    return SUCCESS;
}


/// \brief (Re-)scale radiative source term for loose coupling
int FV_Cell::rescale_Q_rE_rad(void)
{
    // 1. Compute the current scaling factor based on local gas properties
    double T = fs->gas->T[0];
    double f_rad_new = 1.0;
    if ( Q_rad_org <= 0.0 ) {
	// This cell is a net emitter
        f_rad_new = fs->gas->rho * pow(T, 4);
    }
    else if ( Q_rad_org > 0.0 ) {
	// This cell is a net absorber
	f_rad_new = fs->gas->rho / pow(T, 4);
    }
    // 2. (Re-)scale the original source term
    Q_rE_rad = ( f_rad_new / f_rad_org ) * Q_rad_org;
    return SUCCESS;
}


/// \brief (Re-)set radiative source terms to zero
int FV_Cell::reset_Q_rad_to_zero(void)
{
    Q_rE_rad = 0.0;
    return SUCCESS;
}


/// \brief Calculate fabs( f_rad - f_rad_org ) / f_rad_org
double FV_Cell::rad_scaling_ratio(void)
{
    // 1. Compute the current scaling factor based on local gas properties
    double T = fs->gas->T[0];
    double f_rad = 1.0;
    if ( Q_rE_rad <= 0.0 ) {
	// This cell is a net emitter
        f_rad = fs->gas->rho * pow(T, 4);
    }
    else if ( Q_rE_rad > 0.0 ) {
	// This cell is a net absorber
	f_rad = fs->gas->rho / pow(T, 4);
    }
    return fabs( f_rad - f_rad_org ) / f_rad_org;
}

// --------------------------------------------------------------------

int number_of_values_in_cell_copy(int type_of_copy)
{
    int number = 0;
    Gas_model *gmodel = get_gas_model_ptr();
    if (type_of_copy == COPY_ALL_CELL_DATA ||
	type_of_copy == COPY_FLOW_STATE) {
        number += gmodel->number_of_values_in_gas_data_copy();
	number += 8 + 4; // FlowState + cell data
	if ( get_mhd_flag() == 1 ) number += 3;
    }
    if (type_of_copy == COPY_ALL_CELL_DATA ||
        type_of_copy == COPY_CELL_LENGTHS) {
        number += 6;
	if ( get_shock_fitting_flag() ) {
	    number += 24; // Velocity and position for each interface (2D)	
	}
    }
    return number;
}


std::string variable_list_for_cell( void )
{
    // This function needs to be kept consistent with functions
    // FV_Cell::write_values_to_string, FV_Cell::scan_values_from_string
    // (found above) and with the corresponding Python functions
    // write_cell_data and variable_list_for_cell
    // that may be found in app/eilmer3/source/e3_flow.py.
    ostringstream ost;
    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t nmodes = gmodel->get_number_of_modes();
    ost << "\"pos.x\" \"pos.y\" \"pos.z\" \"volume\"";
    ost << " \"rho\" \"vel.x\" \"vel.y\" \"vel.z\" ";
    if ( get_mhd_flag() == 1 ) {
	ost << " \"B.x\" \"B.y\" \"B.z\" ";
    }
    ost << " \"p\" \"a\" \"mu\" \"k[0]\" \"mu_t\" \"k_t\" \"S\"";
    if ( get_radiation_flag() == 1 ) {
	ost << " \"Q_rad_org\" \"f_rad_org\" \"Q_rE_rad\"";
    }
    ost << " \"tke\" \"omega\"";
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	std::string specname = gmodel->species_name(isp);
	size_t found = specname.find(" ");
	while ( found != std::string::npos ) {
	    specname.replace(found, 1, "-");
	    found = specname.find(" ");
	}
	ost << " \"massf[" << isp << "]-" << specname << "\"";
    }
    if ( nsp > 1 ) ost << " \"dt_chem\"";
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	ost << " \"e[" << imode << "]\" \"T[" << imode << "]\"";
    }
    if ( nmodes > 1 ) ost << " \"dt_therm\"";
    return ost.str();
} // end variable_list_for_cell()


/// \brief Reconstruct flow properties at an interface from FV_Cell properties.
///
/// This is essentially a one-dimensional interpolation process.  It needs only
/// the cell-average data and the lengths of the cells in the interpolation direction.
int one_d_interp(FV_Cell &cL1, FV_Cell &cL0,
		 FV_Cell &cR0, FV_Cell &cR1,
		 double cL1Length, double cL0Length,
		 double cR0Length, double cR1Length,
		 FlowState &Lft, FlowState &Rght )
{
    Gas_model *gmodel = get_gas_model_ptr();
    Thermo_interpolator* ti = get_thermo_interpolator_ptr();
    int nsp = gmodel->get_number_of_species();
    int apply_limiter_flag = get_apply_limiter_flag();
    int extrema_clipping_flag = get_extrema_clipping_flag();

    // Low-order reconstruction just copies data from adjacent FV_Cell.
    Lft.copy_values_from(*(cL0.fs));
    Rght.copy_values_from(*(cR0.fs));
    if ( get_Xorder_flag() > 1 ) {
	// High-order reconstruction for some properties.
	one_d_interp_scalar(cL1.fs->vel.x, cL0.fs->vel.x, cR0.fs->vel.x, cR1.fs->vel.x,
			    cL1Length, cL0Length, cR0Length, cR1Length,
			    Lft.vel.x, Rght.vel.x, apply_limiter_flag, extrema_clipping_flag);
	one_d_interp_scalar(cL1.fs->vel.y, cL0.fs->vel.y, cR0.fs->vel.y, cR1.fs->vel.y,
			    cL1Length, cL0Length, cR0Length, cR1Length,
			    Lft.vel.y, Rght.vel.y, apply_limiter_flag, extrema_clipping_flag);
	one_d_interp_scalar(cL1.fs->vel.z, cL0.fs->vel.z, cR0.fs->vel.z, cR1.fs->vel.z,
			    cL1Length, cL0Length, cR0Length, cR1Length,
			    Lft.vel.z, Rght.vel.z, apply_limiter_flag, extrema_clipping_flag);
	if ( get_mhd_flag() == 1 ) {
	    one_d_interp_scalar(cL1.fs->B.x, cL0.fs->B.x, cR0.fs->B.x, cR1.fs->B.x,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.B.x, Rght.B.x, apply_limiter_flag, extrema_clipping_flag);
	    one_d_interp_scalar(cL1.fs->B.y, cL0.fs->B.y, cR0.fs->B.y, cR1.fs->B.y,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.B.y, Rght.B.y, apply_limiter_flag, extrema_clipping_flag);
	    one_d_interp_scalar(cL1.fs->B.z, cL0.fs->B.z, cR0.fs->B.z, cR1.fs->B.z,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.B.z, Rght.B.z, apply_limiter_flag, extrema_clipping_flag);
	}
	if ( get_k_omega_flag() == 1 ) {
	    one_d_interp_scalar(cL1.fs->tke, cL0.fs->tke, cR0.fs->tke, cR1.fs->tke,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.tke, Rght.tke, apply_limiter_flag, extrema_clipping_flag);
	    one_d_interp_scalar(cL1.fs->omega, cL0.fs->omega, cR0.fs->omega, cR1.fs->omega,
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.omega, Rght.omega, apply_limiter_flag, extrema_clipping_flag);
	}
        for ( int isp = 0; isp < nsp; ++isp ) {
	    one_d_interp_scalar(cL1.fs->gas->massf[isp], cL0.fs->gas->massf[isp],
				cR0.fs->gas->massf[isp], cR1.fs->gas->massf[isp],
				cL1Length, cL0Length, cR0Length, cR1Length,
				Lft.gas->massf[isp], Rght.gas->massf[isp],
				apply_limiter_flag, extrema_clipping_flag);
        }

	// Make the thermodynamic properties consistent.
	// Pressure, Local Speed of Sound and Temperature.
        // The value of 1 indicates that the old temperature
        // should be used as an initial guess for the iterative
        // EOS functions.
	// If the EOS call fouls up, just copy the cell data, low-order.
	if ( nsp > 1 ) {
	    if ( scale_mass_fractions( Lft.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp<Lft.gas->massf.size(); ++isp )
		    cL0.fs->gas->massf[isp] = Lft.gas->massf[isp];
	    }
	    if ( scale_mass_fractions( Rght.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp<Rght.gas->massf.size(); ++isp )
		    cR0.fs->gas->massf[isp] = Rght.gas->massf[isp];
	    }
	}

	// Interpolate on two of the thermodynamic quantities, and fill
	// in the rest based on an EOS call.

	int status = ti->one_d_interp(*(cL1.fs->gas), *(cL0.fs->gas), *(cR0.fs->gas), *(cR1.fs->gas),
				      cL1Length, cL0Length, cR0Length, cR1Length,
				      *(Lft.gas), *(Rght.gas));
	if ( status != SUCCESS ) {
	    if ( status == 1 ) {
		// Lft state failed.
		Lft.copy_values_from(*(cL0.fs));
	    }
	    else if ( status == 2 ) {
		// Rght state failed.
		Rght.copy_values_from(*(cR0.fs));
	    }
	    else if ( status == 3 ) {
		// Both failed.
		Lft.copy_values_from(*(cL0.fs));
		Rght.copy_values_from(*(cR0.fs));
	    }
	    else {
		printf("one_d_interp(): Problem in flow state reconstruction.");
		printf("Failure status: %d is unknown.\n", status);
		printf("Bailing out!\n");
		exit(RECONSTRUCTION_ERROR);
	    }
	}

	if ( get_viscous_flag() ) {
	    gmodel->eval_transport_coefficients(*(Lft.gas));
	    gmodel->eval_transport_coefficients(*(Rght.gas));
	}
	if ( get_diffusion_flag() ) {
	    gmodel->eval_diffusion_coefficients(*(Lft.gas));
	    gmodel->eval_diffusion_coefficients(*(Rght.gas));
	}
    } // end of high-order reconstruction
    return SUCCESS;
} // end of one_d_interp()

/// \brief Reconstruct upstream weighted flow properties at an interface from FV_Cell properties.
///
/// This is essentially a one-dimensional interpolation process.  It needs only
/// the cell-average data and the lengths of the cells in the interpolation direction.
int wone_d_interp(FV_Cell &cL1, FV_Cell &cL0,
				 FV_Cell &cR0, FV_Cell &cR1,
				 double cL1Length, double cL0Length,
				 double cR0Length, double cR1Length,
				 FlowState &Lft, FlowState &Rght )
{
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data &gL1 = *(cL1.fs->gas);
    Gas_data &gL0 = *(cL0.fs->gas);
    Gas_data &gR0 = *(cR0.fs->gas);
    Gas_data &gR1 = *(cR1.fs->gas);
    Thermo_interpolator* ti = get_thermo_interpolator_ptr();
    int nsp = gmodel->get_number_of_species();
    int apply_limiter_flag = get_apply_limiter_flag();
    int extrema_clipping_flag = get_extrema_clipping_flag();
    
	double ML = dot(cL0.fs->vel, unit(cR0.pos - cL0.pos)) / cL0.fs->gas->a;
	double MR = dot(cR0.fs->vel, unit(cL0.pos - cR0.pos)) / cR0.fs->gas->a;
	double kL;
	double kR;
	if ( ML > 1.0 ) {
	    kL = -1.0;
	} else {
		kL = -max(0.0, ML);
	}
	if ( MR > 1.0 ) {
	    kR = -1.0;
	} else {
	    kR = -max(0.0, MR);
	}
	// Low-order reconstruction just copies data from adjacent FV_Cell.
	Lft.copy_values_from(*(cL0.fs));
	Rght.copy_values_from(*(cR0.fs));
	if ( get_Xorder_flag() > 1 ) {
	// High-order reconstruction for some properties.
	    wone_d_interp_scalar(cL1.fs->vel.x, cL0.fs->vel.x, cR0.fs->vel.x, cR1.fs->vel.x,
				 cL1Length, cL0Length, cR0Length, cR1Length,
				 Lft.vel.x, Rght.vel.x, kL, kR, apply_limiter_flag, extrema_clipping_flag);
	    wone_d_interp_scalar(cL1.fs->vel.y, cL0.fs->vel.y, cR0.fs->vel.y, cR1.fs->vel.y,
				 cL1Length, cL0Length, cR0Length, cR1Length,
				 Lft.vel.y, Rght.vel.y, kL, kR, apply_limiter_flag, extrema_clipping_flag);
	    wone_d_interp_scalar(cL1.fs->vel.z, cL0.fs->vel.z, cR0.fs->vel.z, cR1.fs->vel.z,
				 cL1Length, cL0Length, cR0Length, cR1Length,
				 Lft.vel.z, Rght.vel.z, kL, kR, apply_limiter_flag, extrema_clipping_flag);
	    if ( get_mhd_flag() == 1 ) {
		wone_d_interp_scalar(cL1.fs->B.x, cL0.fs->B.x, cR0.fs->B.x, cR1.fs->B.x,
				     cL1Length, cL0Length, cR0Length, cR1Length,
				     Lft.B.x, Rght.B.x, kL, kR, apply_limiter_flag, extrema_clipping_flag);
		wone_d_interp_scalar(cL1.fs->B.y, cL0.fs->B.y, cR0.fs->B.y, cR1.fs->B.y,
				     cL1Length, cL0Length, cR0Length, cR1Length,
				     Lft.B.y, Rght.B.y, kL, kR, apply_limiter_flag, extrema_clipping_flag);
		wone_d_interp_scalar(cL1.fs->B.z, cL0.fs->B.z, cR0.fs->B.z, cR1.fs->B.z,
				     cL1Length, cL0Length, cR0Length, cR1Length,
				     Lft.B.z, Rght.B.z, kL, kR, apply_limiter_flag, extrema_clipping_flag);
	    }
	    if ( get_k_omega_flag() == 1 ) {
		wone_d_interp_scalar(cL1.fs->tke, cL0.fs->tke, cR0.fs->tke, cR1.fs->tke,
				     cL1Length, cL0Length, cR0Length, cR1Length,
				     Lft.tke, Rght.tke, kL, kR, apply_limiter_flag, extrema_clipping_flag);
		wone_d_interp_scalar(cL1.fs->omega, cL0.fs->omega, cR0.fs->omega, cR1.fs->omega,
				     cL1Length, cL0Length, cR0Length, cR1Length,
				     Lft.omega, Rght.omega, kL, kR, apply_limiter_flag, extrema_clipping_flag);
	    }
	    for ( int isp = 0; isp < nsp; ++isp ) {
		wone_d_interp_scalar(cL1.fs->gas->massf[isp], cL0.fs->gas->massf[isp],
				     cR0.fs->gas->massf[isp], cR1.fs->gas->massf[isp],
				     cL1Length, cL0Length, cR0Length, cR1Length,
				     Lft.gas->massf[isp], Rght.gas->massf[isp],  kL, kR,
				     apply_limiter_flag, extrema_clipping_flag);
	    }
	    
	    // Make the thermodynamic properties consistent.
	    // Pressure, Local Speed of Sound and Temperature.
	    // The value of 1 indicates that the old temperature
	    // should be used as an initial guess for the iterative
	    // EOS functions.
	    // If the EOS call fouls up, just copy the cell data, low-order.
	    if ( nsp > 1 ) {
		if ( scale_mass_fractions( Lft.gas->massf ) != 0 ) {
		    for ( size_t isp=0; isp<Lft.gas->massf.size(); ++isp )
			cL0.fs->gas->massf[isp] = Lft.gas->massf[isp];
		}
		if ( scale_mass_fractions( Rght.gas->massf ) != 0 ) {
		    for ( size_t isp=0; isp<Rght.gas->massf.size(); ++isp )
			cR0.fs->gas->massf[isp] = Rght.gas->massf[isp];
		}
	    }
	    
	    // Interpolate on two of the thermodynamic quantities, and fill
	    // in the rest based on an EOS call.
	    
	    wone_d_interp_scalar(gL1.rho, gL0.rho, gR0.rho, gR1.rho,
				 cL1Length, cL0Length, cR0Length, cR1Length,
				 Lft.gas->rho, Rght.gas->rho,
				 kL, kR,
				 apply_limiter_flag, extrema_clipping_flag);
	    for ( int i = 0; i < gmodel->get_number_of_modes(); ++i ) {
		wone_d_interp_scalar(gL1.e[i], gL0.e[i], gR0.e[i], gR1.e[i],
				     cL1Length, cL0Length, cR0Length, cR1Length,
				     Lft.gas->e[i], Rght.gas->e[i],
				     kL, kR,
				     apply_limiter_flag, extrema_clipping_flag);
	    }
	    int status = gmodel->eval_thermo_state_rhoe(*(Lft.gas));
	    status = gmodel->eval_thermo_state_rhoe(*(Rght.gas));
	    if ( status != SUCCESS ) {
		if ( status == 1 ) {
		    // Lft state failed.
		    Lft.copy_values_from(*(cL0.fs));
		}
		else if ( status == 2 ) {
		    // Rght state failed.
		    Rght.copy_values_from(*(cR0.fs));
		}
		else if ( status == 3 ) {
		    // Both failed.
		    Lft.copy_values_from(*(cL0.fs));
		    Rght.copy_values_from(*(cR0.fs));
		}
		else {
		    printf("one_d_interp(): Problem in flow state reconstruction.");
		    printf("Failure status: %d is unknown.\n", status);
		    printf("Bailing out!\n");
		    exit(RECONSTRUCTION_ERROR);
		}
	    }
	    
	    if ( get_viscous_flag() ) {
		gmodel->eval_transport_coefficients(*(Lft.gas));
		gmodel->eval_transport_coefficients(*(Rght.gas));
	    }
	    if ( get_diffusion_flag() ) {
		gmodel->eval_diffusion_coefficients(*(Lft.gas));
		gmodel->eval_diffusion_coefficients(*(Rght.gas));
	    }
	} // end of high-order reconstruction
	return SUCCESS;
} // end of wone_d_interp()

/// \brief Reconstruct flow properties at an interface from FV_Cell properties from one side.
///
/// This is essentially a one-sided one-dimensional interpolation process.  It needs only
/// the cell-average data and the lengths of the cells in the interpolation direction.
int onesided_interp(FV_Cell &cL0, FV_Cell &cR0, FV_Cell &cR1,
		    double cL0Length, double cR0Length, double cR1Length,
		    FlowState &Rght )
{
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data &gL0 = *(cL0.fs->gas);
    Gas_data &gR0 = *(cR0.fs->gas);
    Gas_data &gR1 = *(cR1.fs->gas);
    Thermo_interpolator* ti = get_thermo_interpolator_ptr();
    int nsp = gmodel->get_number_of_species();
    int apply_limiter_flag = get_apply_limiter_flag();
    int extrema_clipping_flag = get_extrema_clipping_flag();
    
    // Low-order reconstruction just copies data from adjacent FV_Cell.
    Rght.copy_values_from(*(cR0.fs));
    if ( get_Xorder_flag() > 1 ) {
	// High-order reconstruction for some properties.
	onesided_interp_scalar(cL0.fs->vel.x, cR0.fs->vel.x, cR1.fs->vel.x,
			       cL0Length, cR0Length, cR1Length,
			       Rght.vel.x, apply_limiter_flag, extrema_clipping_flag);
	onesided_interp_scalar(cL0.fs->vel.y, cR0.fs->vel.y, cR1.fs->vel.y,
			       cL0Length, cR0Length, cR1Length,
			       Rght.vel.y, apply_limiter_flag, extrema_clipping_flag);
	onesided_interp_scalar(cL0.fs->vel.z, cR0.fs->vel.z, cR1.fs->vel.z,
			       cL0Length, cR0Length, cR1Length,
			       Rght.vel.z, apply_limiter_flag, extrema_clipping_flag);
	if ( get_mhd_flag() == 1 ) {
	    onesided_interp_scalar(cL0.fs->B.x, cR0.fs->B.x, cR1.fs->B.x,
				   cL0Length, cR0Length, cR1Length,
				   Rght.B.x, apply_limiter_flag, extrema_clipping_flag);
	    onesided_interp_scalar(cL0.fs->B.y, cR0.fs->B.y, cR1.fs->B.y,
				   cL0Length, cR0Length, cR1Length,
				   Rght.B.y, apply_limiter_flag, extrema_clipping_flag);
	    onesided_interp_scalar(cL0.fs->B.z, cR0.fs->B.z, cR1.fs->B.z,
				   cL0Length, cR0Length, cR1Length,
				   Rght.B.z, apply_limiter_flag, extrema_clipping_flag);
	}
	if ( get_k_omega_flag() == 1 ) {
	    onesided_interp_scalar(cL0.fs->tke, cR0.fs->tke, cR1.fs->tke,
				   cL0Length, cR0Length, cR1Length,
				   Rght.tke, apply_limiter_flag, extrema_clipping_flag);
	    onesided_interp_scalar(cL0.fs->omega, cR0.fs->omega, cR1.fs->omega,
				   cL0Length, cR0Length, cR1Length,
				   Rght.omega, apply_limiter_flag, extrema_clipping_flag);
	}
        for ( int isp = 0; isp < nsp; ++isp ) {
	    onesided_interp_scalar(cL0.fs->gas->massf[isp],
				   cR0.fs->gas->massf[isp], cR1.fs->gas->massf[isp],
				   cL0Length, cR0Length, cR1Length,
				   Rght.gas->massf[isp],
				   apply_limiter_flag, extrema_clipping_flag);
        }
	
	// Make the thermodynamic properties consistent.
	// Pressure, Local Speed of Sound and Temperature.
        // The value of 1 indicates that the old temperature
        // should be used as an initial guess for the iterative
        // EOS functions.
	// If the EOS call fouls up, just copy the cell data, low-order.
	if ( nsp > 1 ) {
	    if ( scale_mass_fractions( Rght.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp<Rght.gas->massf.size(); ++isp )
		    cR0.fs->gas->massf[isp] = Rght.gas->massf[isp];
	    }
	}
	
	// Interpolate on two of the thermodynamic quantities, and fill
	// in the rest based on an EOS call.
	
	onesided_interp_scalar(gL0.rho, gR0.rho, gR1.rho,
	                       cL0Length, cR0Length, cR1Length,
	                       Rght.gas->rho, apply_limiter_flag, extrema_clipping_flag);
	for ( int i = 0; i < gmodel->get_number_of_modes(); ++i ) {
	    onesided_interp_scalar(gL0.e[i], gR0.e[i], gR1.e[i],
				   cL0Length, cR0Length, cR1Length,
				   Rght.gas->e[i], apply_limiter_flag, extrema_clipping_flag);
	}
	int status = gmodel->eval_thermo_state_rhoe(*(Rght.gas));
	
	if ( status != SUCCESS ) {
	    // Rght state failed.
	    Rght.copy_values_from(*(cL0.fs));
	}			      
	
	if ( get_viscous_flag() ) {
	    gmodel->eval_transport_coefficients(*(Rght.gas));
	}
	if ( get_diffusion_flag() ) {
	    gmodel->eval_diffusion_coefficients(*(Rght.gas));
	}
    } // end of high-order reconstruction
    return SUCCESS;
} // end of onesided_interp()

/// \brief Reconstruct flow properties at an interface from FV_Cell properties.
///
/// This is essentially a one-dimensional linear interpolation process.  It needs only
/// the cell-average data and the lengths of the cells in the interpolation direction.

int one_d_linear_interp(FV_Cell &cL0, FV_Cell &cR0,
			double cL0Length, double cR0Length,
			FlowState &Lft)
{
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data &gL0 = *(cL0.fs->gas);
    Gas_data &gR0 = *(cR0.fs->gas);
    Thermo_interpolator* ti = get_thermo_interpolator_ptr();
    int nsp = gmodel->get_number_of_species();
    int apply_limiter_flag = get_apply_limiter_flag();
    int extrema_clipping_flag = get_extrema_clipping_flag();
    
    // Low-order reconstruction just copies data from adjacent FV_Cell.
    Lft.copy_values_from(*(cL0.fs));
    if ( get_Xorder_flag() > 1 ) {
	// High-order reconstruction for some properties.
	linear_interp_scalar(cL0.fs->vel.x, cR0.fs->vel.x,
			     cL0Length, cR0Length,
			     Lft.vel.x, apply_limiter_flag, extrema_clipping_flag);
	linear_interp_scalar(cL0.fs->vel.y, cR0.fs->vel.y,
			     cL0Length, cR0Length,
			     Lft.vel.y, apply_limiter_flag, extrema_clipping_flag);
	linear_interp_scalar(cL0.fs->vel.z, cR0.fs->vel.z,
			     cL0Length, cR0Length,
			     Lft.vel.z, apply_limiter_flag, extrema_clipping_flag);
	if ( get_mhd_flag() == 1 ) {
	    linear_interp_scalar(cL0.fs->B.x, cR0.fs->B.x,
				 cL0Length, cR0Length,
				 Lft.B.x, apply_limiter_flag, extrema_clipping_flag);
	    linear_interp_scalar(cL0.fs->B.y, cR0.fs->B.y,
				 cL0Length, cR0Length,
				 Lft.B.y, apply_limiter_flag, extrema_clipping_flag);
	    linear_interp_scalar(cL0.fs->B.z, cR0.fs->B.z,
				 cL0Length, cR0Length,
				 Lft.B.z, apply_limiter_flag, extrema_clipping_flag);
	}
	if ( get_k_omega_flag() == 1 ) {
	    linear_interp_scalar(cL0.fs->tke, cR0.fs->tke,
				 cL0Length, cR0Length,
				 Lft.tke, apply_limiter_flag, extrema_clipping_flag);
	    linear_interp_scalar(cL0.fs->omega, cR0.fs->omega,
				 cL0Length, cR0Length,
				 Lft.omega, apply_limiter_flag, extrema_clipping_flag);
	}
	for ( int isp = 0; isp < nsp; ++isp ) {
	    linear_interp_scalar(cL0.fs->gas->massf[isp], cR0.fs->gas->massf[isp],
				 cL0Length, cR0Length,
				 Lft.gas->massf[isp], apply_limiter_flag, extrema_clipping_flag);
	}
	
	// Make the thermodynamic properties consistent.
	// Pressure, Local Speed of Sound and Temperature.
        // The value of 1 indicates that the old temperature
        // should be used as an initial guess for the iterative
        // EOS functions.
	// If the EOS call fouls up, just copy the cell data, low-order.
	if ( nsp > 1 ) {
	    if ( scale_mass_fractions( Lft.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp<Lft.gas->massf.size(); ++isp )
		    cL0.fs->gas->massf[isp] = Lft.gas->massf[isp];
	    }
	}
	
	// Interpolate on two of the thermodynamic quantities, and fill
	// in the rest based on an EOS call.
	
	//if ( *ti = Rhoe_interpolator ) {
	linear_interp_scalar(gL0.rho, gR0.rho,
	                     cL0Length, cR0Length,
	                     Lft.gas->rho, apply_limiter_flag, extrema_clipping_flag);
	for ( int i = 0; i < gmodel->get_number_of_modes(); ++i ) {
	    linear_interp_scalar(gL0.e[i], gR0.e[i],
		                 cL0Length, cR0Length,
		                 Lft.gas->e[i], apply_limiter_flag, extrema_clipping_flag);
	}
	int status = gmodel->eval_thermo_state_rhoe(*(Lft.gas));		      
	
	if ( status != SUCCESS ) {
	    // Lft state failed.
	    Lft.copy_values_from(*(cL0.fs));
	}
	
	if ( get_viscous_flag() ) {
	    gmodel->eval_transport_coefficients(*(Lft.gas));
	}
	if ( get_diffusion_flag() ) {
	    gmodel->eval_diffusion_coefficients(*(Lft.gas));
	}
    } // end of high-order reconstruction
    return SUCCESS;
} // end of onesided_interp()


