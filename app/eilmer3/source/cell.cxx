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
#include <numeric>
#include <stdexcept>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/geometry2/source/geom.hh"
#include "cell.hh"
#include "kernel.hh"
#include "flux_calc.hh"
#include "diffusion.hh"
#include "bgk.hh"


enum update_scheme_t gasdynamic_update_scheme = PC_UPDATE;

update_scheme_t set_gasdynamic_update_scheme(update_scheme_t my_scheme)
{
    gasdynamic_update_scheme = my_scheme;
    return gasdynamic_update_scheme;
}

update_scheme_t get_gasdynamic_update_scheme()
{
    return gasdynamic_update_scheme;
}

std::string get_name_of_gasdynamic_update_scheme(update_scheme_t my_scheme)
{
    switch ( my_scheme ) {
    case EULER_UPDATE: return "euler";
    case PC_UPDATE: return "predictor-corrector";
    case MIDPOINT_UPDATE: return "midpoint";
    case CLASSIC_RK3_UPDATE: return "classic-rk3";
    case TVD_RK3_UPDATE: return "tvd-rk3";
    case DENMAN_RK3_UPDATE: return "denman-rk3";
    default: return "unknown";
    }
} // end get_name_of_gasdynamic_update_scheme()

size_t number_of_stages_for_update_scheme(update_scheme_t my_scheme)
{
    switch ( my_scheme ) {
    case EULER_UPDATE: return 1;
    case PC_UPDATE: return 2;
    case MIDPOINT_UPDATE: return 2;
    case CLASSIC_RK3_UPDATE: return 3;
    case TVD_RK3_UPDATE: return 3;
    case DENMAN_RK3_UPDATE: return 3;
    default: return 0;
    }
} // end number_of_stages_for_update_scheme()

/*----------------------------------------------------------------*/

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
    return -1; // "not a valid face" index
}

/// \brief Translates an index to a face name.
std::string get_face_name(int index)
{
    return ( index >= 0 && index <= 5 ) ? face_name[index] : "none";
}

//---------------------------------------------------------------------------

FlowState::FlowState(Gas_model *gm)
    : gas(new Gas_data(gm)),
      vel(0.0,0.0,0.0), B(0.0,0.0,0.0), S(0),
      tke(0.0), omega(0.0), mu_t(0.0), k_t(0.0), 
      G(get_velocity_buckets(),0.0), H(get_velocity_buckets(),0.0)
{}

FlowState::FlowState()
    : gas(new Gas_data(get_gas_model_ptr())),
      vel(0.0,0.0,0.0), B(0.0,0.0,0.0), S(0), 
      tke(0.0), omega(0.0), mu_t(0.0), k_t(0.0), 
      G(get_velocity_buckets(),0.0), H(get_velocity_buckets(),0.0)
{}

FlowState::FlowState(const FlowState &fs)
    : gas(new Gas_data(*(fs.gas))),
      vel(fs.vel), B(fs.B), S(fs.S), 
      tke(fs.tke), omega(fs.omega), mu_t(fs.mu_t), k_t(fs.k_t), 
      G(fs.G), H(fs.H)
{}

FlowState & FlowState::operator=(const FlowState &fs)
{
    if ( this != &fs ) { // Avoid aliasing
	gas = new Gas_data(*(fs.gas));
        vel = fs.vel; B = fs.B; S = fs.S; 
	tke = fs.tke; omega=fs.omega; mu_t = fs.mu_t; k_t = fs.k_t;
	G = fs.G; H = fs.H;
    }
    return *this;
}

FlowState::~FlowState()
{
    delete gas;
}

int FlowState::print() const
{
    global_data &gd = *get_global_data_ptr();
    printf("----------- Data for a flow state... ------------\n");
    gas->print_values();
    printf("v.x= %e, v.y= %e, v.z= %e \n", vel.x, vel.y, vel.z);
    if ( gd.MHD ) {
	printf("B.x= %e, B.y=%e, B.z=%e \n", B.x, B.y, B.z);
    }
    if ( get_velocity_buckets() > 0) {
	printf("Velocity Distribution Partial Densities:\n");
	for ( size_t ipd = 0; ipd < G.size(); ++ipd ) {
	    printf("i=%d: G=%e, H=%e\n", static_cast<int>(ipd), G[ipd], H[ipd]);
	}
    }
    printf("tke= %e, omega=%e\n", tke, omega);
    printf("S= %d, mu_t=%e, k_t=%e\n", S, mu_t, k_t);
    return SUCCESS;
}

int FlowState::copy_values_from(const FlowState &src)
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

int FlowState::copy_values_from(const CFlowCondition &src)
{
    global_data &G = *get_global_data_ptr();
    gas->copy_values_from(*(src.gas));
    vel.x = src.u; vel.y = src.v; vel.z = src.w;
    if ( G.MHD ) {
	B.x = src.Bx; B.y = src.By; B.z = src.Bz;
    }
    S = src.S;
    tke = src.tke;
    omega = src.omega;
    mu_t = src.mu_t;
    k_t = src.k_t;
    return SUCCESS;
}

int FlowState::average_values_from(const FlowState &src0, const FlowState &src1,
				   bool with_diff_coeff)
{
    global_data &gd = *get_global_data_ptr();
    gas->average_values_from(*(src0.gas), 0.5, *(src1.gas), 0.5, with_diff_coeff);
    vel.x = 0.5 * (src0.vel.x + src1.vel.x);
    vel.y = 0.5 * (src0.vel.y + src1.vel.y);
    vel.z = 0.5 * (src0.vel.z + src1.vel.z);
    if ( gd.MHD ) {
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
double * FlowState::copy_values_to_buffer(double *buf) const
{
    global_data &gd = *get_global_data_ptr();
    buf = gas->copy_values_to_buffer(buf);
    *buf++ = vel.x;
    *buf++ = vel.y;
    *buf++ = vel.z;
    if ( gd.MHD ) {
	*buf++ = B.x;
	*buf++ = B.y;
	*buf++ = B.z;
    }
    *buf++ = static_cast<double>(S);
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
    global_data &gd = *get_global_data_ptr();
    buf = gas->copy_values_from_buffer(buf);
    vel.x = *buf++;
    vel.y = *buf++;
    vel.z = *buf++;
    if ( gd.MHD ) {
	B.x = *buf++;
	B.y = *buf++;
	B.z = *buf++;
    }
    S = static_cast<int>(*buf++);
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
    : mass(0.0), momentum(0.0,0.0,0.0), B(0.0,0.0,0.0),
      total_energy(0.0),
      massf(gm->get_number_of_species(),0.0),
      energies(gm->get_number_of_modes(),0.0),
      tke(0.0), omega(0.0),
      G(get_velocity_buckets(),0.0), H(get_velocity_buckets(),0.0)
{}

ConservedQuantities::ConservedQuantities()
    : mass(0.0), momentum(0.0,0.0,0.0), B(0.0,0.0,0.0),
      total_energy(0.0),
      massf(get_gas_model_ptr()->get_number_of_species(),0.0),
      energies(get_gas_model_ptr()->get_number_of_modes(),0.0),
      tke(0.0), omega(0.0),
      G(get_velocity_buckets(),0.0), H(get_velocity_buckets(),0.0)
{}

ConservedQuantities::ConservedQuantities(const ConservedQuantities &Q)
    : mass(Q.mass), momentum(Q.momentum), B(Q.B),
      total_energy(Q.total_energy),
      massf(Q.massf),
      energies(Q.energies),
      tke(Q.tke), omega(Q.omega),
      G(Q.G), H(Q.H)
{}

ConservedQuantities & ConservedQuantities::operator=(const ConservedQuantities &Q)
{
    if ( this != &Q ) { // Avoid aliasing
	mass = Q.mass; momentum = Q.momentum; B = Q.B;
	total_energy = Q.total_energy;
	massf = Q.massf;
	energies = Q.energies;
	tke = Q.tke;
	G = Q.G; H = Q.H;
    }
    return *this;
}

ConservedQuantities::~ConservedQuantities()
{}

int ConservedQuantities::print() const
{
    global_data &gd = *get_global_data_ptr();
    cout << "mass= " << mass << endl;
    cout << "momentum.x= " << momentum.x << " .y= " << momentum.y
	 << " .z= " << momentum.z << endl;
    if ( gd.MHD ) {
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

int ConservedQuantities::copy_values_from(const ConservedQuantities &src)
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
    : id(0), status(0), pos(0.0,0.0,0.0), vel(0.0,0.0,0.0),
      Ybar(0.0), length(0.0), area(N_LEVEL, 0.0),
      n(0.0,0.0,0.0), t1(0.0,0.0,0.0), t2(0.0,0.0,0.0),
      fs(new FlowState(gm)), F(new ConservedQuantities(gm))
{
#   if WITH_IMPLICIT == 1
    // Ojas, FIX-ME, may be good to initialize point-implicit variables
    // double Lambda[6][2], R[6][6], R_inv[6][6], sl[6][2], sl_min[6][2];
    // double A[6][6], J[6][6], hl[6][2], gl[6][2];
#   endif
}

FV_Interface::FV_Interface()
    : id(0), status(0), pos(0.0,0.0,0.0), vel(0.0,0.0,0.0),
      Ybar(0.0), length(0.0), area(N_LEVEL, 0.0),
      n(0.0,0.0,0.0), t1(0.0,0.0,0.0), t2(0.0,0.0,0.0),
      fs(new FlowState(get_gas_model_ptr())),
      F(new ConservedQuantities(get_gas_model_ptr()))
{}

FV_Interface::FV_Interface(const FV_Interface &iface)
    : id(iface.id), status(iface.status), pos(iface.pos), vel(iface.vel),
      Ybar(iface.Ybar), length(iface.length), area(iface.area), 
      n(iface.n), t1(iface.t1), t2(iface.t2),
      fs(new FlowState(*iface.fs)), F(new ConservedQuantities(*iface.F))
{}

FV_Interface & FV_Interface::operator=(const FV_Interface &iface)
{
    if ( this != &iface ) { // Avoid aliasing
	id = iface.id; status = iface.status;
	pos = iface.pos; vel = iface.vel;
	Ybar = iface.Ybar; length = iface.length; area = iface.area;
	n = iface.n; t1 = iface.t1; t2 = iface.t2;
	fs = new FlowState(*iface.fs);
	F = new ConservedQuantities(*iface.F);
    }
    return *this;
}

FV_Interface::~FV_Interface()
{
    delete fs;
    delete F;
}

int FV_Interface::print() const
{
    printf( "----------- Begin data for interface -----------\n");
    // printf( "id = %i\n", iface->id);
    printf("x=%e, y=%e, z=%e\n", pos.x, pos.y, pos.z);
    for ( size_t i = 0; i < area.size(); ++i ) {
	printf("area[%d]=%e, ", static_cast<int>(i), area[i]);
    }
    printf("\nYbar=%e, length=%e\n", Ybar, length);
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
int FV_Interface::copy_values_from(const FV_Interface &src, int type_of_copy)
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
	Ybar = src.Ybar; length = src.length; 
	for ( size_t i = 0; i < area.size(); ++i ) area[i] = src.area[i];
    }
    return SUCCESS;
}

int FV_Interface::copy_grid_level_to_level(size_t from_level, size_t to_level)
{
    area[to_level] = area[from_level];
    return SUCCESS;
}

//----------------------------------------------------------------------------

FV_Vertex::FV_Vertex(Gas_model *gm)
    : id(0), pos(N_LEVEL,Vector3(0.0,0.0,0.0)),
      vel(N_LEVEL,Vector3(0.0,0.0,0.0)),
      area(0.0), volume(0.0),
      dudx(0.0), dudy(0.0), dudz(0.0),
      dvdx(0.0), dvdy(0.0), dvdz(0.0),
      dwdx(0.0), dwdy(0.0), dwdz(0.0),
      dTdx(gm->get_number_of_modes(),0.0),
      dTdy(gm->get_number_of_modes(),0.0),
      dTdz(gm->get_number_of_modes(),0.0),
      dtkedx(0.0), dtkedy(0.0), dtkedz(0.0),
      domegadx(0.0), domegady(0.0), domegadz(0.0),
      dfdx(gm->get_number_of_species(),0.0),
      dfdy(gm->get_number_of_species(),0.0),
      dfdz(gm->get_number_of_species(),0.0),
      dpedx(0.0), dpedy(0.0), dpedz(0.0)
{}

FV_Vertex::FV_Vertex()
    : id(0), pos(N_LEVEL,Vector3(0.0,0.0,0.0)),
      vel(N_LEVEL,Vector3(0.0,0.0,0.0)),
      area(0.0), volume(0.0),
      dudx(0.0), dudy(0.0), dudz(0.0),
      dvdx(0.0), dvdy(0.0), dvdz(0.0),
      dwdx(0.0), dwdy(0.0), dwdz(0.0),
      dTdx(get_gas_model_ptr()->get_number_of_modes(),0.0),
      dTdy(get_gas_model_ptr()->get_number_of_modes(),0.0),
      dTdz(get_gas_model_ptr()->get_number_of_modes(),0.0),
      dtkedx(0.0), dtkedy(0.0), dtkedz(0.0),
      domegadx(0.0), domegady(0.0), domegadz(0.0),
      dfdx(get_gas_model_ptr()->get_number_of_species(),0.0),
      dfdy(get_gas_model_ptr()->get_number_of_species(),0.0),
      dfdz(get_gas_model_ptr()->get_number_of_species(),0.0),
      dpedx(0.0), dpedy(0.0), dpedz(0.0)
{}

FV_Vertex::FV_Vertex(const FV_Vertex &vtx)
    : id(vtx.id), pos(vtx.pos), vel(vtx.vel),
      area(vtx.area), volume(vtx.volume),
      dudx(vtx.dudx), dudy(vtx.dudy), dudz(vtx.dudz),
      dvdx(vtx.dvdx), dvdy(vtx.dvdy), dvdz(vtx.dvdz),
      dwdx(vtx.dwdx), dwdy(vtx.dwdy), dwdz(vtx.dwdz),
      dTdx(vtx.dTdx), dTdy(vtx.dTdy), dTdz(vtx.dTdz),
      dtkedx(vtx.dtkedx), dtkedy(vtx.dtkedy), dtkedz(vtx.dtkedz),
      domegadx(vtx.domegadx), domegady(vtx.domegady), domegadz(vtx.domegadz),
      dfdx(vtx.dfdx), dfdy(vtx.dfdy), dfdz(vtx.dfdz),
      dpedx(vtx.dpedx), dpedy(vtx.dpedy), dpedz(vtx.dpedz)
{}

FV_Vertex::~FV_Vertex()
{}

int FV_Vertex::copy_values_from(const FV_Vertex &src)
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
    dpedx = src.dpedx; dpedy = src.dpedy; dpedz = src.dpedz;
    return SUCCESS;
}

int FV_Vertex::copy_grid_level_to_level(size_t from_level, size_t to_level)
{
    pos[to_level].x = pos[from_level].x;
    pos[to_level].y = pos[from_level].y;
    pos[to_level].z = pos[from_level].z;
    vel[to_level].x = vel[from_level].x;
    vel[to_level].y = vel[from_level].y;
    vel[to_level].z = vel[from_level].z;
    return SUCCESS;
}


//----------------------------------------------------------------------------

FV_Cell::FV_Cell(Gas_model *gm)
    : id(0), status(NORMAL_CELL), fr_reactions_allowed(false),
      dt_chem(0.0), dt_therm(0.0), in_turbulent_zone(false),
      base_qdot(0.0), pos(N_LEVEL,Vector3(0.0,0.0,0.0)),
      volume(N_LEVEL,0.0), area(N_LEVEL,0.0), uf(0.0),
      iLength(0.0), jLength(0.0), kLength(0.0), L_min(0.0),
      distance_to_nearest_wall(0.0), half_cell_width_at_wall(0.0),
      cell_at_nearest_wall(NULL), iface(N_INTERFACE,NULL), vtx(N_VERTEX,NULL),
      fs(new FlowState(gm)), Q(new ConservedQuantities(gm)),
      Q_rad_org(0.0), f_rad_org(0.0), Q_rE_rad(0.0),
      rho_at_start_of_step(0.0), rE_at_start_of_step(0.0)
{
    // Maybe we could put the following into the init section as
    // dUdt(N_LEVEL,new ConservedQuantities(gm))
    for ( size_t i = 0; i < N_LEVEL; ++i ) {
	U.push_back(new ConservedQuantities(gm));
	dUdt.push_back(new ConservedQuantities(gm));
    }
#   if WITH_IMPLICIT == 1
    // Ojas, it may be good to initialize point-implicit variables.
    // double piM[6][6], pir[6][2], qL[6][2], g_Ll[6][2], gjvec[6][2], gjmtx[6][6];
#   endif
}

FV_Cell::FV_Cell()
    : id(0), status(NORMAL_CELL), fr_reactions_allowed(false),
      dt_chem(0.0), dt_therm(0.0), in_turbulent_zone(false),
      base_qdot(0.0), pos(N_LEVEL,Vector3(0.0,0.0,0.0)),
      volume(N_LEVEL,0.0), area(N_LEVEL,0.0), uf(0.0),
      iLength(0.0), jLength(0.0), kLength(0.0), L_min(0.0),
      distance_to_nearest_wall(0.0), half_cell_width_at_wall(0.0),
      cell_at_nearest_wall(NULL), iface(N_INTERFACE,NULL), vtx(N_VERTEX,NULL),
      fs(new FlowState(get_gas_model_ptr())),
      Q(new ConservedQuantities(get_gas_model_ptr())),
      Q_rad_org(0.0), f_rad_org(0.0), Q_rE_rad(0.0),
      rho_at_start_of_step(0.0), rE_at_start_of_step(0.0)
{
    Gas_model *gm = get_gas_model_ptr();
    for ( size_t i = 0; i < N_LEVEL; ++i ) {
	U.push_back(new ConservedQuantities(gm));
	dUdt.push_back(new ConservedQuantities(gm));
    }
}

FV_Cell::FV_Cell(const FV_Cell &cell)
    : id(cell.id), status(cell.status), 
      fr_reactions_allowed(cell.fr_reactions_allowed),
      dt_chem(cell.dt_chem), dt_therm(cell.dt_therm),
      in_turbulent_zone(cell.in_turbulent_zone),
      base_qdot(cell.base_qdot), pos(cell.pos),
      volume(cell.volume), area(cell.area), uf(cell.uf),
      iLength(cell.iLength), jLength(cell.jLength), kLength(cell.kLength),
      L_min(cell.L_min),
      distance_to_nearest_wall(cell.distance_to_nearest_wall),
      half_cell_width_at_wall(cell.half_cell_width_at_wall),
      cell_at_nearest_wall(cell.cell_at_nearest_wall),
      iface(cell.iface), vtx(cell.vtx),
      fs(new FlowState(*cell.fs)), Q(cell.Q),
      Q_rad_org(cell.Q_rad_org), f_rad_org(cell.f_rad_org), Q_rE_rad(cell.Q_rE_rad),
      rho_at_start_of_step(cell.rho_at_start_of_step),
      rE_at_start_of_step(cell.rE_at_start_of_step)
{
    for ( size_t i = 0; i < N_LEVEL; ++i ) {
	U.push_back(new ConservedQuantities(*cell.U[i]));
	dUdt.push_back(new ConservedQuantities(*cell.dUdt[i]));
    }
#   if WITH_IMPLICIT == 1
    // Ojas, it may be good to initialize point-implicit variables.
    // double piM[6][6], pir[6][2], qL[6][2], g_Ll[6][2], gjvec[6][2], gjmtx[6][6];
#   endif
}

FV_Cell & FV_Cell::operator=(const FV_Cell &cell)
{
    if ( this != &cell ) { // Avoid aliasing
	id = cell.id; status = cell.status;
	fr_reactions_allowed = cell.fr_reactions_allowed;
	dt_chem = cell.dt_chem; dt_therm = cell.dt_therm;
	in_turbulent_zone = cell.in_turbulent_zone;
	base_qdot = cell.base_qdot;
	pos = cell.pos; volume = cell.volume;
	area = cell.area; uf = cell.uf;
	iLength = cell.iLength; jLength = cell.jLength;
	kLength = cell.kLength; L_min = cell.L_min;
	distance_to_nearest_wall = cell.distance_to_nearest_wall;
	half_cell_width_at_wall = cell.half_cell_width_at_wall;
	cell_at_nearest_wall = cell.cell_at_nearest_wall;
	iface =cell.iface; vtx = cell.vtx;
	delete fs; fs = new FlowState(*cell.fs);
	for ( size_t i = 0; i < dUdt.size(); ++i ) {
	    delete U[i];
	    delete dUdt[i];
	}
	for ( size_t i = 0; i < N_LEVEL; ++i ) {
	    U.push_back(new ConservedQuantities(*cell.U[i]));
	    dUdt.push_back(new ConservedQuantities(*cell.dUdt[i]));
	}
	Q = cell.Q;
	Q_rad_org = cell.Q_rad_org;
	f_rad_org = cell.f_rad_org;
	Q_rE_rad = cell.Q_rE_rad;
	rho_at_start_of_step = cell.rho_at_start_of_step;
	rE_at_start_of_step = cell.rE_at_start_of_step;
    }
    return *this;
}

FV_Cell::~FV_Cell()
{
    delete fs;
    for ( size_t i = 0; i < N_LEVEL; ++i ) {
	delete U[i];
	delete dUdt[i];
    }
    delete Q;
}

int FV_Cell::print() const
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    printf("----------- Begin data for cell -----------\n");
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();
    printf("nsp=%d, nmodes=%d\n", nsp, nmodes);
    // We'll just report grid-level 0 position.
    printf("x=%e, y=%e, z=%e, base_qdot=%e\n", pos[0].x, pos[0].y, pos[0].z, base_qdot);
    fs->print();
    if ( G.radiation ) {
	printf("radiation source Q_rE_rad=%e\n", Q_rE_rad);
    }
    printf("Conserved quantities: \n");
    for ( size_t i = 0; i < U.size(); ++i ) {
	printf("stage %d:\n", static_cast<int>(i));
	U[i]->print();
    }
    printf("status= %d \n", status);
    printf( "----------- End data for cell -----------\n");
    return SUCCESS;
}

int FV_Cell::point_is_inside(Vector3 &p, int dimensions, size_t gtl) const
/// \brief Returns 1 if the point p is inside or on the cell surface.
{
    if ( dimensions == 2 ) {
	// In 2 dimensions,
	// we split the x,y-plane into half-planes and check which side p is on.
	double xA = vtx[1]->pos[gtl].x; double yA = vtx[1]->pos[gtl].y;
	double xB = vtx[1]->pos[gtl].x; double yB = vtx[2]->pos[gtl].y;
	double xC = vtx[3]->pos[gtl].x; double yC = vtx[3]->pos[gtl].y;
	double xD = vtx[0]->pos[gtl].x; double yD = vtx[0]->pos[gtl].y;
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
	if ( tetrahedron_volume(vtx[2]->pos[gtl], vtx[3]->pos[gtl], vtx[7]->pos[gtl], p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[7]->pos[gtl], vtx[6]->pos[gtl], vtx[2]->pos[gtl], p) > 0.0 ) return 0;
	// East
	if ( tetrahedron_volume(vtx[1]->pos[gtl], vtx[2]->pos[gtl], vtx[6]->pos[gtl], p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[6]->pos[gtl], vtx[5]->pos[gtl], vtx[1]->pos[gtl], p) > 0.0 ) return 0;
	// South
	if ( tetrahedron_volume(vtx[0]->pos[gtl], vtx[1]->pos[gtl], vtx[5]->pos[gtl], p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[5]->pos[gtl], vtx[4]->pos[gtl], vtx[0]->pos[gtl], p) > 0.0 ) return 0;
	// West
	if ( tetrahedron_volume(vtx[3]->pos[gtl], vtx[0]->pos[gtl], vtx[4]->pos[gtl], p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[4]->pos[gtl], vtx[7]->pos[gtl], vtx[3]->pos[gtl], p) > 0.0 ) return 0;
	// Bottom
	if ( tetrahedron_volume(vtx[1]->pos[gtl], vtx[0]->pos[gtl], vtx[3]->pos[gtl], p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[3]->pos[gtl], vtx[2]->pos[gtl], vtx[1]->pos[gtl], p) > 0.0 ) return 0;
	// Top
	if ( tetrahedron_volume(vtx[4]->pos[gtl], vtx[5]->pos[gtl], vtx[6]->pos[gtl], p) > 0.0 ) return 0;
	if ( tetrahedron_volume(vtx[6]->pos[gtl], vtx[7]->pos[gtl], vtx[4]->pos[gtl], p) > 0.0 ) return 0;
	// If we arrive here, we haven't determined that the point is outside...
	return 1;
    } // end dimensions != 2
} // end point_is_inside()

int FV_Cell::copy_values_from(const CFlowCondition &src)
{
    global_data &G = *get_global_data_ptr();
    fs->gas->copy_values_from(*(src.gas));
    fs->vel.x = src.u; fs->vel.y = src.v; fs->vel.z = src.w;
    if ( G.MHD ) {
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
int FV_Cell::copy_values_from(const FV_Cell &src, int type_of_copy, size_t gtl)
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
	pos[0].x = src.pos[gtl].x;
	pos[0].y = src.pos[gtl].y;
	pos[0].z = src.pos[gtl].z;
    }
    if (type_of_copy == COPY_INTERFACE_DATA) {
	for ( size_t j = 0; j < N_INTERFACE; ++j ) {
	    if ( src.iface[j] == 0 || iface[j] == 0 ) {
		// When copying from ghost cell which may
		// not have initialised interfaces, or when 2D.
		continue;
	    }
	    // [todo] FIX-ME -- I'm not happy with conditional copies here -- PJ 23-Mar-2013.
	    iface[j]->copy_values_from(*(src.iface[j]), COPY_ALL_CELL_DATA);
	}
    }
    return SUCCESS;
}


/// \brief Copy the cell data into a linear data buffer.
///
/// \param buf : pointer to the current element somewhere in buffer
/// \returns a pointer to the next available location in the data buffer.
double * FV_Cell::copy_values_to_buffer(double *buf, int type_of_copy, size_t gtl) const
{
    if (type_of_copy == COPY_ALL_CELL_DATA ||
	type_of_copy == COPY_FLOW_STATE) {
        buf = fs->copy_values_to_buffer(buf);
	*buf++ = static_cast<double>(status);
	*buf++ = Q_rE_rad;
    }
    if (type_of_copy == COPY_ALL_CELL_DATA ||
        type_of_copy == COPY_CELL_LENGTHS) {
        *buf++ = iLength; *buf++ = jLength; *buf++ = kLength;
        *buf++ = pos[gtl].x; *buf++ = pos[gtl].y; *buf++ = pos[gtl].z;
    }
    if (type_of_copy == COPY_INTERFACE_DATA) {
	for ( size_t j = 0; j < N_INTERFACE; ++j ) {
	    if ( iface[j] == 0 ) { // When copying from ghost cell which may
		continue;          // not have initialised interfaces, or when 2D.
	    }
	    // [todo] FIX-ME -- I'm not happy with conditional copies here -- PJ 23-Mar-2013.
	    buf = iface[j]->fs->copy_values_to_buffer(buf);
	    *buf++ = iface[j]->pos.x; *buf++ = iface[j]->pos.y; *buf++ = iface[j]->pos.z;
	    *buf++ = iface[j]->vel.x; *buf++ = iface[j]->vel.y; *buf++ = iface[j]->vel.z;
	    *buf++ = iface[j]->length;
	}
    }
    return buf;
}

/// \brief Copy the data from a linear data buffer into the cell.
/// \param buf : pointer to the current element somewhere in buffer
/// \returns a pointer to the next available location in the data buffer.
double * FV_Cell::copy_values_from_buffer(double *buf, int type_of_copy, size_t gtl)
{
    if (type_of_copy == COPY_ALL_CELL_DATA ||
	type_of_copy == COPY_FLOW_STATE) {
	buf = fs->copy_values_from_buffer(buf);
	status = static_cast<int>(*buf++);
	Q_rE_rad = *buf++;
    }
    if (type_of_copy == COPY_ALL_CELL_DATA ||
        type_of_copy == COPY_CELL_LENGTHS) {
        iLength = *buf++; jLength = *buf++; kLength = *buf++;
        pos[gtl].x = *buf++; pos[gtl].y = *buf++; pos[gtl].z = *buf++;
    }
    if (type_of_copy == COPY_INTERFACE_DATA) {
	for ( size_t j = 0; j < N_INTERFACE; ++j ) {
	    if ( iface[j] == 0 ) { // When copying from ghost cell which may
		continue;          // not have initialised interfaces, or when 2D.
	    }
	    // [todo] FIX-ME -- I'm not happy with conditional copies here -- PJ 23-Mar-2013.
	    buf = iface[j]->fs->copy_values_from_buffer(buf);
	    iface[j]->pos.x = *buf++; iface[j]->pos.y = *buf++; iface[j]->pos.z = *buf++;
	    iface[j]->vel.x = *buf++; iface[j]->vel.y = *buf++; iface[j]->vel.z = *buf++;
	    iface[j]->length = *buf++;
	}
    }
    return buf;
}


int FV_Cell::copy_grid_level_to_level(size_t from_level, size_t to_level)
{
    pos[to_level].x = pos[from_level].x;
    pos[to_level].y = pos[from_level].y;
    pos[to_level].z = pos[from_level].z;
    // When working over all cells in a block, the following copies
    // will no doubt do some doubled-up work, but it should be otherwise benign.
    for ( FV_Interface *face : iface ) {
	if ( face ) face->copy_grid_level_to_level(from_level, to_level);
    }
    for ( FV_Vertex *v : vtx ) {
	if ( v ) v->copy_grid_level_to_level(from_level, to_level);
    }
    return SUCCESS;
}


/// \brief Replace the flow data in a cell with the average from neighbour cells.
int FV_Cell::replace_flow_data_with_average(std::vector<FV_Cell *> src)
{
    global_data &G = *get_global_data_ptr();
    size_t ncell = src.size(); 
    size_t ii = 0;
    if ( ncell < 1 ) {
	throw std::runtime_error("FV_Cell::replace_flow_data_with_average(): "
				 "The list of cells was empty.");
    }

    /* First, replace the crappy data with that from the first cell. */
    ii = 0;
    fs->gas->accumulate_values_from(*(src[ii]->fs->gas), 0.0);
    fs->vel.x = src[ii]->fs->vel.x;
    fs->vel.y = src[ii]->fs->vel.y;
    fs->vel.z = src[ii]->fs->vel.z;
    if ( G.MHD ) {
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
	    if ( G.MHD ) {
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
	double dncell = static_cast<double>(ncell);
	fs->gas->accumulate_values_from(*(fs->gas), (1.0-ncell)/dncell );
	fs->vel.x /= dncell;
	fs->vel.y /= dncell;
	fs->vel.z /= dncell;
	if ( G.MHD ) {
	    fs->B.x /= dncell;
	    fs->B.y /= dncell;
	    fs->B.z /= dncell;
	}
	fs->mu_t /= dncell;
	fs->k_t /= dncell;
	fs->tke /= dncell;
	fs->omega /= dncell;
	Q_rE_rad /= dncell;
    }
    // The following calls are expensive but getting to this point should be very rare.
    // If it is common, we have debugging to do...
    Gas_model *gmodel = get_gas_model_ptr();
    gmodel->eval_thermo_state_pT(*(fs->gas));
    if ( G.viscous ) gmodel->eval_transport_coefficients(*(fs->gas));
    if ( G.diffusion ) gmodel->eval_diffusion_coefficients(*(fs->gas));

    return SUCCESS;
} // end function replace_cell_data_with_average()


/// \brief Scan a string, extracting the flow solution (i.e. the primary variables).
int FV_Cell::scan_values_from_string(char *bufptr)
// There isn't any checking of the file content.
// If anything gets out of place, the result is wrong data.
{
    global_data &G = *get_global_data_ptr();
    size_t nmodes = fs->gas->T.size();
    size_t nsp = fs->gas->massf.size();
    // Look for a new-line character and truncate the string there.
    char *cptr = strchr(bufptr, '\n');
    if ( cptr != NULL ) *cptr = '\0';
    // Now, we should have a string with only numbers separated by spaces.
    pos[0].x = atof(strtok( bufptr, " " )); // Note grid-level 0.
    pos[0].y = atof(strtok( NULL, " " ));
    pos[0].z = atof(strtok( NULL, " " ));
    volume[0] = atof(strtok( NULL, " " ));
    fs->gas->rho = atof(strtok( NULL, " " ));
    fs->vel.x = atof(strtok( NULL, " " ));
    fs->vel.y = atof(strtok( NULL, " " ));
    fs->vel.z = atof(strtok( NULL, " " ));
    if ( G.MHD ) {
	fs->B.x = atof(strtok( NULL, " " ));
	fs->B.y = atof(strtok( NULL, " " ));
	fs->B.z = atof(strtok( NULL, " " ));
    }
    fs->gas->p = atof(strtok( NULL, " " ));
    fs->gas->a = atof(strtok( NULL, " " ));
    fs->gas->mu = atof(strtok( NULL, " " ));
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	fs->gas->k[imode] = atof(strtok( NULL, " " ));
    }
    fs->mu_t = atof(strtok( NULL, " " ));
    fs->k_t = atof(strtok( NULL, " " ));
    fs->S = atoi(strtok( NULL, " " ));
    if ( G.radiation ) {
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
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	fs->gas->massf[isp] = atof(strtok( NULL, " " ));
    }
    if ( nsp > 1 ) dt_chem = atof(strtok( NULL, " " ));
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	fs->gas->e[imode] = atof(strtok( NULL, " " ));
	fs->gas->T[imode] = atof(strtok( NULL, " " ));
    }
    if ( nmodes > 1 ) dt_therm = atof(strtok( NULL, " " ));
    return SUCCESS;
} // end scan_values_from_string()

/// \brief Write the flow solution (i.e. the primary variables) to a string.
std::string FV_Cell::write_values_to_string() const
{
    global_data &G = *get_global_data_ptr();
    size_t nsp = fs->gas->massf.size();
    size_t nmodes = fs->gas->T.size();
    // The new format for Elmer3 puts everything onto one line.
    ostringstream ost;
    ost.setf(ios_base::scientific);
    ost.precision(12);
    ost << pos[0].x << " " << pos[0].y << " " << pos[0].z // Note grid-level 0.
	<< " " << volume[0] << " " <<  fs->gas->rho
	<< " " << fs->vel.x << " " << fs->vel.y << " " << fs->vel.z;
    if ( G.MHD ) {
	ost << " " << fs->B.x << " " << fs->B.y << " " << fs->B.z;
    }
    ost << " " << fs->gas->p << " " << fs->gas->a << " " << fs->gas->mu;
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	ost << " " << fs->gas->k[imode];
    }
    ost << " " << fs->mu_t << " " << fs->k_t << " " << fs->S;
    if ( G.radiation ) {
	ost << " " << Q_rad_org << " " << f_rad_org << " " << Q_rE_rad;
    }
    ost << " " << fs->tke << " " << fs->omega;
    // Species mass fractions.
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	ost << " " << fs->gas->massf[isp];
    }
    if ( nsp > 1 ) ost << " " << dt_chem;
    // Individual energies (in e, T pairs)
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
    if ( cptr != NULL ) *cptr = '\0';
    // Now, we should have a string with only numbers separated by spaces.
    // include the position data, duplicate of "flow", as insurance against 
    // finding this file in isolation
    // Note that grid-level zero is the destination of the data.
    pos[0].x = atof(strtok( bufptr, " " ));
    pos[0].y = atof(strtok( NULL, " " ));
    pos[0].z = atof(strtok( NULL, " " ));
    volume[0] = atof(strtok( NULL, " " ));
    // values of G and H are interleaved
    for (size_t iGH = 0; iGH < get_velocity_buckets(); ++iGH) {
	fs->G[iGH] = atof(strtok( NULL, " " ));
	fs->H[iGH] = atof(strtok( NULL, " " ));
    }
    return SUCCESS;
} // end scan_BGK_from_string()

/// \brief Write the discrete samples of the BGK velocity distribution function to a string.
std::string FV_Cell::write_BGK_to_string() const
{
    // The new format for Elmer3 puts everything onto one line.
    ostringstream ost;
    ost.setf(ios_base::scientific);
    ost.precision(12);
    // Note that grid-level is zero.
    ost << pos[0].x << " " << pos[0].y << " " << pos[0].z << " " << volume[0];
    // BGK discrete samples of velocity distribution function
    // interleave G and H
    for ( size_t iGH = 0; iGH < get_velocity_buckets(); ++iGH ) {
	ost << " " << fs->G[iGH];
	ost << " " << fs->H[iGH];
    }
    // Don't put the newline char on the end.
    return ost.str();
} // end of write_BGK_to_string()


int FV_Cell::encode_conserved(size_t gtl, size_t ftl, double omegaz, bool with_k_omega)
{
    global_data &G = *get_global_data_ptr();
    ConservedQuantities &myU = *(U[ftl]);

    myU.mass = fs->gas->rho;
    // X-, Y- and Z-momentum per unit volume.
    myU.momentum.x = fs->gas->rho * fs->vel.x;
    myU.momentum.y = fs->gas->rho * fs->vel.y;
    myU.momentum.z = fs->gas->rho * fs->vel.z;
    // Magnetic field
    myU.B.x = fs->B.x;
    myU.B.y = fs->B.y;
    myU.B.z = fs->B.z;
    // Total Energy / unit volume = density
    // (specific internal energy + kinetic energy/unit mass).
    double e = accumulate(fs->gas->e.begin(), fs->gas->e.end(), 0.0);
    double ke = 0.5 * (fs->vel.x * fs->vel.x
		       + fs->vel.y * fs->vel.y
		       + fs->vel.z * fs->vel.z);
    if ( with_k_omega ) {
	myU.tke = fs->gas->rho * fs->tke;
	myU.omega = fs->gas->rho * fs->omega;
	myU.total_energy = fs->gas->rho * (e + ke + fs->tke);
    } else {
	myU.tke = 0.0;
	myU.omega = fs->gas->rho * 1.0;
	myU.total_energy = fs->gas->rho * (e + ke);
    }
    if ( G.MHD ) {
	double me = 0.5 * (fs->B.x * fs->B.x
			   + fs->B.y * fs->B.y
			   + fs->B.z * fs->B.z);
	myU.total_energy += me;
    }
    // Species densities: mass of species is per unit volume.
    for ( size_t isp = 0; isp < myU.massf.size(); ++isp ) {
	myU.massf[isp] = fs->gas->rho * fs->gas->massf[isp];
    }
    // Individual energies: energy in mode per unit volume
    for ( size_t imode = 0; imode < myU.energies.size(); ++imode ) {
	myU.energies[imode] = fs->gas->rho * fs->gas->e[imode];
    }
    
    if ( omegaz != 0.0 ) {
	// Rotating frame.
	// Finally, we adjust the total energy to make rothalpy.
	// We do this last because the gas models don't know anything
	// about rotating frames and we don't want to mess their
	// energy calculations around.
	double rho = fs->gas->rho;
	double x = pos[gtl].x;
	double y = pos[gtl].y;
	double rsq = x*x + y*y;
	// The conserved quantity is rothalpy. I = E - (u**2)/2
	// where rotating frame velocity  u = omegaz * r.
	myU.total_energy -= rho * 0.5 * omegaz * omegaz * rsq;
    }
    return SUCCESS;
} // end of encode_conserved()


int FV_Cell::decode_conserved(size_t gtl, size_t ftl, double omegaz, bool with_k_omega)
{
    global_data &G = *get_global_data_ptr();
    ConservedQuantities &myU = *(U[ftl]);
    Gas_model *gmodel = get_gas_model_ptr();
    double e, ke, dinv, rE, me;

    // Mass / unit volume = Density
    double rho = myU.mass;
    fs->gas->rho = rho; // This is limited to nonnegative and finite values.
    if ( rho <= 0.0 ) {
	cout << "FV_Cell::decode_conserved(): Density is below minimum rho= " 
	     << rho << endl;
	cout << "id= " << id << " x= " << pos[gtl].x << " y= " << pos[gtl].y 
	     << " z= " << pos[gtl].z << endl;
	fs->gas->print_values();
    }
    dinv = 1.0 / rho;
    if ( omegaz != 0.0 ) {
	// Rotating frame.
	// The conserved quantity is rothalpy so we need to convert
	// back to enthalpy to do the rest of the decode.
	double x = pos[gtl].x;
	double y = pos[gtl].y;
	double rsq = x*x + y*y;
	rE = myU.total_energy + rho * 0.5 * omegaz * omegaz * rsq;
    } else {
	// Non-rotating frame.
	rE = myU.total_energy;
    }
    // Velocities from momenta.
    fs->vel.x = myU.momentum.x * dinv;
    fs->vel.y = myU.momentum.y * dinv;
    fs->vel.z = myU.momentum.z * dinv;
    // Magnetic field
    fs->B.x = myU.B.x;
    fs->B.y = myU.B.y;
    fs->B.z = myU.B.z;
    // Specific internal energy from total energy per unit volume.
    ke = 0.5 * (fs->vel.x * fs->vel.x + fs->vel.y * fs->vel.y + fs->vel.z * fs->vel.z);
    if ( G.MHD ) {
        me = 0.5*(fs->B.x*fs->B.x + fs->B.y*fs->B.y + fs->B.z*fs->B.z);
    } else {
        me = 0.0;
    }
    if ( with_k_omega ) {
        fs->tke = myU.tke * dinv;
        fs->omega = myU.omega * dinv;
	e = (rE - myU.tke - me) * dinv - ke;
    } else {
        fs->tke = 0.0;
        fs->omega = 1.0;
	e = (rE - me) * dinv - ke;
    }
    size_t nsp = myU.massf.size();
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	fs->gas->massf[isp] = myU.massf[isp] * dinv;
    }
    if ( nsp > 1 ) scale_mass_fractions( fs->gas->massf );
    size_t nmodes = myU.energies.size();
    for ( size_t imode = 1; imode < nmodes; ++imode ) {
	fs->gas->e[imode] = myU.energies[imode] * dinv;
    }
    // We can recompute e[0] from total energy and component
    // modes NOT in translation.
    if ( nmodes > 1 ) {
	double e_tmp = accumulate(fs->gas->e.begin()+1, fs->gas->e.end(), 0.0);
	fs->gas->e[0] = e - e_tmp;
    }
    else {
	fs->gas->e[0] = e;
    }
    // Fill out the other variables; P, T, a and
    // check the species mass fractions.
    // Update the viscous transport coefficients.
    gmodel->eval_thermo_state_rhoe(*(fs->gas));
    if ( G.viscous ) gmodel->eval_transport_coefficients(*(fs->gas));
    if ( G.diffusion ) gmodel->eval_diffusion_coefficients(*(fs->gas));

    return SUCCESS;
} // end of decode_conserved()


/// \brief Check the primary flow data for a specified cell.
/// \returns true for valid data, false for bad data.
bool FV_Cell::check_flow_data(void)
{
    bool is_data_valid = true;
    is_data_valid = fs->gas->check_values(true);
    constexpr double MAXVEL = 30000.0;
    if (fabs(fs->vel.x) > MAXVEL || fabs(fs->vel.y) > MAXVEL || fabs(fs->vel.z) > MAXVEL) {
	cout << "Velocity bad " << fs->vel.x << " " << fs->vel.y << " " << fs->vel.z << endl;
	is_data_valid = false;
    }
    if ( !isfinite(fs->tke) ) {
	cout << "Turbulence KE invalid number " << fs->tke << endl;
	is_data_valid = false;
    }
    if ( fs->tke < 0.0 ) {
	cout << "Turbulence KE negative " << fs->tke << endl;
	is_data_valid = false;
    }
    if ( !isfinite(fs->omega) ) {
	cout << "Turbulence frequency invalid number " << fs->omega << endl;
	is_data_valid = false;
    }
    if ( fs->omega <= 0.0 ) {
	cout << "Turbulence frequency nonpositive " << fs->omega << endl;
	is_data_valid = false;
    }
    if ( !is_data_valid ) {
	cout << "cell pos=(" << pos[0].x << "," << pos[0].y << "," << pos[0].z << ")" << endl;
	fs->print();
	cout << "----------------------------------------------------------" << endl;
    }
    return is_data_valid;
} // end of check_flow_data()


/// \brief Compute the time derivatives for the conserved quantities.
///
/// These are the spatial (RHS) terms in the semi-discrete governing equations.
/// \param gtl : (grid-time-level) flow derivatives are evaluated at this grid level
/// \param ftl : (flow-time-level) specifies where computed derivatives are to be stored.
///              0: Start of stage-1 update.
///              1: End of stage-1.
///              2: End of stage-2.
/// \param dimensions : number of space dimensions (2 or 3)
int FV_Cell::time_derivatives(size_t gtl, size_t ftl, size_t dimensions, bool with_k_omega)
{
    global_data &G = *get_global_data_ptr();
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
    double vol_inv = 1.0 / volume[gtl];
    double integral;

    if (status == MASKED_CELL || status == SHADOWED_CELL) {
	dUdt[ftl]->clear_values();
    	return SUCCESS; // do not update cell if it is covered with a piston
    }
    
    // Time-derivative for Mass/unit volume.
    // Note that the unit normals for the interfaces are oriented
    // such that the unit normals for the east, north and top faces
    // are outward and the unit normals for the south, west and
    // bottom faces are inward.
    integral = -IFe->F->mass * IFe->area[gtl] - IFn->F->mass * IFn->area[gtl]
	+ IFw->F->mass * IFw->area[gtl] + IFs->F->mass * IFs->area[gtl];
    if ( dimensions == 3 )
	integral += IFb->F->mass * IFb->area[gtl] - IFt->F->mass * IFt->area[gtl];
    dUdt[ftl]->mass = vol_inv * integral + Q->mass;

    // Time-derivative for X-Momentum/unit volume.
    integral = -IFe->F->momentum.x * IFe->area[gtl] - IFn->F->momentum.x * IFn->area[gtl]
	+ IFw->F->momentum.x * IFw->area[gtl] + IFs->F->momentum.x * IFs->area[gtl];
    if ( dimensions == 3 )
	integral += IFb->F->momentum.x * IFb->area[gtl] - IFt->F->momentum.x * IFt->area[gtl];
    dUdt[ftl]->momentum.x = vol_inv * integral + Q->momentum.x;
    // Time-derivative for Y-Momentum/unit volume.
    integral = -IFe->F->momentum.y * IFe->area[gtl] - IFn->F->momentum.y * IFn->area[gtl]
	+ IFw->F->momentum.y * IFw->area[gtl] + IFs->F->momentum.y * IFs->area[gtl];
    if ( dimensions == 3 )
	integral += IFb->F->momentum.y * IFb->area[gtl] - IFt->F->momentum.y * IFt->area[gtl];
    dUdt[ftl]->momentum.y = vol_inv * integral + Q->momentum.y;
    
    // we require the z-momentum for MHD even in 2D
    if ((dimensions == 3) || ( G.MHD )) {
	// Time-derivative for Z-Momentum/unit volume.
	integral = -IFe->F->momentum.z * IFe->area[gtl] - IFn->F->momentum.z * IFn->area[gtl]
	    + IFw->F->momentum.z * IFw->area[gtl] + IFs->F->momentum.z * IFs->area[gtl];
    }
    if ( dimensions == 3) {
	integral += IFb->F->momentum.z * IFb->area[gtl] - IFt->F->momentum.z * IFt->area[gtl];
    }
    if ((dimensions == 3) || ( G.MHD )) {
	dUdt[ftl]->momentum.z = vol_inv * integral + Q->momentum.z;
    } else {
	dUdt[ftl]->momentum.z = 0.0;
    }
    
    if ( G.MHD ) {
	// Time-derivative for X-Magnetic Field/unit volume.
	integral = -IFe->F->B.x * IFe->area[gtl] - IFn->F->B.x * IFn->area[gtl]
	    + IFw->F->B.x * IFw->area[gtl] + IFs->F->B.x * IFs->area[gtl];
	if ( dimensions == 3 )
	    integral += IFb->F->B.x * IFb->area[gtl] - IFt->F->B.x * IFt->area[gtl];
	dUdt[ftl]->B.x = vol_inv * integral + Q->B.x;
	// Time-derivative for Y-Magnetic Field/unit volume.
	integral = -IFe->F->B.y * IFe->area[gtl] - IFn->F->B.y * IFn->area[gtl]
	    + IFw->F->B.y * IFw->area[gtl] + IFs->F->B.y * IFs->area[gtl];
	if ( dimensions == 3 )
	    integral += IFb->F->B.y * IFb->area[gtl] - IFt->F->B.y * IFt->area[gtl];
	dUdt[ftl]->B.y = vol_inv * integral + Q->B.y;
	// Time-derivative for Z-Magnetic Field/unit volume.
	integral = -IFe->F->B.z * IFe->area[gtl] - IFn->F->B.z * IFn->area[gtl]
	    + IFw->F->B.z * IFw->area[gtl] + IFs->F->B.z * IFs->area[gtl];
	if ( dimensions == 3 ) {
	    integral += IFb->F->B.z * IFb->area[gtl] - IFt->F->B.z * IFt->area[gtl];
	}
	dUdt[ftl]->B.z = vol_inv * integral + Q->B.z;
    }
    else {
	dUdt[ftl]->B.x = 0.0;
	dUdt[ftl]->B.y = 0.0;
	dUdt[ftl]->B.z = 0.0;
    }

    // Time-derivative for Total Energy/unit volume.
    integral = -IFe->F->total_energy * IFe->area[gtl] - IFn->F->total_energy * IFn->area[gtl]
	+ IFw->F->total_energy * IFw->area[gtl] + IFs->F->total_energy * IFs->area[gtl];
    if ( dimensions == 3 )
	integral += IFb->F->total_energy * IFb->area[gtl] - IFt->F->total_energy * IFt->area[gtl];
    dUdt[ftl]->total_energy = vol_inv * integral + Q->total_energy;
    
    if ( with_k_omega ) {
	integral = -IFe->F->tke * IFe->area[gtl] - IFn->F->tke * IFn->area[gtl]
	    + IFw->F->tke * IFw->area[gtl] + IFs->F->tke * IFs->area[gtl];
	if ( dimensions == 3 )
	    integral += IFb->F->tke * IFb->area[gtl] - IFt->F->tke * IFt->area[gtl];
	dUdt[ftl]->tke = vol_inv * integral + Q->tke;
	
	integral = -IFe->F->omega * IFe->area[gtl] - IFn->F->omega * IFn->area[gtl]
	    + IFw->F->omega * IFw->area[gtl] + IFs->F->omega * IFs->area[gtl];
	if ( dimensions == 3 )
	    integral += IFb->F->omega * IFb->area[gtl] - IFt->F->omega * IFt->area[gtl];
	dUdt[ftl]->omega = vol_inv * integral + Q->omega;
    } else {
	dUdt[ftl]->tke = 0.0;
	dUdt[ftl]->omega = 0.0;
    }
    // Time-derivative for individual species.
    // The conserved quantity is the mass per unit
    // volume of species isp and
    // the fluxes are mass/unit-time/unit-area.
    // Units of DmassfDt are 1/sec.
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	integral =
	    -IFe->F->massf[isp] * IFe->area[gtl]
	    - IFn->F->massf[isp] * IFn->area[gtl]
	    + IFw->F->massf[isp] * IFw->area[gtl]
	    + IFs->F->massf[isp] * IFs->area[gtl];
	if ( dimensions == 3 )
	    integral += IFb->F->massf[isp] * IFb->area[gtl] - IFt->F->massf[isp] * IFt->area[gtl];
	dUdt[ftl]->massf[isp] = vol_inv * integral + Q->massf[isp];
    }
    // Individual energies.
    // We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
    // Instead we get this from the conservation of total energy
    for ( size_t imode = 1; imode < nmodes; ++imode ) {
	integral =
	    -IFe->F->energies[imode] * IFe->area[gtl]
	    - IFn->F->energies[imode] * IFn->area[gtl]
	    + IFw->F->energies[imode] * IFw->area[gtl]
	    + IFs->F->energies[imode] * IFs->area[gtl];
	if ( dimensions == 3 )
	    integral += IFb->F->energies[imode] * IFb->area[gtl] - IFt->F->energies[imode] * IFt->area[gtl];
	dUdt[ftl]->energies[imode] = vol_inv * integral + Q->energies[imode];
    }
    return SUCCESS;
} // end of time_derivatives()


int FV_Cell::stage_1_update_for_flow_on_fixed_grid(double dt, bool force_euler, bool with_k_omega)
{
    global_data &G = *get_global_data_ptr();
    ConservedQuantities &dUdt0 = *(dUdt[0]);
    ConservedQuantities &U0 = *(U[0]);
    ConservedQuantities &U1 = *(U[1]);
    double gamma_1 = 1.0; // for normal Predictor-Corrector or Euler update.
    // In some parts of the code (viscous updates, k-omega updates)
    // we use this function as an Euler update even when the main
    // gasdynamic_update_scheme is of higher order.
    if ( !force_euler ) {
	switch ( get_gasdynamic_update_scheme() ) {
	case EULER_UPDATE:
	case PC_UPDATE: gamma_1 = 1.0; break;
	case MIDPOINT_UPDATE: gamma_1 = 0.5; break;
	case CLASSIC_RK3_UPDATE: gamma_1 = 0.5; break;
	case TVD_RK3_UPDATE: gamma_1 = 1.0; break;
	case DENMAN_RK3_UPDATE: gamma_1 = 8.0/15.0; break;
	default:
	    throw std::runtime_error("FV_Cell::stage_1_update_for_flow_on_fixed_grid(): "
				     "should not be here!");
	}
    }
    U1.mass = U0.mass + dt * gamma_1 * dUdt0.mass;
    // Side note: 
    // It would be convenient (codewise) for the updates of these Vector3 quantities to
    // be done with the Vector3 arithmetic operators but I suspect that the implementation
    // of those oerators is such that a whole lot of Vector3 temporaries would be created.
    U1.momentum.x = U0.momentum.x + dt * gamma_1 * dUdt0.momentum.x;
    U1.momentum.y = U0.momentum.y + dt * gamma_1 * dUdt0.momentum.y;
    U1.momentum.z = U0.momentum.z + dt * gamma_1 * dUdt0.momentum.z;
    if ( G.MHD ) {
	// Magnetic field
	U1.B.x = U0.B.x + dt * gamma_1 * dUdt0.B.x;
	U1.B.y = U0.B.y + dt * gamma_1 * dUdt0.B.y;
	U1.B.z = U0.B.z + dt * gamma_1 * dUdt0.B.z;
    }
    U1.total_energy = U0.total_energy + dt * gamma_1 * dUdt0.total_energy;
    if ( with_k_omega ) {
	U1.tke = U0.tke + dt * gamma_1 * dUdt0.tke;
	U1.tke = max(U1.tke, 0.0);
	U1.omega = U0.omega + dt * gamma_1 * dUdt0.omega;
	U1.omega = max(U1.omega, U0.mass);
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
	U1.tke = U0.tke;
	U1.omega = U0.omega;
    }
    for ( size_t isp = 0; isp < U1.massf.size(); ++isp ) {
	U1.massf[isp] = U0.massf[isp] + dt * gamma_1 * dUdt0.massf[isp];
    }
    // We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
    // Instead we get this from the conservation of total energy
    for ( size_t imode = 1; imode < U1.energies.size(); ++imode ) {
	U1.energies[imode] = U0.energies[imode] + dt * gamma_1 * dUdt0.energies[imode];
    }
    return SUCCESS;
} // end of stage_1_update_for_flow_on_fixed_grid()


int FV_Cell::stage_2_update_for_flow_on_fixed_grid(double dt, bool with_k_omega)
{
    global_data &G = *get_global_data_ptr();
    ConservedQuantities &dUdt0 = *(dUdt[0]);
    ConservedQuantities &dUdt1 = *(dUdt[1]);
    ConservedQuantities *U_old = U[0];
    if ( get_gasdynamic_update_scheme() == DENMAN_RK3_UPDATE ) U_old = U[1];
    ConservedQuantities &U2 = *(U[2]);
    double gamma_1 = 0.5; // Presume predictor-corrector.
    double gamma_2 = 0.5;
    switch ( get_gasdynamic_update_scheme() ) {
    case PC_UPDATE: gamma_1 = 0.5, gamma_2 = 0.5; break;
    case MIDPOINT_UPDATE: gamma_1 = 0.0; gamma_2 = 1.0; break;
    case CLASSIC_RK3_UPDATE: gamma_1 = -1.0; gamma_2 = 2.0; break;
    case TVD_RK3_UPDATE: gamma_1 = 0.25; gamma_2 = 0.25; break;
    case DENMAN_RK3_UPDATE: gamma_1 = -17.0/60.0; gamma_2 = 5.0/12.0; break;
    default:
	throw std::runtime_error("FV_Cell::stage_2_update_for_flow_on_fixed_grid(): "
				 "should not be here!");
    }
    U2.mass = U_old->mass + dt * (gamma_1 * dUdt0.mass + gamma_2 * dUdt1.mass);
    U2.momentum.x = U_old->momentum.x + dt * (gamma_1 * dUdt0.momentum.x + gamma_2 * dUdt1.momentum.x);
    U2.momentum.y = U_old->momentum.y + dt * (gamma_1 * dUdt0.momentum.y + gamma_2 * dUdt1.momentum.y);
    U2.momentum.z = U_old->momentum.z + dt * (gamma_1 * dUdt0.momentum.z + gamma_2 * dUdt1.momentum.z);
    if ( G.MHD ) {
	// Magnetic field
	U2.B.x = U_old->B.x + dt * (gamma_1 * dUdt0.B.x + gamma_2 * dUdt1.B.x);
	U2.B.y = U_old->B.y + dt * (gamma_1 * dUdt0.B.y + gamma_2 * dUdt1.B.y);
	U2.B.z = U_old->B.z + dt * (gamma_1 * dUdt0.B.z + gamma_2 * dUdt1.B.z);
    }
    U2.total_energy = U_old->total_energy + 
	dt * (gamma_1 * dUdt0.total_energy + gamma_2 * dUdt1.total_energy);
    if ( with_k_omega ) {
	U2.tke = U_old->tke + dt * (gamma_1 * dUdt0.tke + gamma_2 * dUdt1.tke);
	U2.tke = max(U2.tke, 0.0);
	U2.omega = U_old->omega + dt * (gamma_1 * dUdt0.omega + gamma_2 * dUdt1.omega);
	U2.omega = max(U2.omega, U_old->mass);
    } else {
	U2.tke = U_old->tke;
	U2.omega = U_old->omega;
    }
    for ( size_t isp = 0; isp < U2.massf.size(); ++isp ) {
	U2.massf[isp] = U_old->massf[isp] + dt * (gamma_1 * dUdt0.massf[isp] + gamma_2 * dUdt1.massf[isp]);
    }
    // We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
    // Instead we get this from the conservation of total energy
    for ( size_t imode = 1; imode < U2.energies.size(); ++imode ) {
	U2.energies[imode] = U_old->energies[imode] + 
	    dt * (gamma_1 * dUdt0.energies[imode] + gamma_2 * dUdt1.energies[imode]);
    }
    return SUCCESS;
} // end of stage_2_update_for_flow_on_fixed_grid()


int FV_Cell::stage_3_update_for_flow_on_fixed_grid(double dt, bool with_k_omega)
{
    global_data &G = *get_global_data_ptr();
    ConservedQuantities &dUdt0 = *(dUdt[0]);
    ConservedQuantities &dUdt1 = *(dUdt[1]);
    ConservedQuantities &dUdt2 = *(dUdt[2]);
    ConservedQuantities *U_old = U[0];
    if ( get_gasdynamic_update_scheme() == DENMAN_RK3_UPDATE ) U_old = U[2];
    ConservedQuantities &U3 = *(U[3]);
    double gamma_1 = 1.0/6.0; // presume TVD_RK3 scheme.
    double gamma_2 = 1.0/6.0;
    double gamma_3 = 4.0/6.0;
    switch ( get_gasdynamic_update_scheme() ) {
    case CLASSIC_RK3_UPDATE: gamma_1 = 1.0/6.0; gamma_2 = 4.0/6.0; gamma_3 = 1.0/6.0; break;
    case TVD_RK3_UPDATE: gamma_1 = 1.0/6.0; gamma_2 = 1.0/6.0; gamma_3 = 4.0/6.0; break;
    // FIX-ME: Really don't think that we have Andrew Denman's scheme ported correctly.
    case DENMAN_RK3_UPDATE: gamma_1 = 0.0; gamma_2 = -5.0/12.0; gamma_3 = 3.0/4.0; break;
    default:
	throw std::runtime_error("FV_Cell::stage_3_update_for_flow_on_fixed_grid(): "
				 "should not be here!");
    }
    U3.mass = U_old->mass + dt * (gamma_1*dUdt0.mass + gamma_2*dUdt1.mass + gamma_3*dUdt2.mass);
    U3.momentum.x = U_old->momentum.x +
	dt * (gamma_1*dUdt0.momentum.x + gamma_2*dUdt1.momentum.x + gamma_3*dUdt2.momentum.x);
    U3.momentum.y = U_old->momentum.y +
	dt * (gamma_1*dUdt0.momentum.y + gamma_2*dUdt1.momentum.y + gamma_3*dUdt2.momentum.y);
    U3.momentum.z = U_old->momentum.z + 
	dt * (gamma_1*dUdt0.momentum.z + gamma_2*dUdt1.momentum.z + gamma_3*dUdt2.momentum.z);
    if ( G.MHD ) {
	// Magnetic field
	U3.B.x = U_old->B.x + dt * (gamma_1*dUdt0.B.x + gamma_2*dUdt1.B.x + gamma_3*dUdt2.B.x);
	U3.B.y = U_old->B.y + dt * (gamma_1*dUdt0.B.y + gamma_2*dUdt1.B.y + gamma_3*dUdt2.B.y);
	U3.B.z = U_old->B.z + dt * (gamma_1*dUdt0.B.z + gamma_2*dUdt1.B.z + gamma_3*dUdt2.B.z);
    }
    U3.total_energy = U_old->total_energy + 
	dt * (gamma_1*dUdt0.total_energy + gamma_2*dUdt1.total_energy + gamma_3*dUdt2.total_energy);
    if ( with_k_omega ) {
	U3.tke = U_old->tke + dt * (gamma_1*dUdt0.tke + gamma_2*dUdt1.tke + gamma_3*dUdt2.tke);
	U3.tke = max(U3.tke, 0.0);
	U3.omega = U_old->omega + dt * (gamma_1*dUdt0.omega + gamma_2*dUdt1.omega + gamma_3*dUdt2.omega);
	U3.omega = max(U3.omega, U_old->mass);
    } else {
	U3.tke = U_old->tke;
	U3.omega = U_old->omega;
    }
    for ( size_t isp = 0; isp < U3.massf.size(); ++isp ) {
	U3.massf[isp] = U_old->massf[isp] +
	    dt * (gamma_1*dUdt0.massf[isp] + gamma_2*dUdt1.massf[isp] + gamma_3*dUdt2.massf[isp]);
    }
    // We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
    // Instead we get this from the conservation of total energy
    for ( size_t imode = 1; imode < U3.energies.size(); ++imode ) {
	U3.energies[imode] = U_old->energies[imode] +
	    dt * (gamma_1*dUdt0.energies[imode] + gamma_2*dUdt1.energies[imode] +
		  gamma_3*dUdt2.energies[imode]);
    }
    return SUCCESS;
} // end of stage_3_update_for_flow_on_fixed_grid()


int FV_Cell::stage_1_update_for_flow_on_moving_grid(double dt, bool with_k_omega)
{
    global_data &G = *get_global_data_ptr();
    ConservedQuantities &dUdt0 = *(dUdt[0]);
    ConservedQuantities &U0 = *(U[0]);
    ConservedQuantities &U1 = *(U[1]);
    double gamma_1 = 1.0;
    double vr = volume[0] / volume[1];

    U1.mass = vr * (U0.mass + dt * gamma_1 * dUdt0.mass);
    U1.momentum.x = vr * (U0.momentum.x + dt * gamma_1 * dUdt0.momentum.x);
    U1.momentum.y = vr * (U0.momentum.y + dt * gamma_1 * dUdt0.momentum.y);
    U1.momentum.z = vr * (U0.momentum.z + dt * gamma_1 * dUdt0.momentum.z);
    if ( G.MHD ) {
	// Magnetic field
	U1.B.x = vr * (U0.B.x + dt * gamma_1 * dUdt0.B.x);
	U1.B.y = vr * (U0.B.y + dt * gamma_1 * dUdt0.B.y);
	U1.B.z = vr * (U0.B.z + dt * gamma_1 * dUdt0.B.z);
    }
    U1.total_energy = vr * (U0.total_energy + dt * gamma_1 * dUdt0.total_energy);
    if ( with_k_omega ) {
	U1.tke = vr * (U0.tke + dt * gamma_1 * dUdt0.tke);
	U1.tke = max(U1.tke, 0.0);
	U1.omega = vr * (U0.omega + dt * gamma_1 * dUdt0.omega);
	U1.omega = max(U1.omega, U0.mass);
    } else {
	U1.tke = U0.tke;
	U1.omega = U0.omega;
    }
    for ( size_t isp = 0; isp < U1.massf.size(); ++isp ) {
	U1.massf[isp] = vr * (U0.massf[isp] + dt * gamma_1 * dUdt0.massf[isp]);
    }
    // We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
    // Instead we get this from the conservation of total energy
    for ( size_t imode = 1; imode < U1.energies.size(); ++imode ) {
	U1.energies[imode] = vr * (U0.energies[imode] + dt * gamma_1 * dUdt0.energies[imode]);
    }
    return SUCCESS;
} // end of stage_1_update_for_flow_on_moving_grid()


int FV_Cell::stage_2_update_for_flow_on_moving_grid(double dt, bool with_k_omega)
{
    global_data &G = *get_global_data_ptr();
    ConservedQuantities &dUdt0 = *(dUdt[0]);
    ConservedQuantities &dUdt1 = *(dUdt[1]);
    ConservedQuantities &U0 = *(U[0]);
    // ConservedQuantities &U1 = *(U[1]);
    ConservedQuantities &U2 = *(U[2]);
    double gamma_2 = 0.5;
    double gamma_1 = 0.5;
    double v_old = volume[0];
    double vol_inv = 1.0 / volume[2];
    gamma_1 *= volume[0]; gamma_2 *= volume[1]; // Roll-in the volumes for convenience below. 
    
    U2.mass = vol_inv * (v_old * U0.mass + dt * (gamma_1 * dUdt0.mass + gamma_2 * dUdt1.mass));
    U2.momentum.x = vol_inv * (v_old * U0.momentum.x + 
			       dt * (gamma_1 * dUdt0.momentum.x + gamma_2 * dUdt1.momentum.x));
    U2.momentum.y = vol_inv * (v_old * U0.momentum.y + 
			       dt * (gamma_1 * dUdt0.momentum.y + gamma_2 * dUdt1.momentum.y));
    U2.momentum.z = vol_inv * (v_old * U0.momentum.z + 
			       dt * (gamma_1 * dUdt0.momentum.z + gamma_2 * dUdt1.momentum.z));
    if ( G.MHD ) {
	// Magnetic field
	U2.B.x = vol_inv * (v_old * U0.B.x + dt * (gamma_1 * dUdt0.B.x + gamma_2 * dUdt1.B.x));
	U2.B.y = vol_inv * (v_old * U0.B.y + dt * (gamma_1 * dUdt0.B.y + gamma_2 * dUdt1.B.y));
	U2.B.z = vol_inv * (v_old * U0.B.z + dt * (gamma_1 * dUdt0.B.z + gamma_2 * dUdt1.B.z));
    }
    U2.total_energy = vol_inv * (v_old * U0.total_energy + 
				 dt * (gamma_1 * dUdt0.total_energy + gamma_2 * dUdt1.total_energy));
    if ( with_k_omega ) {
	U2.tke = vol_inv * (v_old * U0.tke + dt * (gamma_1 * dUdt0.tke + gamma_2 * dUdt1.tke));
	U2.tke = max(U2.tke, 0.0);
	U2.omega = vol_inv * (v_old * U0.omega + dt * (gamma_1 * dUdt0.omega + gamma_2 * dUdt1.omega));
	U2.omega = max(U2.omega, U0.mass);
    } else {
	U2.tke = vol_inv * (v_old * U0.tke);
	U2.omega = vol_inv * (v_old * U0.omega);
    }
    for ( size_t isp = 0; isp < U2.massf.size(); ++isp ) {
	U2.massf[isp] = vol_inv * (v_old * U0.massf[isp] +
				   dt * (gamma_1 * dUdt0.massf[isp] + gamma_2 * dUdt1.massf[isp]));
    }
    // We will not put anything meaningful in imode = 0 (RJG & DFP : 22-Apr-2013)
    // Instead we get this from the conservation of total energy
    for ( size_t imode = 1; imode < U2.energies.size(); ++imode ) {
	U2.energies[imode] = vol_inv * (v_old * U0.energies[imode] +
					dt * (gamma_1 * dUdt0.energies[imode] + gamma_2 * dUdt1.energies[imode]));
    }
    return SUCCESS;
} // end of stage_2_update_for_flow_on_moving_grid()


/// \brief Apply the chemistry update for a specified cell.
///
/// Use the finite-rate chemistry module to update the
/// species fractions and the other thermochemical properties.
/// \param dt   : size of the time-step
int FV_Cell::chemical_increment(double dt, double T_frozen)
{
    if ( !fr_reactions_allowed or fs->gas->T[0] <= T_frozen ) return SUCCESS;
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    Reaction_update *rupdate = get_reaction_update_ptr();
    constexpr bool copy_gas_in_case_of_failure = false;
    Gas_data *gcopy;
    if ( copy_gas_in_case_of_failure ) {
	// Make a copy so that we can print out if things go wrong.
	gcopy = new Gas_data(*(fs->gas));
    }
    double T_save = fs->gas->T[0];
    if ( G.ignition_zone_active ) {
	// When active, replace gas temperature with an effective ignition temperature
	for ( CIgnitionZone &iz : G.ignition_zone ) {
	    if ( pos[0].x >= iz.x0 && pos[0].x <= iz.x1 && pos[0].y >= iz.y0 && pos[0].y <= iz.y1 &&
		 (G.dimensions == 2 || (pos[0].z >= iz.z0 && pos[0].z <= iz.z1)) ) {
		fs->gas->T[0] = iz.Tig;
	    }
	}
    }
    int flag = rupdate->update_state(*(fs->gas), dt, dt_chem, gmodel);
    if ( G.ignition_zone_active ) {
	// Restore actual gas temperature
	fs->gas->T[0] = T_save;
    }
    if ( flag != SUCCESS ) {
	cout << "The chemical_increment() failed for cell: " << id << endl;
	if ( copy_gas_in_case_of_failure ) {
	    cout << "The gas state before the update was:\n";
	    gcopy->print_values();
	    delete gcopy;
	}
	cout << "The gas state after the update was:\n";
	fs->gas->print_values();
	return flag;
    }

    // The update only changes mass fractions, we need to impose
    // a thermodynamic constraint based on a call to the equation
    // of state.
    gmodel->eval_thermo_state_rhoe(*(fs->gas));

    // If we are doing a viscous sim, we'll need to ensure
    // viscous properties are up-to-date
    if ( G.viscous ) gmodel->eval_transport_coefficients(*(fs->gas));
    if ( G.diffusion ) gmodel->eval_diffusion_coefficients(*(fs->gas));
    // ...but we have to manually update the conservation quantities
    // for the gas-dynamics time integration.
    // Species densities: mass of species isp per unit volume.
    for ( size_t isp = 0; isp < fs->gas->massf.size(); ++isp )
	U[0]->massf[isp] = fs->gas->rho * fs->gas->massf[isp];
    return flag;
} // end of chemical_increment()


/// \brief Apply the thermal update for a specified cell.
///
/// Use the nonequilibrium multi-Temperature module to update the
/// energy values and the other thermochemical properties.
/// We are assuming that this is done after a successful gas-dynamic update
/// and that the current conserved quantities are held in U[0].
///
/// \param dt   : size of the time-step
int FV_Cell::thermal_increment(double dt, double T_frozen_energy)
{
    if ( !fr_reactions_allowed or fs->gas->T[0] <= T_frozen_energy ) return SUCCESS;
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    Energy_exchange_update *eeupdate = get_energy_exchange_update_ptr();

    int flag = eeupdate->update_state(*(fs->gas), dt, dt_therm, gmodel);

    // The update only changes modal energies, we need to impose
    // a thermodynamic constraint based on a call to the equation
    // of state.
    gmodel->eval_thermo_state_rhoe(*(fs->gas));

    // If we are doing a viscous sim, we'll need to ensure
    // viscous properties are up-to-date
    if ( G.viscous ) gmodel->eval_transport_coefficients(*(fs->gas));
    if ( G.diffusion ) gmodel->eval_diffusion_coefficients(*(fs->gas));
    // ...but we have to manually update the conservation quantities
    // for the gas-dynamics time integration.
    // Independent energies energy: Joules per unit volume.
    for ( size_t imode = 0; imode < U[0]->energies.size(); ++imode ) {
	U[0]->energies[imode] = fs->gas->rho * fs->gas->e[imode];
    }
    return flag;
} // end of thermal_increment()


/// \brief Compute the signal frequency (1/seconds) for a cell.
///
/// \param dimensions: number of spatial dimensions: 2 for mbcns2, 3 for eilmer
/// The North and East faces are taken as the representative lengths the cells.
double FV_Cell::signal_frequency(size_t dimensions, bool with_k_omega)
{
    double signal;
    double un_N, un_E, un_T, u_mag;
    double Bn_N = 0.0;
    double Bn_E = 0.0;
    double Bn_T = 0.0;
    double B_mag = 0.0;
    double ca2 = 0.0;
    double cfast = 0.0;
    double gam_eff;
    int statusf;
    global_data &G = *get_global_data_ptr();
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
	u_mag = sqrt(fs->vel.x*fs->vel.x + fs->vel.y*fs->vel.y + fs->vel.z*fs->vel.z);
    }  else {
	un_T = 0.0;
	u_mag = sqrt(fs->vel.x*fs->vel.x + fs->vel.y*fs->vel.y);
    }
    if ( G.MHD ) {
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
    if ( G.stringent_cfl ) {
	// Make the worst case.
	if ( G.MHD ) {
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
	if ( G.MHD ) {
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
	    } else {
		signalN = (un_N + cfast) / jLength;
		signalE = (un_E + cfast) / iLength;
		signal = max(signalN, signalE);
	    }
	} else if ( dimensions == 3 ) {
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
	    signalN = (un_N + fs->gas->a) / jLength;
	    signalE = (un_E + fs->gas->a) / iLength;
	    signal = max(signalN, signalE);
	}
    }
    if ( G.viscous && G.implicit_mode == 0 && fs->gas->mu > 10.0e-23) {
	// Factor for the viscous time limit.
	// This factor is not included if viscosity is zero.
	const int VISCOUS_TIME_LIMIT_MODEL = 0; // ==0 original Swanson model,
                                                // ==1 Ramshaw model
	if ( VISCOUS_TIME_LIMIT_MODEL == 0 ) {
	    // See Swanson, Turkel and White (1991)
	    gam_eff = gmodel->gamma(*(fs->gas), statusf);
	    // Need to sum conductivities for TNE
	    double k_total = 0.0;
	    for ( size_t i=0; i<fs->gas->k.size(); ++i ) k_total += fs->gas->k[i];
	    double Prandtl = fs->gas->mu * gmodel->Cp(*(fs->gas), statusf) / k_total;
	    if ( dimensions == 3 ) {
		signal += 4.0 * G.viscous_factor * (fs->gas->mu + fs->mu_t)
		    * gam_eff / (Prandtl * fs->gas->rho)
		    * (1.0/(iLength*iLength) + 1.0/(jLength*jLength) + 1.0/(kLength*kLength));
	    } else {
		signal += 4.0 * G.viscous_factor * (fs->gas->mu + fs->mu_t) 
		    * gam_eff / (Prandtl * fs->gas->rho)
		    * (1.0/(iLength*iLength) + 1.0/(jLength*jLength));
	    }
	} else if ( VISCOUS_TIME_LIMIT_MODEL == 1 ) {
	    // A viscous time limit model incorporating diffusion effects
	    // See Ramshaw and Chang PCPP V.12 n.3 1992 p314
	    // 1. Viscosity signal frequency
	    double D_a = 4.0 * G.viscous_factor * (fs->gas->mu + -0.66667 * fs->gas->mu +
						   fs->mu_t - 0.66667 * fs->mu_t ) / fs->gas->rho;
	    // 1. Conductivity signal frequency
	    double D_b = fs->gas->k[0];
	    size_t nmodes = gmodel->get_number_of_modes();
	    for ( size_t imode=1; imode < nmodes; ++imode )
		D_b += fs->gas->k[imode];
	    double c_v = gmodel->Cv(*(fs->gas), statusf);
	    D_b *= G.viscous_factor / ( fs->gas->rho * c_v );
	    // 3. calculate the largest mixture diffusivity
	    double D_c = 0.0;
	    if ( G.diffusion ) {
		size_t nsp = gmodel->get_number_of_species();
		vector<double> x(nsp);
		vector<double> M(nsp);
		vector<double> DAV_im(nsp);
		for ( size_t isp=0; isp<nsp; ++isp )
		    M[isp] = gmodel->molecular_weight(isp);
		// [todo] FIX-ME Rowan, please
		// fill_in_x and fill_in_DAV_im don't seem to be available.
		// fill_in_x(fs->gas->rho, fs->gas->T, fs->gas->massf, M, x);
		// fill_in_DAV_im(fs->gas->D_AB, x, DAV_im);
		throw std::runtime_error("Ramshaw and Chang viscous time-limit is broken.");
		for ( size_t isp=0; isp<nsp; ++isp )
		    if ( DAV_im[isp] > D_c ) D_c = DAV_im[isp];
	    }
	    D_c *= G.diffusion_factor;
	    // 4. Find the maximum effective diffusion
	    double D_max_i = max( D_a, D_b );
	    double D_max_ii = max( D_b, D_c );
	    double D_max = max( D_max_i, D_max_ii );
	    // 5. Add signal frequency contribution
	    if ( dimensions == 3 ) {
		signal += D_max * (1.0/(iLength*iLength)
				   + 1.0/(jLength*jLength)
				   + 1.0/(kLength*kLength));
	    } else {
		signal += D_max * (1.0/(iLength*iLength)
				   + 1.0/(jLength*jLength));
	    }
	    // end if VISCOUS_TIME_LIMIT_MODEL == 1
	} else {
	    throw std::runtime_error("Viscous effects are on but a viscous"
				     " time-step-limit model is not active.");
	} // end if VISCOUS_TIME_LIMIT_MODEL
    }
    if ( with_k_omega == 1 ) {
	if ( fs->omega > signal ) signal = fs->omega;
    }
    return signal;
} // end of signal_frequency()


int FV_Cell::turbulence_viscosity_zero()
{
    fs->mu_t = 0.0;
    fs->k_t = 0.0;
    return SUCCESS;
}


int FV_Cell::turbulence_viscosity_zero_if_not_in_zone()
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
    fs->mu_t = min(fs->mu_t, factor * fs->gas->mu);
    fs->k_t = min(fs->k_t , factor * fs->gas->k[0]); // ASSUMPTION re k[0]
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
int FV_Cell::turbulence_viscosity_k_omega()
{
    global_data &G = *get_global_data_ptr();
    if ( G.turbulence_model != TM_K_OMEGA ) {
	// FIX-ME may have to do something better if another turbulence model is active.
	fs->mu_t = 0.0;
	fs->k_t = 0.0;
	return SUCCESS;
    }
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
        if ( G.axisymmetric ) {
            // 2D axisymmetric
            double v_over_y = fs->vel.y / pos[0].y;
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
    S_bar_squared = max(0.0, S_bar_squared);
    double omega_t = max(fs->omega, C_lim*sqrt(2.0*S_bar_squared/beta_star));
    fs->mu_t = fs->gas->rho * fs->tke / omega_t;
    double Pr_t = G.turbulence_prandtl;
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
    if ( !in_turbulent_zone ) return SUCCESS;

    constexpr bool KOMEGA_IMPLICIT_UPDATE = true;
    if ( KOMEGA_IMPLICIT_UPDATE ) {
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
	U[0]->tke = fs->gas->rho * fs->tke;
	U[0]->omega = fs->gas->rho * fs->omega;

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
		    U[0]->tke = fs->gas->rho * tke_updated;
		} else {
		    // Next estimate for rtke.
		    U[0]->tke = delta_rtke + fs->gas->rho * tke_updated;
		}
		if (delta_romega + fs->gas->rho * omega_updated < 0.0) {
		    // Don't let romega go negative.
		    U[0]->omega = fs->gas->rho * omega_updated;
		} else {
		    // Next estimate for romega.
		    U[0]->omega = delta_romega + fs->gas->rho * omega_updated;
		}
		// Decode for the next step of the Newton-solve loop
		tke_updated = U[0]->tke / fs->gas->rho;
		omega_updated = U[0]->omega / fs->gas->rho;
	    }
	} // End of Newton-solve loop for implicit update scheme
    } else {
	// Explicit updating scheme.
	double DrtkeDt, DromegaDt, rtke_increment, romega_increment;
	int n_little_steps = 20;  // Make this as large as needed to get stability.
	double dt_little = dt / n_little_steps;
    
	// Encode conserved quantities for cell.
	U[0]->tke = fs->gas->rho * fs->tke;
	U[0]->omega = fs->gas->rho * fs->omega;
	for ( int i = 1; i <= n_little_steps; ++i ) {
	    this->k_omega_time_derivatives(&DrtkeDt, &DromegaDt, fs->tke, fs->omega);
	    rtke_increment = dt_little * DrtkeDt;
	    romega_increment = dt_little * DromegaDt;
	    if ( U[0]->tke + rtke_increment < 0.0 ||
		 (rtke_increment > 0.0 && U[0]->tke + rtke_increment > 0.5 * U[0]->total_energy) ) {
		// Don't let rtke go negative and don't let it grow too large.
		rtke_increment = 0.0;
	    }
	    if ( U[0]->omega + romega_increment < 0.0 ) {
		// Don't let romega go negative.
		romega_increment = 0.0;
            }
	    U[0]->tke += rtke_increment;
	    U[0]->omega += romega_increment;
	    // Decode conserved quantities.
	    fs->tke = U[0]->tke / fs->gas->rho;
	    fs->omega = U[0]->omega / fs->gas->rho;
	} // End of for-loop for explicit update scheme
    } // end if

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
    global_data &G = *get_global_data_ptr();
    if ( G.turbulence_model != TM_K_OMEGA ) {
	// FIX-ME may need to do something better is another turbulence model is active.
	*Q_rtke = 0.0;
	*Q_romega = 0.0;
	return SUCCESS;
    }
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
        if ( G.axisymmetric ) {
            // 2D axisymmetric
            double v_over_y = fs->vel.y / pos[0].y;
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
    constexpr double P_OVER_D_LIMIT = 25.0;
    P_K = min(P_K, P_OVER_D_LIMIT*D_K);

    if ( cross_diff > 0 ) sigma_d = 0.125;
    P_W = alpha * omega / max(tke,SMALL_TKE) * P_K +
          sigma_d * fs->gas->rho / max(omega,SMALL_OMEGA) * cross_diff;

    X_w = fabs(WWS / pow(beta_star*omega, 3)) ;
    f_beta = (1.0 + 85.0 * X_w) / (1.0 + 100.0 * X_w) ;
    beta = beta_0 * f_beta;
    D_W = beta * fs->gas->rho * omega * omega;

    *Q_rtke = P_K - D_K;
    *Q_romega = P_W - D_W;
    return SUCCESS;
} // end k_omega_time_derivatives()


int FV_Cell::clear_source_vector()
// When doing the gasdynamic update stages, the source vector values
// are accumulated for the inviscid and then viscous terms, so we 
// have to start with a clean slate, so to speak.
{
    Q->mass = 0.0;
    Q->momentum.x = 0.0;
    Q->momentum.y = 0.0;
    Q->momentum.z = 0.0;
    Q->B.x = 0.0;
    Q->B.y = 0.0;
    Q->B.z = 0.0;
    Q->total_energy = 0.0;
    Q->tke = 0.0;
    Q->omega = 0.0;
    for ( size_t isp = 0; isp < Q->massf.size(); ++isp )
	Q->massf[isp] = 0.0;
    for ( size_t imode = 0; imode < Q->energies.size(); ++imode )
	Q->energies[imode] = 0.0;
    Q_rE_rad = 0.0;
    return SUCCESS;
} // end FV_Cell::clear_source_vector()


/// \brief Add the components of the source vector, Q, for inviscid flow.
///
/// Currently, the axisymmetric equations include the
/// pressure contribution to the y-momentum equation
/// here rather than in the boundary fluxes.
/// By default, assume 2D-planar, or 3D-Cartesian flow.
int FV_Cell::add_inviscid_source_vector(int gtl, double omegaz)
{
    global_data &G = *get_global_data_ptr();
    Q->total_energy += G.heat_factor * base_qdot;
    if ( omegaz != 0.0 ) {
	// Rotating frame.
	double rho = fs->gas->rho;
	double x = pos[gtl].x;
	double y = pos[gtl].y;
        double wx = fs->vel.x;
	double wy = fs->vel.y;
	// Coriolis and centrifugal forces contribute to momenta.
	Q->momentum.x += rho * (omegaz*omegaz*x + 2.0*omegaz*wy);
	Q->momentum.y += rho * (omegaz*omegaz*y - 2.0*omegaz*wx);
	// There is no contribution to the energy equation in the rotating frame
        // because it is implicit in the use of rothalpy as the conserved quantity.
    }
    if ( G.axisymmetric ) {
	// For axisymmetric flow:
	// pressure contribution from the Front and Back (radial) interfaces.
	Q->momentum.y += fs->gas->p * area[gtl] / volume[gtl];
    }
    // Species production (other than chemistry).
    // For the chemistry, see chemical_increment().
    // Individual energies (other than energy exchange)
    // For the energy exchange, see thermal_increment()
    // Radiation can potentially be removed from both the electronic and
    // total energy source terms.
    if ( G.radiation ) {
	// Radiative source term should be already calculated
	// Add value to total energy
	// FIX-ME: - assuming electronic mode is the last in the vector of energies
	//         - what about Q_renergies[0]?
	Q->total_energy += Q_rE_rad;
	Q->energies.back() += Q_rE_rad;
    }
    return SUCCESS;
} // end FV_Cell::add_inviscid_source_vector()


/// \brief Add the components of the source vector, Q, for viscous flow.
int FV_Cell::add_viscous_source_vector(bool with_k_omega)
{
    global_data &G = *get_global_data_ptr();
    if ( G.axisymmetric ) {
	// For viscous, axisymmetric flow:
	double v_over_y = fs->vel.y / pos[0].y;
	double dudx = 0.25 * (vtx[0]->dudx + vtx[1]->dudx + vtx[2]->dudx + vtx[3]->dudx);
	double dvdy = 0.25 * (vtx[0]->dvdy + vtx[1]->dvdy + vtx[2]->dvdy + vtx[3]->dvdy);
	double mu = 0.25 * (iface[EAST]->fs->gas->mu + iface[WEST]->fs->gas->mu +
			    iface[NORTH]->fs->gas->mu + iface[SOUTH]->fs->gas->mu) +
	            0.25 * (iface[EAST]->fs->mu_t + iface[WEST]->fs->mu_t +
			    iface[NORTH]->fs->mu_t + iface[SOUTH]->fs->mu_t);
	mu *= G.viscous_factor;
	double lmbda = -2.0/3.0 * mu;
	double tau_00 = 2.0 * mu * v_over_y + lmbda * (dudx + dvdy + v_over_y);
	// Y-Momentum; viscous stress contribution from the front and Back interfaces.
	// Note that these quantities are approximated at the
	// mid-point of the cell face and so should never be
	// singular -- at least I hope that this is so.
	Q->momentum.y -= tau_00 * area[0] / volume[0];
    } // end if ( G.axisymmetric )

    if ( with_k_omega ) {
	double Q_tke = 0.0; double Q_omega = 0.0;
        if ( in_turbulent_zone ) {
	    this->k_omega_time_derivatives(&Q_tke, &Q_omega, fs->tke, fs->omega);
        }
	Q->tke += Q_tke; Q->omega += Q_omega;
    }

    if ( G.electric_field_work ) {
        // Work done on electrons due to electric field induced by charge separation
	// on scales less than the Debye length
	// FIXME: Only consistent with ambipolar diffusion. Currently this is up to
	//        the user to enforce.
	double udivpe = 0.0;
	if ( G.dimensions==2 ) {
	    // Estimate electron pressure gradient as average of all vertices
	    double dpedx = 0.25 * (vtx[0]->dpedx + vtx[1]->dpedx + vtx[2]->dpedx + vtx[3]->dpedx);
	    double dpedy = 0.25 * (vtx[0]->dpedy + vtx[1]->dpedy + vtx[2]->dpedy + vtx[3]->dpedy);
	    // Approximation for work done on electrons: u dot div(pe)
            udivpe = fs->vel.x * dpedx + fs->vel.y * dpedy;
	} else {
	    // Estimate electron pressure gradient as average of all vertices
	    double dpedx = 0.125 * (vtx[0]->dpedx + vtx[1]->dpedx + vtx[2]->dpedx + vtx[3]->dpedx
			          + vtx[4]->dpedx + vtx[5]->dpedx + vtx[6]->dpedx + vtx[7]->dpedx);
	    double dpedy = 0.125 * (vtx[0]->dpedy + vtx[1]->dpedy + vtx[2]->dpedy + vtx[3]->dpedy
			          + vtx[4]->dpedy + vtx[5]->dpedy + vtx[6]->dpedy + vtx[7]->dpedy);
	    double dpedz = 0.125 * (vtx[0]->dpedz + vtx[1]->dpedz + vtx[2]->dpedz + vtx[3]->dpedz
			          + vtx[4]->dpedz + vtx[5]->dpedz + vtx[6]->dpedz + vtx[7]->dpedz);
	    // Approximation for work done on electrons: u dot div(pe)
	    udivpe = fs->vel.x * dpedx + fs->vel.y * dpedy + fs->vel.z * dpedz;
	}
	// FIXME: Assuming the free electron energy is included in the last mode
        Q->energies.back() += udivpe * G.diffusion_factor;
    } // end if ( G.electric_field_work )

    return SUCCESS;
} // end FV_Cell::add_viscous_source_vector()


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
} // end FV_Cell::calculate_wall_Reynolds_number()


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
double FV_Cell::rad_scaling_ratio(void) const
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
// This function must match the copy-to/from-buffer methods above.
// The buffers are typically used for communication between MPI processes.
{
    global_data &G = *get_global_data_ptr();
    int number = 0;
    Gas_model *gmodel = get_gas_model_ptr();
    if (type_of_copy == COPY_ALL_CELL_DATA ||
	type_of_copy == COPY_FLOW_STATE) {
        number += gmodel->number_of_values_in_gas_data_copy();
	number += 8 + 4; // FlowState + cell data
	if ( G.MHD ) number += 3;
    }
    if (type_of_copy == COPY_ALL_CELL_DATA ||
        type_of_copy == COPY_CELL_LENGTHS) {
        number += 6;
    }
    // [todo] PJ Do we really need all this interface data carried around?
    // I think that Andrew P. used it for a while and it is left over from the
    // early moving-grid work.
    if (type_of_copy == COPY_ALL_CELL_DATA ||
	type_of_copy == COPY_INTERFACE_DATA) {
	number += N_INTERFACE * gmodel->number_of_values_in_gas_data_copy();
	number += N_INTERFACE * 8; // FlowState data for each interface
	if ( G.MHD ) number += N_INTERFACE * 3;
	number += N_INTERFACE * 7; // Velocity, position and length for each interface
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
    global_data &G = *get_global_data_ptr();
    ostringstream ost;
    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t nmodes = gmodel->get_number_of_modes();
    ost << "\"pos.x\" \"pos.y\" \"pos.z\" \"volume\"";
    ost << " \"rho\" \"vel.x\" \"vel.y\" \"vel.z\" ";
    if ( G.MHD ) {
	ost << " \"B.x\" \"B.y\" \"B.z\" ";
    }
    ost << " \"p\" \"a\" \"mu\"";
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	ost << " \"k[" << imode << "]\"";
    }
    ost << " \"mu_t\" \"k_t\" \"S\"";
    if ( G.radiation ) {
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

