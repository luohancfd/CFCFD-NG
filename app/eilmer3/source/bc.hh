/** \file bc.hh
 * \ingroup eilmer3
 * \brief Prototypes for boundary-condition functions.
 */


#ifndef BC_HH
#define BC_HH

#include <string>
#include <vector>

#include "../../../lib/util/source/config_parser.hh"
#include "../../../lib/geometry2/source/geom.hh"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"
#include "../../../lib/nm/source/zero_system.hh"
#include "../../../lib/nm/source/zero_finders.hh"
#include "c-flow-condition.hh"
#include "cell.hh"
#include "flux_calc.hh"
#include "block.hh"

//-----------------------------------------------------------------

enum bc_t {
    ADJACENT, 
    SUP_IN, 
    EXTRAPOLATE_OUT, 
    SLIP_WALL, 
    ADIABATIC, 
    FIXED_T,
    SUBSONIC_IN, 
    SUBSONIC_OUT, 
    TRANSIENT_UNI, 
    TRANSIENT_PROF, 
    STATIC_PROF,
    FIXED_P_OUT,
    TRANSIENT_T_WALL,
    SEB,
    USER_DEFINED,
    ADJACENT_PLUS_UDF,
    ABLATING,
    SLIDING_T,
    FSTC,
    SHOCK_FITTING_IN,
    NON_CATALYTIC,
    EQUIL_CATALYTIC,
    SUPER_CATALYTIC,
    PARTIALLY_CATALYTIC
};
std::string get_bc_name(bc_t bc);

//-----------------------------------------------------------------

/** \brief Reflects normal velocity with respect to the cell interface.
 *
 * The process is to rotate the velocity vector into the local frame of
 * the interface, negate the normal (local x-component) velocity and
 * rotate back to the global frame.
 *
 * \param cell : pointer to a FV_cell object
 * \param IFace : pointer to a FV_Interface (contains direction vectors)
 *
 * The 3D frame-changing functions are in lib/geometry2/source/geom.*
 */
inline int reflect_normal_velocity(FV_Cell *cell, FV_Interface *IFace)
{
    cell->fs->vel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
    cell->fs->vel.x = -(cell->fs->vel.x);
    cell->fs->vel.transform_to_global(IFace->n, IFace->t1, IFace->t2);
    return SUCCESS;
}

inline int reflect_normal_magnetic_field(FV_Cell *cell, FV_Interface *IFace)
{
    cell->fs->B.transform_to_local(IFace->n, IFace->t1, IFace->t2);
    cell->fs->B.x = -(cell->fs->B.x);
    cell->fs->B.transform_to_global(IFace->n, IFace->t1, IFace->t2);
    return SUCCESS;
}

//-----------------------------------------------------------------
// Class-based boundary conditions in which all information about 
// each boundary condition is contained within a single class 
// rather than scattered across separate functions.
// The trade-off is that we have to have boundary selection within 
// the new functions.

int apply_convective_bc( Block &bd, double t, size_t dimensions );
int apply_viscous_bc( Block &bd, double t, size_t dimensions );

class CatalyticWallBC; // forward declaration needed below

class BoundaryCondition {
public:
    // Configuration data that are to be supplied for all BCs.
    Block *bdp;         // reference to the relevant block
    int which_boundary; // identity of the relevant boundary
    bc_t type_code;     // unique code used in switch statements and maps

    // Configuration data that can usually take default values.
    // Even if these data are specific to a only a few BCs,
    // they need to be here so that they are accessible through
    // pointers to objects of this base-class.
    bool is_wall_flag;
    bool ghost_cell_data_available;
    int xforce_flag;
    bool sets_conv_flux_flag;
    bool sets_visc_flux_flag;
    int neighbour_block;
    int neighbour_face;
    int neighbour_orientation;

    // *** FIX-ME ** catalytic configuration that needs to be worked on.
    bc_t wc_bc;
    CatalyticWallBC * cw;

    // Working space for heat-flux calculations.
    std::vector<double> q_cond;		// 1D vectors representing heat flux fields
    std::vector<double> q_diff;		// at surfaces (which are 2D for 3D grids)
    std::vector<double> q_rad;
    size_t imin, imax, jmin, jmax, kmin, kmax;	// boundary cell indice limits

    // Radiative emissivity (used for radiation transport)
    double emissivity;


public:
    BoundaryCondition(Block *bdp, int which_boundary, bc_t type_code);
    BoundaryCondition(); // Shouldn't have one without referring to a Block.
    BoundaryCondition(const BoundaryCondition &bc);
    BoundaryCondition & operator=(const BoundaryCondition &bc);
    virtual ~BoundaryCondition();
    virtual void print_info(std::string lead_in);

    virtual int apply_convective(double t); // reflect normal velocity
    virtual int apply_viscous(double t);  // does nothing

    bool is_wall()
    { return is_wall_flag; }
    // a true value indicates that this boundary is like a (solid) wall
    // some parts of the code (e.g. the turbulence models) need this information
    bool sets_conv_flux()
    { return sets_conv_flux_flag; } 
    // a true value indicates that the boundary condition is responsible for
    // setting the convective flux directly. When this is true, the internally
    // calculated flux is skipped and the b.c. flux is applied.
    bool sets_visc_flux()
    { return sets_visc_flux_flag; }
    // a true value indicates that the boundary condition is responsible for
    // setting the viscous flux directly. When this is true, the internally
    // calculated flux (based on spatial derivatives and interface propertiers)
    // is written over by the b.c. flux.
    double get_emissivity()
    { return emissivity; }

    int write_vertex_velocities(std::string filename, double sim_time,
				size_t dimensions, size_t gtl=0);
    // Heat-flux functions
    int compute_surface_heat_flux(void);
    int compute_cell_interface_surface_heat_flux(FV_Interface * IFace, FV_Cell * cell_one, 
    						 size_t index, size_t gtl=0);
    int write_surface_heat_flux(std::string filename, double sim_time);
    int write_fstc_heat_flux(string filename, double sim_time);
    double read_surface_heat_flux(std::string filename, size_t dimensions, int zip_files=1);
    int get_heat_flux_index( size_t i, size_t j, size_t k )
    { return (jmax-jmin+1)*(imax-imin+1)*(k-kmin) + (imax-imin+1)*(j-jmin) + (i-imin); }
};

BoundaryCondition *create_BC(Block *bdp, int which_boundary, bc_t type_of_BC,
			     ConfigParser &dict, std::string section);

int check_connectivity(void);

int scan_string_for_surface_heat_flux(double &q_cond, double &q_diff, double &q_rad, char *bufptr);

#endif
