/** \file bc.hh
 * \ingroup eilmer3
 * \brief Prototypes for boundary-condition functions.
 */


#ifndef BC_HH
#define BC_HH

#include <string>
#include <vector>

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
#include "bc_defs.hh"

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

//-----------------------------------------------------------------
// Class-based boundary conditions in which all information about 
// each boundary condition is contained within a single class 
// rather than scattered across separate functions.
// The trade-off is that we have to have boundary selection within 
// the new functions.

int apply_inviscid_bc( Block &bdp, double t, int dimensions );
int apply_viscous_bc( Block &bdp, double t, int dimensions );

class CatalyticWallBC; // forward declaration needed below

class BoundaryCondition {
public:
    Block &bdp; // reference to the relevant block
    int which_boundary; // identity of the relevant boundary
    int type_code; // value matches one of the integer type codes in bc_defs.hh
    string name_of_BC;
    bool is_wall_flag;
    bool use_udf_flux_flag;
    int neighbour_block;
    int neighbour_face;
    int neighbour_orientation;
    int wc_bc;
    int sponge_flag;
    int xforce_flag;
    std::vector<double> q_conv;		// 1D vectors representing heat flux fields
    std::vector<double> q_diff;		// at surfaces (which are 2D for 3D grids)
    std::vector<double> q_rad;
    int imin, imax, jmin, jmax, kmin, kmax;	// boundary cell indice limits
    CatalyticWallBC * cw;

public:
    BoundaryCondition( Block &bdp, int which_boundary, int type_code,
		       std::string name_of_BC="Unspecified", 
		       bool is_wall=false, bool use_udf_flux=false,
		       int neighbour_block=-1, int neighbour_face=-1,
		       int neighbour_orientation=0,
		       int wc_bc=NON_CATALYTIC, int sponge_flag=0, 
		       int xforce_flag=0 );
    BoundaryCondition( const BoundaryCondition &bc );
    virtual ~BoundaryCondition();
    virtual int apply_inviscid( double t ); // reflect normal velocity
    virtual int apply_viscous( double t );  // does nothing
    bool is_wall(); 
    // a true value indicates that this boundary is like a (solid) wall
    // some parts of the code (e.g. the turbulence models) need this information
    bool use_udf_flux(); 
    // a true value indicates that we want the UDF flux 
    // instead of the internally calculated flux
    void print_info( std::string lead_in );
    // Heat-flux functions
    int compute_surface_heat_flux( void );
    int compute_cell_interface_surface_heat_flux( FV_Interface * IFace, 
    						  FV_Cell * cell_one, 
    						  int index );
    int write_surface_heat_flux( std::string filename, double sim_time );
    int write_fstc_heat_flux( string filename, double sim_time );
    double read_surface_heat_flux( std::string filename, int dimensions, int zip_files=1 );
    int get_heat_flux_index( int i, int j, int k )
    { return (jmax-jmin+1)*(imax-imin+1)*(k-kmin) + (imax-imin+1)*(j-jmin) + (i-imin); }
};

BoundaryCondition *create_BC( Block &bdp, int which_boundary, int type_of_BC,
			      int inflow_condition_id, std::string filename, int n_profile,
			      double Twall, double Pout, int is_wall, int use_udf_flux,
			      int other_block, int other_face, int neighbour_orientation,
			      int sponge_flag, int xforce_flag,
			      std::vector<double> &mdot, double epsilon,
			      int wc_bc, std::string wcbc_fname, std::vector<double> f_wall,
			      double Twall_i, double Twall_f, double t_i, double t_f,
			      int assume_ideal);

int check_connectivity( void );

int scan_string_for_surface_heat_flux( double &q_conv, double &q_diff, double &q_rad, char *bufptr );

#endif
