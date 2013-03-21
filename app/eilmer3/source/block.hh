/// \file block.hh
/// \ingroup eilmer3
/// \brief Header file for the block data structure within Elmer3.
///
/// \author PJ
/// \version 05-Aug-04 extracted from mb_cns.h
///\version 02-Mar-08 Elmer3 port

#ifndef BLOCK_HH
#define BLOCK_HH

#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/string_util.hh"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "c-flow-condition.hh"
#include "flux_calc.hh"
#include "cell.hh"

class BoundaryCondition; // We need this forward declaration below.
struct global_data; // ...and this

// These definitions to simplify the definition of the apply functions,
// a little further below.
// As per the notes on http://www.parashift.com/c++-faq-lite/
typedef int (FV_Cell::*FV_Cell_MemberFunction_void)(void);
typedef int (FV_Cell::*FV_Cell_MemberFunction_double)(double);
typedef int (FV_Cell::*FV_Cell_MemberFunction_double_double)(double,double);
typedef int (FV_Cell::*FV_Cell_MemberFunction_int_double)(int,double);
typedef int (FV_Cell::*FV_Cell_MemberFunction_int)(int);
typedef int (FV_Cell::*FV_Cell_MemberFunction_int_int)(int,int);
#define CALL_MEMBER_FN(object,ptrToMember) ((object).*(ptrToMember))

/// \brief Number of ghost cells surrounding the active cells.
///
/// This sets the size of the ghost-cell buffer around a block.
const size_t NGHOST = 2;

/// \brief Parameters for cell checking...
///
/// When decoding the array of conserved quantities, 
/// the temperature or the density may try to go negative.  
/// If it does and adjust_invalid_cell_data == true, the cell data
/// is adjusted to make it reasonable.
/// If prefer_copy_from_left == true when adjusting data,
/// the new values are obtained from the cell to the left, else
/// the values are obtained by averaging the properties in nearby cells.
/// These parameters affect the behaviour of function count_invalid_cells().
const bool adjust_invalid_cell_data = true;
const bool prefer_copy_from_left = false;

const bool check_array_bounds = true;

/*----------------------------------------------------------------*/

/** \brief Single-block data structure.
 *
 * This data structure should contain everything needed for
 * a single-block solution -- both geometry and flow data
 */
class Block {
public:
    size_t id; // block identifier: assumed to be the same as the block number.
    int active; // =1: active block; =0: inactive block

    double omegaz; // Angular velocity (in rad/s) of the rotating frame.
                   // There is only one component, about the z-axis.

    double dt_allow;            // Allowable time step
    double cfl_min, cfl_max;    // estimates of CFL number

    double mass_residual, energy_residual; // monitor these for steady state
    Vector3 mass_residual_loc, energy_residual_loc; // locations of worst case

    size_t hncell;              // number of sample cells
    std::vector<size_t> hicell, hjcell, hkcell; // location of sample cells

    // Total number of cells in each direction for this block.
    // these will be used in the array allocation routines.
    size_t nidim, njdim, nkdim;

    // Number of active cells in the i,j index-directions
    size_t nni, nnj, nnk;

    // These index limits are set to allow convenient access 
    // to the arrays without having to worry how many buffer 
    // cells are present.
    // imin <= i <= imax,  jmin <= j <= jmax
    // Typically imin = jmin = 2.
    size_t imin, imax;
    size_t jmin, jmax;
    size_t kmin, kmax;

    // Flag to indicate if the Baldwin-Lomax turbulence model 
    // is active for the current block. 1==active; 0==inactive;
    int baldwin_lomax_iturb;

    // boundary-condition object pointers.
    std::vector<BoundaryCondition *> bcp;

    // Positions of interfaces marked as shocks
    // std::vector<Vector3 *> shock_iface_pos;

private:
    // Most of the data is stored in the following arrays.
    // ctr = cell center values
    // ifi = east-facing face properties and fluxes (unit normal in the i-index direction)
    // ifj = north-facing face properties and fluxes (normal in the j-index direction)
    // ifk = top-facing 
    // vtx = cell vertex values (used for the viscous terms, mostly)
    // sifi, sifj and sifk are secondary-cell faces (whose corner nodes are the
    //                     the primary-cell centres.
    std::vector<FV_Cell *> ctr_;
    std::vector<FV_Interface *> ifi_;
    std::vector<FV_Interface *> ifj_;
    std::vector<FV_Interface *> ifk_;
    std::vector<FV_Vertex *> vtx_;
    std::vector<FV_Interface *> sifi_;
    std::vector<FV_Interface *> sifj_;
    std::vector<FV_Interface *> sifk_;

public:
    Block();
    Block(const Block &b);
    ~Block();
    Block & operator=(const Block &b);

    int array_alloc(size_t dimensions);
    int array_cleanup(size_t dimensions);

    FV_Cell *get_cell(size_t i, size_t j, size_t k=0) 
    {
	if ( check_array_bounds && (k >= nkdim || j >= njdim || i >= nidim) ) {
	    throw std::runtime_error("Block::get_cell: index out of bounds: i="+
				     tostring(i)+" j="+tostring(j)+" k="+tostring(k));
	}
	return ctr_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_ifi(size_t i, size_t j, size_t k=0)
    {
	if ( check_array_bounds && (k >= nkdim || j >= njdim || i >= nidim) ) {
	    throw std::runtime_error("Block::get_ifi: index out of bounds: i="+
				     tostring(i)+" j="+tostring(j)+" k="+tostring(k));
	}
	return ifi_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_ifj(size_t i, size_t j, size_t k=0)
    {
	if ( check_array_bounds && (k >= nkdim || j >= njdim || i >= nidim) ) {
	    throw std::runtime_error("Block::get_ifj: index out of bounds: i="+
				     tostring(i)+" j="+tostring(j)+" k="+tostring(k));
	}
	return ifj_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_ifk(size_t i, size_t j, size_t k=0)
    {
	if ( check_array_bounds && (k >= nkdim || j >= njdim || i >= nidim) ) {
	    throw std::runtime_error("Block::get_ifk: index out of bounds: i="+
				     tostring(i)+" j="+tostring(j)+" k="+tostring(k));
	}
	return ifk_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Vertex *get_vtx(size_t i, size_t j, size_t k=0)
    {
	if ( check_array_bounds && (k >= nkdim || j >= njdim || i >= nidim) ) {
	    throw std::runtime_error("Block::get_vtx: index out of bounds: i="+
				     tostring(i)+" j="+tostring(j)+" k="+tostring(k));
	}
	return vtx_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_sifi(size_t i, size_t j, size_t k=0)
    {
	if ( check_array_bounds && (k >= nkdim || j >= njdim || i >= nidim) ) {
	    throw std::runtime_error("Block::get_sifi: index out of bounds: i="+
				     tostring(i)+" j="+tostring(j)+" k="+tostring(k));
	}
	return sifi_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_sifj(size_t i, size_t j, size_t k=0)
    {
	if ( check_array_bounds && (k >= nkdim || j >= njdim || i >= nidim) ) {
	    throw std::runtime_error("Block::get_sifj: index out of bounds: i="+
				     tostring(i)+" j="+tostring(j)+" k="+tostring(k));
	}
	return sifj_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_sifk(size_t i, size_t j, size_t k=0)
    {
	if ( check_array_bounds && (k >= nkdim || j >= njdim || i >= nidim) ) {
	    throw std::runtime_error("Block::get_sifk: index out of bounds: i="+
				     tostring(i)+" j="+tostring(j)+" k="+tostring(k));
	}
	return sifk_[k*(njdim*nidim)+j*nidim+i]; 
    }

    int apply(FV_Cell_MemberFunction_void f, string failure_message_header);
    int apply(FV_Cell_MemberFunction_double f, double param1, string failure_message_header);
    int apply(FV_Cell_MemberFunction_double_double f, double param1, double param2, 
	      string failure_message_header);
    int apply(FV_Cell_MemberFunction_int_double f, int param1, double param2, 
	      string failure_message_header);
    int apply(FV_Cell_MemberFunction_int f, int param1, string failure_message_header);
    int apply(FV_Cell_MemberFunction_int_int f, int param1, int param2, string failure_message_header);

    int apply(int (*f)(FV_Cell *cellp), 
	      string failure_message_header);
    int apply(int (*f)(FV_Cell *cellp, double param1), 
	      double param1, string failure_message_header);
    int apply(int (*f)(FV_Cell *cellp, int param1), 
	      int param1, string failure_message_header);
    int apply(int (*f)(FV_Cell *cellp, int param1, int param2), 
	      int param1, int param2, string failure_message_header);
    int apply(int (*f)(FV_Cell *cellp, int param1, double param2), 
	      int param1, double param2, string failure_message_header);

    int bind_interfaces_to_cells( size_t dimensions );
    int set_base_qdot( global_data &gdp ); 
    int identify_reaction_zones( global_data &gdp );
    int identify_turbulent_zones( global_data &gdp );
    int clear_fluxes_of_conserved_quantities( size_t dimensions );
    int propagate_data_west_to_east( size_t dimensions );

    int compute_initial_primary_cell_geometric_data( size_t dimensions );
    int compute_primary_cell_geometric_data( size_t dimensions, size_t time_level );
    int compute_distance_to_nearest_wall_for_all_cells( size_t dimensions );
    int compute_secondary_cell_geometric_data( size_t dimensions );
    int calc_initial_volumes_2D( void );
    int calc_volumes_2D( size_t time_level );
    int secondary_areas_2D( void );
    int calc_initial_faces_2D( void );
    int calc_faces_2D( size_t time_level );
    int calc_initial_ghost_cell_geom_2D( void );
    int calc_ghost_cell_geom_2D( size_t time_level );
    int predict_vertex_positions( size_t dimensions, double dt );
    int correct_vertex_positions( size_t dimensions, double dt );
    int rk3_vertex_positions( size_t dimensions, double dt );
    int set_geometry_velocities( size_t dimensions, size_t time_level );
    int set_vertex_velocities2D( size_t time_level );
    int set_gcl_test_vertex_velocities2D( size_t time_level );
    int set_gcl_test_vertex_velocities3D( size_t time_level );
    int set_gcl_test_random_vertex_velocities2D( size_t time_level );
    int set_gcl_interface_properties( size_t dimensions, size_t time_level, double dt );
    int set_gcl_interface_properties2D( size_t time_level, double dt );
    int set_gcl_interface_properties3D( size_t time_level, double dt );
    int set_interface_velocities2D( size_t time_level );
    int set_vertex_velocities3D( size_t time_level );
    int set_interface_velocities3D( size_t time_level );
    int calc_boundary_vertex_velocity(FV_Interface &IFace1, FV_Interface &IFace2,     
				      FV_Vertex &vtx, Vector3 trv, size_t time_level);
    int calc_boundary_vertex_velocity(FV_Interface &IFace1, FV_Interface &IFace2,     
				      FV_Interface &IFace3, FV_Interface &IFace4,
				      FV_Vertex &vtx, Vector3 trv, size_t time_level);
    int velocity_weighting_factor(FV_Interface &IFace, Vector3 vp, double &w, Vector3 &ws);
    int diffuse_vertex_velocities(double mu, int npass, size_t dimensions, size_t time_level);
    int anti_diffuse_vertex_velocities(double mu, int npass, size_t dimensions, size_t time_level);
    int compute_boundary_flux(FV_Interface *IFaceL, FV_Interface *IFaceR, double omegaz);
    
    void compute_x_forces( char *text_string, int ibndy, size_t dimensions );
    int print_forces( FILE *fp, double t, size_t dimensions );

    int read_grid(std::string filename, size_t dimensions, int zip_file=1);
    int read_solution(std::string filename, double *sim_time, 
		      size_t dimensions, int zip_file=1);
    int write_solution(std::string filename, double sim_time, 
		       size_t dimensions, int zip_file=1 );
    int write_block(std::string filename, double sim_time, 
		       size_t dimensions, int zip_file=1 );
    int read_BGK(std::string filename, double *sim_time, 
		      size_t dimensions, int zip_file=1);
    int initialise_BGK_equilibrium( void );
    int write_BGK(std::string filename, double sim_time, 
		       size_t dimensions, int zip_file=1 );
    int write_history( std::string filename, double sim_time, int write_header=0 );

    int count_invalid_cells( size_t dimensions );
    int init_residuals( size_t dimensions );
    int compute_residuals( size_t dimensions );

    int determine_time_step_size( double cfl_target, size_t dimensions );
    int detect_shock_points( size_t dimensions );
    int apply_spatial_filter_diffusion( double alpha, size_t npass, size_t dimensions );
    int apply_spatial_filter_anti_diffusion( double alpha, size_t npass, size_t dimensions );
    double calc_anti_diffusive_flux(double m2, double m1, double p1, double p2, double mu);

    // The implementation for the following function is in invs.cxx
    int inviscid_flux( size_t dimensions );
   
    
};  /* end of the (single-)block data structure definition */

int find_nearest_cell( double x, double y, double z, 
		       size_t *jb_near, size_t *i_near, size_t *j_near, size_t *k_near );
int locate_cell(double x, double y, double z,
		size_t *jb_found, size_t *i_found, size_t *j_found, size_t *k_found);


/** Indexing of the data in 2D.
 *
 * \verbatim
 * The following figure shows cell [i,j] and its associated
 * vertices and faces. 
 * (New arrangement, planned August 2006, implemented Nov 2006)
 *
 *
 *
 *     Vertex 3         North face           Vertex 2 
 *   vtx[i,j+1]         ifj[i,j+1]           vtx[i+1,j+1]
 *             +--------------x--------------+
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *   West      |         cell center         |  East 
 *   face      |          ctr[i,j]           |  face
 *   ifi[i,j]  x              o              x  ifi[i+1,j]
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             +--------------x--------------+
 *     Vertex 0           South face         Vertex 1
 *     vtx[i,j]           ifj[i,j]           vtx[i+1,j]
 *
 *
 * Thus...
 * ----
 * Active cells are indexed as ctr[i][i], where
 * imin <= i <= imax, jmin <= j <= jmax.
 *
 * Active east-facing interfaces are indexed as ifi[i][j], where
 * imin <= i <= imax+1, jmin <= j <= jmax.
 *
 * Active north-facing interfaces are indexed as ifj[i][j], where
 * imin <= i <= imax, jmin <= j <= jmax+1.
 *
 * Active vertices are indexed as vtx[i][j], where
 * imin <= i <= imax+1, jmin <= j <= jmax+1.
 *
 * Space for ghost cells is available outside these ranges.
 * \endverbatim
 */


/// Indexing for the 3D data -- see page 8 in 3D CFD workbook.
///
/// These macros should make indexing over the vertices
/// and over the interfaces more readable.
/// We shall use VTK notation for a hexahedral cell and
/// define the macros with respect to the [i][j][k] cell.
/// 0 i   j   k
/// 1 i+1 j   k
/// 2 i+1 j+1 k
/// 3 i   j+1 k
/// 4 i   j   k+1
/// 5 i+1 j   k+1
/// 6 i+1 j+1 k+1
/// 7 i   j+1 k+1
 
/*-------------------------------------------------------------------*/

#endif
