/** \file block.hh
 * \ingroup eilmer3
 * \brief Header file for the block data structure within Elmer3.
 * 
 * \author PJ
 * \version 05-Aug-04 extracted from mb_cns.h
 * \version 02-Mar-08 Elmer3 port
 */


#ifndef BLOCK_HH
#define BLOCK_HH

#include <string>
#include "../../../lib/util/source/useful.h"
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
typedef int (FV_Cell::*FV_Cell_MemberFunction_int)(int);
typedef int (FV_Cell::*FV_Cell_MemberFunction_int_int)(int,int);
#define CALL_MEMBER_FN(object,ptrToMember) ((object).*(ptrToMember))

/// \brief Number of ghost cells surrounding the active cells.
///
/// This sets the size of the ghost-cell buffer around a block.
#define NGHOST 2

#define MAX_HNCELL 200

/// \brief Parameters for cell checking...
///
/// When decoding the array of conserved quantities, 
/// the temperature or the density may try to go negative.  
/// If it does and ADJUST_INVALID_CELL_DATA == 1, the cell data
/// is adjusted to make it reasonable.
/// If PREFER_COPY_FROM_LEFT == 1 when adjusting data,
/// the new values are obtained from the cell to the left, else
/// the values are obtained by averaging the properties in nearby cells.
/// These parameters affect the behaviour of function count_invalid_cells().
#define ADJUST_INVALID_CELL_DATA 1
#define PREFER_COPY_FROM_LEFT    0

#define CHECK_ARRAY_BOUNDS 1

/*----------------------------------------------------------------*/

/** \brief Single-block data structure.
 *
 * This data structure should contain everything needed for
 * a single-block solution -- both geometry and flow data
 */
class Block {
public:
    int id; // block identifier: assumed to be the same as the block number.
    int active; // =1: active block; =0: inactive block

    double omegaz; // Angular velocity (in rad/s) of the rotating frame.
                   // There is only one component, about the z-axis.

    double dt_allow;            // Allowable time step
    double cfl_min, cfl_max;    // estimates of CFL number

    double mass_residual, energy_residual; // monitor these for steady state
    Vector3 mass_residual_loc, energy_residual_loc; // locations of worst case

    int hncell;                 // number of sample cells
    int hicell[MAX_HNCELL], hjcell[MAX_HNCELL], hkcell[MAX_HNCELL]; // location of sample cell

    // Total number of cells in each direction for this block.
    // these will be used in the array allocation routines.
    int nidim, njdim, nkdim;

    // Number of active cells in the i,j index-directions
    int nni, nnj, nnk;

    // These index limits are set to allow convenient access 
    // to the arrays without having to worry how many buffer 
    // cells are present.
    // imin <= i <= imax,  jmin <= j <= jmax
    // Typically imin = jmin = 2.
    int imin, imax;
    int jmin, jmax;
    int kmin, kmax;

    // Flag to indicate if the Baldwin-Lomax turbulence model 
    // is active for the current block. 1==active; 0==inactive;
    int baldwin_lomax_iturb;

    // boundary-condition object pointers.
    BoundaryCondition *bcp[6];

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
    ~Block();
    int array_alloc(int dimensions);
    int array_cleanup(int dimensions);

    FV_Cell *get_cell(int i, int j, int k=0) 
    {
#       if CHECK_ARRAY_BOUNDS == 1
	if ( (k < 0) || (k >= nkdim) || 
	     (j < 0) || (j >= njdim) || 
	     (i < 0) || (i >= nidim) ) {
		 printf("Block::get_cell: index out of bounds: i=%d, j=%d, k=%d\n", i, j, k);
	         printf("                 nidim=%d, njdim=%d, nkdim=%d\n", nidim, njdim, nkdim);
		 exit( INDEX_OUT_OF_RANGE );
	}
#       endif 
	return ctr_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_ifi(int i, int j, int k=0)
    {
#       if CHECK_ARRAY_BOUNDS == 1
	if ( (k < 0) || (k >= nkdim) || 
	     (j < 0) || (j >= njdim) || 
	     (i < 0) || (i >= nidim) ) {
		 printf("Block::get_ifi: index out of bounds: i=%d, j=%d, k=%d\n", i, j, k);
	         printf("                nidim=%d, njdim=%d, nkdim=%d\n", nidim, njdim, nkdim);
		 exit( INDEX_OUT_OF_RANGE );
	}
#       endif 
	return ifi_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_ifj(int i, int j, int k=0)
    {
#       if CHECK_ARRAY_BOUNDS == 1
	if ( (k < 0) || (k >= nkdim) || 
	     (j < 0) || (j >= njdim) || 
	     (i < 0) || (i >= nidim) ) {
		 printf("Block::get_ifj: index out of bounds: i=%d, j=%d, k=%d\n", i, j, k);
	         printf("                nidim=%d, njdim=%d, nkdim=%d\n", nidim, njdim, nkdim);
		 exit( INDEX_OUT_OF_RANGE );
	}
#       endif 
	return ifj_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_ifk(int i, int j, int k=0)
    {
#       if CHECK_ARRAY_BOUNDS == 1
	if ( (k < 0) || (k >= nkdim) || 
	     (j < 0) || (j >= njdim) || 
	     (i < 0) || (i >= nidim) ) {
		 printf("Block::get_ifk: index out of bounds: i=%d, j=%d, k=%d\n", i, j, k);
	         printf("                nidim=%d, njdim=%d, nkdim=%d\n", nidim, njdim, nkdim);
		 exit( INDEX_OUT_OF_RANGE );
	}
#       endif 
	return ifk_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Vertex *get_vtx(int i, int j, int k=0)
    {
#       if CHECK_ARRAY_BOUNDS == 1
	if ( (k < 0) || (k >= nkdim) || 
	     (j < 0) || (j >= njdim) || 
	     (i < 0) || (i >= nidim) ) {
		 printf("Block::get_vtx: index out of bounds: i=%d, j=%d, k=%d\n", i, j, k);
	         printf("                nidim=%d, njdim=%d, nkdim=%d\n", nidim, njdim, nkdim);
		 exit( INDEX_OUT_OF_RANGE );
	}
#       endif 
	return vtx_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_sifi(int i, int j, int k=0)
    {
#       if CHECK_ARRAY_BOUNDS == 1
	if ( (k < 0) || (k >= nkdim) || 
	     (j < 0) || (j >= njdim) || 
	     (i < 0) || (i >= nidim) ) {
		 printf("Block::get_sifi: index out of bounds: i=%d, j=%d, k=%d\n", i, j, k);
	         printf("                nidim=%d, njdim=%d, nkdim=%d\n", nidim, njdim, nkdim);
		 exit( INDEX_OUT_OF_RANGE );
	}
#       endif 
	return sifi_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_sifj(int i, int j, int k=0)
    {
#       if CHECK_ARRAY_BOUNDS == 1
	if ( (k < 0) || (k >= nkdim) || 
	     (j < 0) || (j >= njdim) || 
	     (i < 0) || (i >= nidim) ) {
		 printf("Block::get_sifj: index out of bounds: i=%d, j=%d, k=%d\n", i, j, k);
	         printf("                nidim=%d, njdim=%d, nkdim=%d\n", nidim, njdim, nkdim);
		 exit( INDEX_OUT_OF_RANGE );
	}
#       endif 
	return sifj_[k*(njdim*nidim)+j*nidim+i]; 
    }
    FV_Interface *get_sifk(int i, int j, int k=0)
    {
#       if CHECK_ARRAY_BOUNDS == 1
	if ( (k < 0) || (k >= nkdim) || 
	     (j < 0) || (j >= njdim) || 
	     (i < 0) || (i >= nidim) ) {
		 printf("Block::get_sifk: index out of bounds: i=%d, j=%d, k=%d\n", i, j, k);
	         printf("                nidim=%d, njdim=%d, nkdim=%d\n", nidim, njdim, nkdim);
		 exit( INDEX_OUT_OF_RANGE );
	}
#       endif 
	return sifk_[k*(njdim*nidim)+j*nidim+i]; 
    }

    int apply(FV_Cell_MemberFunction_void f, string failure_message_header);
    int apply(FV_Cell_MemberFunction_double f, double param1, string failure_message_header);
    int apply(FV_Cell_MemberFunction_double_double f, double param1, double param2, 
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

    int bind_interfaces_to_cells( int dimensions );
    int set_base_qdot( global_data &gdp ); 
    int identify_reaction_zones( global_data &gdp );
    int identify_turbulent_zones( global_data &gdp );
    int clear_fluxes_of_conserved_quantities( int dimensions );
    int propagate_data_west_to_east( int dimensions );

    int compute_primary_cell_geometric_data( int dimensions );
    int compute_distance_to_nearest_wall_for_all_cells( int dimensions );
    int compute_secondary_cell_geometric_data( int dimensions );
    int calc_volumes_2D( void );
    int secondary_areas_2D( void );
    int calc_faces_2D( void );

    void compute_x_forces( char *text_string, int ibndy, int dimensions );
    int print_forces( FILE *fp, double t, int dimensions );

    int read_grid( std::string filename, int dimensions, int zip_file=1 );
    double read_solution( std::string filename, int dimensions, int zip_file=1 );
    int write_solution( std::string filename, double sim_time, int dimensions, int zip_file=1 );
    int write_history( std::string filename, double sim_time, int write_header=0 );

    int count_invalid_cells( int dimensions );
    int init_residuals( int dimensions );
    int compute_residuals( int dimensions );

    int determine_time_step_size( double cfl_target, int dimensions );
    int detect_shock_points( int dimensions );
    int apply_spatial_filter( double alpha, int npass, int dimensions );

    // The implementation for the following function is in invs.cxx
    int inviscid_flux( int dimensions );
    
};  /* end of the (single-)block data structure definition */

int find_nearest_cell( double x, double y, double z, 
		       int *jb_near, int *i_near, int *j_near, int *k_near );
int locate_cell(double x, double y, double z,
		int *jb_found, int *i_found, int *j_found, int *k_found);


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


/// Indexing Macros for the 3D data -- see page 8 in 3D CFD workbook.
///
/// These macros should make indexing over the vertices
/// and over the interfaces more readable.
/// We shall use VTK notation for a hexahedral cell and
/// define the macros with respect to the [i][j][k] cell.
#define IVTX0 (i)
#define JVTX0 (j)
#define KVTX0 (k)

#define IVTX1 (i+1)
#define JVTX1 (j)
#define KVTX1 (k)

#define IVTX2 (i+1)
#define JVTX2 (j+1)
#define KVTX2 (k)

#define IVTX3 (i)
#define JVTX3 (j+1)
#define KVTX3 (k)

#define IVTX4 (i)
#define JVTX4 (j)
#define KVTX4 (k+1)

#define IVTX5 (i+1)
#define JVTX5 (j)
#define KVTX5 (k+1)

#define IVTX6 (i+1)
#define JVTX6 (j+1)
#define KVTX6 (k+1)

#define IVTX7 (i)
#define JVTX7 (j+1)
#define KVTX7 (k+1)

/// Corners of the 2D secondary cell in 2D, located at vtx[i][j].
/// Note that these are indices of cell-centres wrt vtx[i][j].
#define  IADSH  (i)
#define  JADSH  (j-1)
#define  IBDSH  (i)
#define  JBDSH  (j)
#define  ICDSH  (i-1)
#define  JCDSH  (j)
#define  IDDSH  (i-1)
#define  JDDSH  (j-1)

 
/*-------------------------------------------------------------------*/

#endif
