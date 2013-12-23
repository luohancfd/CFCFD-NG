/// \file piston.cxx
/// \ingroup eilmer3
/// \brief Class for handling pistons.
/// \author Brendan O'Flaherty and Rowan J. Gollan
///
/// \version 16-Mar-2007 -- Began copying from mbcns2 piston
/// \version 31-Mar-2007 -- Both codes consistant, friction added along with two test cases. 
/// \version 12-Jun-2007 -- Piston faces now composed of multiple component faces.
/// \version 07-Jan-2008 -- Ported from two_phase code
///
/// ---------------------------------------------------------

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "../../../lib/util/source/useful.h"
#include "kernel.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "piston.hh"

constexpr int DB0 = 0; // function calls
constexpr int DB1 = 0; // index checking
constexpr int DB2 = 0; // initialising
constexpr int DB3 = 0; // derivatives
constexpr int DB4 = 0; // piston information
constexpr int DB5 = 0; // cell information

using namespace std;

// NW vtx[bd[].imin,bd[].jmax+1]
// NE vtx[bd[].imax+1,bd[].jmax+1]
// SE vtx[bd[].imax+1,bd[].jmin]
// SW vtx[bd[].imin,bd[].jmin]

// ----------------------------------------------------------------------------
// Piston face implementation

// default constructor
PistonFace::PistonFace() 
    : id_(0),
      io_(0), bo_(0), ii_(0), bi_(0),
      jmin_(0), jmax_(0),
      nni_(0), nnj_(0),
      nvi_(0), nvj_(0),
      nii_(0), nij_(0),
      njcells_(0),
      dFg_(0.0), dR_(0.0),
      n_(), vtx_C_(), vtx_L_(), vtx_R_(),
      type_(UNSPECIFIED) {}

// copy constructor
PistonFace::PistonFace(const PistonFace &p) 
    : id_(p.id_),
      io_(p.io_), bo_(p.bo_), ii_(p.ii_), bi_(p.bi_),
      jmin_(p.jmin_), jmax_(p.jmax_),
      nni_(p.nni_), nnj_(p.nnj_),
      nvi_(p.nvi_), nvj_(p.nvj_),
      nii_(p.nii_), nij_(p.nij_),
      njcells_(p.njcells_),
      dFg_(p.dFg_), dR_(p.dR_),
      n_(p.n_), vtx_C_(p.vtx_C_), vtx_L_(p.vtx_L_), vtx_R_(p.vtx_R_),
      vel_(p.vel_),
      type_(p.type_)
{
    if (DB0) cout << "PistonFace copy constructor called" << endl;

    Gas_model *gm = get_gas_model_ptr();

    ifi_.resize(p.ifi_.size());
    for (size_t i=0; i<ifi_.size(); ++i) {
	ifi_[i] = new FV_Interface(gm);
	ifi_[i]->copy_values_from(*(p.ifi_[i]), COPY_ALL_CELL_DATA);
    }
    ifj_.resize(p.ifj_.size());
    for (size_t i=0; i<ifj_.size(); ++i) {
	ifj_[i] = new FV_Interface(gm);
	ifj_[i]->copy_values_from(*(p.ifj_[i]), COPY_ALL_CELL_DATA);
    }
    vtx_.resize(p.vtx_.size());
    for (size_t i=0; i<vtx_.size(); ++i) {
	vtx_[i] = new FV_Vertex(gm); 
	vtx_[i]->copy_values_from(*(p.vtx_[i]));
    }
    flex_cells_.resize(p.flex_cells_.size());
    for (size_t i=0; i<flex_cells_.size(); ++i) {
	flex_cells_[i] = new FV_Cell(gm);
	flex_cells_[i]->copy_values_from(*(p.flex_cells_[i]), COPY_ALL_CELL_DATA, 0);
    }
}

PistonFace::~PistonFace() 
{
    for( size_t i = 0; i < flex_cells_.size(); ++i ) delete flex_cells_[i];
    for( size_t i = 0; i < ifi_.size(); ++i ) delete ifi_[i];
    for( size_t i = 0; i < ifj_.size(); ++i ) delete ifj_[i];
    for( size_t i = 0; i < vtx_.size(); ++i ) delete vtx_[i];
}

// assignment operator
PistonFace& PistonFace::operator=(const PistonFace &pf)
{ 
    if (this == &pf) // avoid aliasing
	return *this; 
    else {
	// We need to be careful with memory.
	// 1. First clean out any existing pointers and
	//    associated memory
	for( size_t i = 0; i < flex_cells_.size(); ++i ) delete flex_cells_[i];
	for( size_t i = 0; i < ifi_.size(); ++i ) delete ifi_[i];
	for( size_t i = 0; i < ifj_.size(); ++i ) delete ifj_[i];
	for( size_t i = 0; i < vtx_.size(); ++i ) delete vtx_[i];
	// 2. Do all the simple assigments to preserve
	//    semantic equivalence
	id_ = pf.id_;
	io_ = pf.io_; bo_ = pf.bo_; ii_ = pf.ii_; bi_ = pf.bi_;
	jmin_ = pf.jmin_; jmax_ = pf.jmax_;
	nni_ = pf.nni_; nnj_ = pf.nnj_;
	nvi_ = pf.nvi_; nvj_ = pf.nvj_;
	nii_ = pf.nii_; nij_ = pf.nij_;
	njcells_ = pf.njcells_;
	dFg_ = pf.dFg_; dR_ = pf.dR_;
	n_ = pf.n_;  vtx_C_ = pf.vtx_C_;  vtx_L_ = pf.vtx_L_; vtx_R_ = pf.vtx_R_;
	vel_ = pf.vel_;
	type_ = pf.type_;
	// 3. Do the "trickier" assignment of memory-allocated objects.
	Gas_model *gm = get_gas_model_ptr();
	ifi_.resize(pf.ifi_.size());
	for (size_t i=0; i<ifi_.size(); ++i) {
	    ifi_[i] = new FV_Interface(gm);
	    ifi_[i]->copy_values_from(*(pf.ifi_[i]), COPY_ALL_CELL_DATA);
	}
	ifj_.resize(pf.ifj_.size());
	for (size_t i=0; i<ifj_.size(); ++i) {
	    ifj_[i] = new FV_Interface(gm);
	    ifj_[i]->copy_values_from(*(pf.ifj_[i]), COPY_ALL_CELL_DATA);
	}
	vtx_.resize(pf.vtx_.size());
	for (size_t i=0; i<vtx_.size(); ++i) {
	    vtx_[i] = new FV_Vertex(gm); 
	    vtx_[i]->copy_values_from(*(pf.vtx_[i]));
	}
	flex_cells_.resize(pf.flex_cells_.size());
	for (size_t i=0; i<flex_cells_.size(); ++i) {
	    flex_cells_[i] = new FV_Cell(gm);
	    flex_cells_[i]->copy_values_from(*(pf.flex_cells_[i]), COPY_ALL_CELL_DATA, 0);
	}
    }
    return (*this);
}

// ----------------------------------------------------------------------------

void PistonFace:: move_face(Vector3 &pos)
{
    vtx_C_ = pos;
    vtx_L_.x = pos.x;
    vtx_R_.x = pos.x;
}

int PistonFace:: set_values(double n_x, Vector3 vtxL, Vector3 vtxR, Vector3* vel, PistonFaceType type)
{
    n_.x = n_x;
    vtx_C_ = (vtxL + vtxR)/2.0;
    vtx_L_ = vtxL;
    vtx_R_ = vtxR;
    vel_ = vel;
    type_ = type;

    return SUCCESS;
}

int PistonFace:: set_indices_for_flex_cells(Block &bd)
{
    global_data &G = *get_global_data_ptr();
    if (DB0) cout << "entering set_indices_for_flex_cells()" << endl;
    nni_ = 1;
    nvi_ = nni_ + 1;
    nii_ = nvi_;
    
    nnj_ = bd.njdim;
    nvj_ = nnj_+ 1 + 2*G.nghost;
    nij_ = nvj_;
    njcells_ = nnj_ + 2*G.nghost;
    
    jmin_ = bd.jmin;
    jmax_ = bd.jmax;

    if (DB1) {
	cout << "nni_ = " << nni_ << endl;
	cout << "nvi_ = " << nvi_ << endl;
	cout << "nii_ = " << nii_ << endl;

	cout << "nnj_ = " << nnj_ << endl;
	cout << "nvj_ = " << nvj_ << endl;
	cout << "nij_ = " << nij_ << endl;
	cout << "njcells_ = " << njcells_ << endl;
    }
    
    if (DB0) cout << "leaving set_indices_for_flex_cells()" << endl;
    return SUCCESS;
}

static const double min_allowable_flex_fraction = 0.5;
static const double WIDTH_TOL = 1.0e-12;

int PistonFace:: locate_shadow_indices(global_data &G, Block *(bd[]))
{
    if (DB0) cout << "entering locate_shadow_indices()" << endl;

    size_t i_near, b_index;
    double x_inside, x_outside;
    double covered_cell_width, allowable_cell_width;

    if (type_ == WEST_FACE) {
	// pos is the position of the piston surface
	if (locate_west_block_index(b_index, vtx_L_, G, bd) != SUCCESS) return FAILURE;
	linear_search_from_west(i_near, vtx_L_.x, *bd[b_index]);
	x_outside = bd[b_index]->get_ifi(i_near-1,bd[b_index]->jmin)->pos.x;
	x_inside = bd[b_index]->get_ifi(i_near,bd[b_index]->jmin)->pos.x;
	
	covered_cell_width = fabs(vtx_L_.x - x_outside);
	allowable_cell_width = min_allowable_flex_fraction*fabs(x_inside - x_outside);
	
	if ( (covered_cell_width > allowable_cell_width) ||
	     (fabs(covered_cell_width - allowable_cell_width) <= WIDTH_TOL) ) {
	    // Case A:
	    // Flex-cell is greater than half a cell.
	    // Given we know the cell interface, return index of cell to the west
	    covers_two_cells_ = false;
	    io_ = i_near - 1;
	}
	else {
	    // Case B:
	    // if this index results in a flex_cell that's less than min_allowable_flex_fraction*grid_cell,
	    // shuffle left
	    covers_two_cells_ = true;
	    io_ = i_near - 2;
	}
	bo_ = bi_ = b_index;
	ii_ = io_ + 1;

	// At this point we have assumed the piston is somewhere in the
	// block away from edges.  Now test to see if we are too close to the edge.
	if( covers_two_cells_ ) {
	    // Check that the outer cell is within this block.
	    if( io_ < bd[bo_]->imin ) {
		// Ww are the edge.
		// Question: does this block have a connection on the west edge.
		if( bd[bo_]->bcp[WEST]->neighbour_block == -1 ) {
		    // There is no connection.
		    // We cannot sensibly proceed.
		    cout << "During locate_cell_index() when searching for the\n";
		    cout << "WEST face, the projectile approached the end of block.\n";
		    cout << "It is sensible to terminate at this point.\n";
		    return END_OF_BLOCK;
		}
		else {
		    // We are connected to something.
		    bo_ = bd[bo_]->bcp[WEST]->neighbour_block;
		    io_ = bd[bo_]->imax;
		}
	    }
	}
    }
    else {
	// pos is the position of the piston surface
	if (locate_east_block_index(b_index, vtx_R_, G, bd) != SUCCESS) return FAILURE;
	linear_search_from_east(i_near, vtx_R_.x, *bd[b_index]);
	x_inside = bd[b_index]->get_ifi(i_near,bd[b_index]->jmin)->pos.x;
	x_outside = bd[b_index]->get_ifi(i_near+1,bd[b_index]->jmin)->pos.x;

	i_near += 0; // i_near is the interface, and now it points
	             // to the neighbouring cell (which has the same index)

	covered_cell_width = fabs(x_outside - vtx_R_.x);
	allowable_cell_width = min_allowable_flex_fraction*fabs(x_inside - x_outside);
	    
	if ( (covered_cell_width > allowable_cell_width) ||
	     (fabs(covered_cell_width - allowable_cell_width) <= WIDTH_TOL) ) {
	    // Case A:
	    // Flex-cell is greater than half a cell.
	    // Given we know the cell interface, return the index of cell to the east
	    covers_two_cells_ = false;
	    io_ = i_near + 0;
	}
	else {
	    // Case B:
	    // Flex-cell is less than half a cell.
	    // Given we know the cell interface, return the index of the cell to the east + 1
	    covers_two_cells_ = true;
	    io_ = i_near + 1;
	}
	bo_ = bi_ = b_index;
	ii_ = io_ - 1;

	// At this point we have assumed the piston is somewhere in the
	// block away from edges.  Now test to see if we are too close to the edge.
	if( covers_two_cells_ ) {
	    // Check that the outer cell is within this block.
	    if( io_ > bd[bo_]->imax ) {
		// Ww are the edge.
		// Question: does this block have a connection on the east edge.
		if( bd[bo_]->bcp[EAST]->neighbour_block == -1 ) {
		    // There is no connection.
		    // We cannot sensibly proceed.
		    cout << "During locate_cell_index() when searching for the\n";
		    cout << "EAST face, the projectile approached the end of block.\n";
		    cout << "It is sensible to terminate at this point.\n";
		    return END_OF_BLOCK;
		}
		else {
		    // We are connected to something.
		    bo_ = bd[bo_]->bcp[EAST]->neighbour_block;
		    io_ = bd[bo_]->imin;
		}
	    }
	}
    }

    if (DB0) cout << "leaving locate_shadow_indices()" << endl;
    
    return SUCCESS;
}

/// \brief Allocates memory for the vectors of pointers to structures
///
/// Given jmin, jmax are set, thie function allocates memory
/// for a column of flex_cells.  This includes:
///  - vertices
///  - i-interfaces
///  - j-interfaces 
///  - flex_cells (composed of aforementioned vertices and interfaces)
///
int PistonFace::allocate_memory()
{

    if (DB0) cout << "entering allocate_memory()" << endl;
    size_t i, j, k, indx;
    Gas_model *gm = get_gas_model_ptr();
   
    if (DB2) cout << "allocating vertices" << endl;
    // vertices
    vtx_.resize(nvi_*nvj_, 0);
    for (k = 0; k < nvi_; ++k) {
	for (j = jmin_; j <= jmax_+1; ++j) {
	    indx = nvi_*j+k;
	    // Purely memory allocation.
	    // These are fixed when grid data is pasted into flex_cells
	    vtx_[indx] = new FV_Vertex(gm); 
	    vtx_[indx]->id = indx;
	}
    }
    
    if (DB2) cout << "allocating i-interfaces" << endl;
    // i-interfaces
    ifi_.resize(nii_*nij_, 0);
    for (i = 0; i < nii_; ++i) {
	for (j = jmin_; j <= jmax_; ++j) {
	    indx = j*nii_ + i;
	    ifi_[indx] = new FV_Interface(gm);
	    ifi_[indx]->id = indx;
	}
    }
    
    if (DB2) cout << "allocating j-interfaces" << endl;
    // j-interfaces
    ifj_.resize(nii_*nij_, 0);
    for (i = 0; i < (nii_-1); ++i) {
	for (j = jmin_; j <= jmax_+1; ++j) {
	    indx = j*(nii_-1) + i;
	    ifj_[indx] = new FV_Interface(gm);
	    ifj_[indx]->id = indx;
	}
    }
    
    if (DB2) cout << "allocating cells" << endl;
    // cells
    flex_cells_.resize(nni_*njcells_, 0);
    for (i = 0; i < nni_; ++i) {
	for (j = jmin_; j <= jmax_; ++j) {
	    indx = j*nni_ + i;
	    flex_cells_[indx] = new FV_Cell(gm);
	}
    }
    
    bind_interfaces_to_flex_cells();

    if (DB0) cout << "leaving allocate_memory()" << endl;
    return SUCCESS;
}

int PistonFace::bind_interfaces_to_flex_cells()
{
/**
\brief
Here we point the addresses for vtx and iface in each of the flex_cells 
to the appropriate structs stored in the PistonFace. This should be
done once after memory allocation for each component PistonFace.
**/
    if (DB0) cout << "entering bind_interfaces_to_flex_cells()" << endl;
    size_t i, j, indx;

    FV_Cell *fcp;
    for (i = 0; i < nni_; ++i) {
	for (j = jmin_; j <= jmax_; ++j) {
	    indx = j*nni_ + i;
	    fcp = flex_cells_[indx]; 

	    fcp->vtx[0] = get_vertex(i,j);
	    fcp->vtx[1] = get_vertex(i+1,j);
	    fcp->vtx[2] = get_vertex(i+1,j+1);
	    fcp->vtx[3] = get_vertex(i,j+1);
	    fcp->vtx[4] = NULL;
	    fcp->vtx[5] = NULL;
	    fcp->vtx[6] = NULL;
	    fcp->vtx[7] = NULL;
	}
    }

    for (i = 0; i < nni_; ++i) {
	for (j = jmin_; j <= jmax_; ++j) {
	    indx = j*nni_ + i;
	    fcp = flex_cells_[indx]; 

	    fcp->iface[NORTH] = get_ifj(j+1);
	    fcp->iface[EAST] = get_ifi(i+1,j);
	    fcp->iface[SOUTH] = get_ifj(j);
	    fcp->iface[WEST] = get_ifi(i,j);
	    fcp->iface[TOP] = NULL;
	    fcp->iface[BOTTOM] = NULL;
	}
    }
    
    if (DB0) cout << "leaving bind_interfaces_to_flex_cells()" << endl;
    return SUCCESS;
}

int PistonFace:: update_west_flex_cell_after_motion()
{
    if (DB0) cout << "entering update_west_flex_cell_after_motion()" << endl;
    size_t j;
    for (j = jmin_; j <= jmax_+1; ++j) {
	vtx_[j*nvi_+1]->pos[0].x = vtx_C_.x;
    }
    for (j = jmin_; j <= jmax_; ++j) {
        update_cell_geometry_based_on_vertices(flex_cells_[j]);
	update_gas_primaries_from_extensive(flex_cells_[j], flex_cells_[j]->volume[0]);
    }
    if (DB0) cout << "leaving update_west_flex_cell_after_motion()" << endl;
    return SUCCESS;
}

int PistonFace:: update_east_flex_cell_after_motion()
{
    if (DB0) cout << "entering update_east_flex_cell_after_motion()" << endl;
    size_t j;
    for (j = jmin_; j <= jmax_+1; ++j) {
	vtx_[j*nvi_+0]->pos[0].x = vtx_C_.x;
    }
    for (j = jmin_; j <= jmax_; ++j) {
	update_cell_geometry_based_on_vertices(flex_cells_[j]);
	update_gas_primaries_from_extensive(flex_cells_[j], flex_cells_[j]->volume[0]);
    }
    if (DB0) cout << "leaving update_east_flex_cell_after_motion()" << endl;
    return SUCCESS;
}

int PistonFace:: update_flex_cell_states(double dt)
{
    size_t j;
    for (j = jmin_; j <= jmax_; ++j) {
        update_extensive_gas_time_derivatives_from_fluxes(flex_cells_[j]);
    }
    if (DB3) {    
	for (j = jmin_; j <= jmax_; ++j) {
	    cout << "j = " << j << endl;
	    print_gas_fluxes(flex_cells_[j]);
	}
    }
    for (j = jmin_; j <= jmax_; ++j) {
	update_gas_conserved_from_time_derivatives(flex_cells_[j], dt);
    }
    return SUCCESS;
}

int PistonFace:: compute_forces(vector<double> bore_resistance_f, vector<double> bore_resistance_x)
{
    size_t j;
    double local_f;
    FV_Cell *cp;

    dFg_ = 0.0;
    dR_ = 0.0;
    
    // linearly interpolate two lists of values
    linear_eval(vtx_C_.x, local_f, bore_resistance_x, bore_resistance_f);

    if (DB4) {
	cout << "x = (";
	for (size_t i=0; i<bore_resistance_x.size(); ++i) {
	    cout << bore_resistance_x[i] << ", ";
	}
	cout << ") f = (";
	for (size_t i=0; i<bore_resistance_f.size(); ++i) {
	    cout << bore_resistance_f[i] << ", ";
	}
	cout << ")" << endl;
	cout << "interpolated f = " << local_f << endl;
    }
    
    if (type_ == WEST_FACE) {
	for (j = jmin_; j <= jmax_; ++j) {
	    cp = flex_cells_[j];
	    // Calculate gas component of force
	    dFg_ += cp->fs->gas->p * ifi_[j*nii_+1]->area[0];
	}
    }
    else {
	for (j = jmin_; j <= jmax_; ++j) {
	    cp = flex_cells_[j];
	    // Calculate gas component of force
	    dFg_ += cp->fs->gas->p * ifi_[j*nii_+0]->area[0];
	}
    }

    // Calculate resistance
    dR_ = local_f;

    return SUCCESS;
}

int PistonFace:: copy_geometry_data_to_flex_cells(Block *(bd[]))
{
/** \brief 
**/    
    if (DB0) cout << "entering copy_geometry_data_to_flex_cells()" << endl;
    size_t j, indx;
    FV_Cell *fcp;

    for (j = jmin_; j <= jmax_; ++j) {
	fcp = flex_cells_[j]; 
	
	fcp->vtx[0]->pos = bd[bo_]->get_vtx(io_,j)->pos;
	fcp->vtx[1]->pos = bd[bo_]->get_vtx(io_+1,j)->pos;
	fcp->vtx[2]->pos = bd[bo_]->get_vtx(io_+1,j+1)->pos;
	fcp->vtx[3]->pos = bd[bo_]->get_vtx(io_,j+1)->pos;
    }
    
    if (DB4) cout << "setting flex-cell x-positions to " << vtx_C_.x << endl;
    // set x-positions to be flush with the piston face.
    if (type_ == WEST_FACE) {
	for (j = jmin_; j <= jmax_+1; ++j) {
	    indx = nvi_*j+1;
	    vtx_[indx]->pos[0].x = vtx_C_.x; 
	}
    }
    else {
	for (j = jmin_; j <= jmax_+1; ++j) {
	    indx = nvi_*j+0;
	    vtx_[indx]->pos[0].x = vtx_C_.x;
	}
    }
    
    if (DB2) {
	for (j = jmin_; j <= jmax_; ++j) {
	    fcp = flex_cells_[j]; 
	    cout << "flex-cell vertices:" << endl;
	    cout << "SW " << fcp->vtx[0]->pos[0].str() << endl;
	    cout << "SE " << fcp->vtx[1]->pos[0].str() << endl;
	    cout << "NE " << fcp->vtx[2]->pos[0].str() << endl;
	    cout << "NW " << fcp->vtx[3]->pos[0].str() << endl;
	}
    }

    for (j = jmin_; j <= jmax_; ++j) {
	update_cell_geometry_based_on_vertices(flex_cells_[j]);
    }
    
    if (DB0) cout << "leaving copy_geometry_data_to_flex_cells()" << endl;
    return SUCCESS;
}

int PistonFace:: copy_flow_data_to_flex_cells(Block *(bd[]))
{
    if (DB0) cout << "entering copy_flow_data_to_flex_cells()" << endl;
    size_t j;
    double uf;
    double vol_o, vol_i, total_vol;

    if( ! covers_two_cells_ ) {
	// only one grid cell is associated with the flex cell
	for (j = jmin_; j <= jmax_; ++j) {
	    flex_cells_[j]->copy_values_from(*(bd[bo_]->get_cell(io_,j)), COPY_FLOW_STATE, 0);
	    uf = bd[bo_]->get_cell(io_,j)->uf;
	    vol_o = bd[bo_]->get_cell(io_,j)->volume[0];
	    update_gas_extensive_conserved(flex_cells_[j], uf*vol_o);
	}
    }
    else {
	for (j = jmin_; j <= jmax_; ++j) {
	    // Two grid cells associated with the flex cell.
	    uf = bd[bi_]->get_cell(ii_,j)->uf;
	    vol_o = bd[bo_]->get_cell(io_,j)->volume[0];
	    vol_i = bd[bi_]->get_cell(ii_,j)->volume[0];
	    total_vol = 1.0*vol_o + uf*vol_i;
	    average_two_flow_states(flex_cells_[j],
				    bd[bo_]->get_cell(io_,j), vol_o,
				    bd[bi_]->get_cell(ii_,j), uf*vol_i);
	    update_gas_extensive_conserved(flex_cells_[j], total_vol);
	}
    }
    
    if (DB0) cout << "leaving copy_flow_data_to_flex_cells()" << endl;

    return SUCCESS;
}

int PistonFace:: copy_flow_data_to_interfaces(Block *(bd[]))
{
    if (DB0) cout << "entering copy_flow_data_to_interfaces()" << endl;
    size_t j, indx;
    FV_Interface* bif;

    // for WEST_FACE ifi
    if (type_ == WEST_FACE) {
	for (j = jmin_; j <= jmax_; ++j) {
	    indx = j*nii_ + 0;
	    bif = bd[bo_]->get_ifi(io_,j);

	    if (DB2) {
		cout << "copying data to "
		     << "WEST ifi(" << io_ << ", " << j << ") >= "
		     << "block ifi(" << bd[bo_]->imin << ", " << j << ")" << endl;
	    }
	    ifi_[indx]->copy_values_from(*bif, COPY_FLOW_STATE);
	}
    }

    // for EAST_FACE ifi
    else { 
	for (j = jmin_; j <= jmax_; ++j) {
	    indx = j*nii_ + 1;
	    bif = bd[bo_]->get_ifi(io_+1,j);
	    
	    if (DB2) {
		cout << "copying data to "
		     << "EAST ifi(" << io_+1 << ", " << j << ") <= "
		     << "block ifi(" << bd[bo_]->imax+1 << ", " << j << ")" << endl;
	    }
	    ifi_[indx]->copy_values_from(*bif, COPY_FLOW_STATE);
	}
    }

    // for both WEST_FACE and EAST_FACE ifj faces
    for (j = jmin_; j <= jmax_+1; ++j) {
	bif = bd[bo_]->get_ifj(io_,j);
	ifj_[j]->copy_values_from(*bif, COPY_FLOW_STATE);
    }
    
    if (DB0) cout << "leaving copy_flow_data_to_interfaces()" << endl;
    return SUCCESS;
}

int PistonFace:: copy_flow_data_to_grid(Block *(bd[])) 
{
    if (DB0) cout << "entering copy_flow_data_to_grid()" << endl;
    size_t j;
    FV_Cell *cp;

    for (j = jmin_; j <= jmax_; ++j) {
        cp = flex_cells_[j];

	bd[bo_]->get_cell(io_,j)->copy_values_from(*cp, COPY_FLOW_STATE, 0);
	bd[bo_]->get_cell(io_,j)->fs->gas->rho /= cp->volume[0];
	update_gas_conserved(bd[bo_]->get_cell(io_,j));

	bd[bi_]->get_cell(ii_,j)->copy_values_from(*cp, COPY_FLOW_STATE, 0);
	bd[bi_]->get_cell(ii_,j)->fs->gas->rho /= cp->volume[0];
	update_gas_conserved(bd[bi_]->get_cell(ii_,j));
    }
    
    if (DB0) cout << "leaving copy_flow_data_to_grid()" << endl;
    return SUCCESS;
}

int PistonFace:: shadow_cells(Block *(bd[])) 
{
    if (DB0) cout << "entering shadow_cells()" << endl;

    double uf = 0.0;
    
    for (size_t j = jmin_; j <= jmax_; ++j) {
	bd[bo_]->get_cell(io_,j)->status = SHADOWED_CELL;
	bd[bi_]->get_cell(ii_,j)->status = SHADOWED_CELL;
	if( covers_two_cells_ ) {
	    size_t i = ii_;
	    if (DB5) cout << "cell ii_ is partially covered." << endl;
	    if( type_ == WEST_FACE ) {
		uf = (vtx_C_.x - bd[bi_]->get_vtx(i,j)->pos[0].x) / bd[bi_]->get_cell(ii_,j)->iLength;
		if (DB5) printf("vtx_C_.x = %g, cell-x = %g, uf = %g\n", \
				vtx_C_.x, bd[bi_]->get_vtx(i,j)->pos[0].x, uf);
	    }
	    else {
		uf = (bd[bi_]->get_vtx(i+1,j)->pos[0].x - vtx_C_.x)/bd[bi_]->get_cell(ii_,j)->iLength;
		if (DB5) printf("vtx_C_.x = %g, cell-x = %g, uf = %g\n", \
				vtx_C_.x, bd[bi_]->get_vtx(i+1,j)->pos[0].x, uf);
	    }
	    bd[bo_]->get_cell(io_,j)->uf = 1.0; // redundant
	    bd[bi_]->get_cell(ii_,j)->uf = uf;
	}
	else {
	    size_t i = io_;
	    if (DB5) cout << "cell io_ is partially covered." << endl;
	    if( type_ == WEST_FACE ) {
		uf = fabs(vtx_C_.x - bd[bo_]->get_vtx(i,j)->pos[0].x)/bd[bo_]->get_cell(ii_,j)->iLength;
	    }
	    else {
		uf = fabs(vtx_C_.x - bd[bo_]->get_vtx(i+1,j)->pos[0].x)/bd[bo_]->get_cell(io_,j)->iLength;
	    }
	    bd[bo_]->get_cell(io_,j)->uf = uf;
	    bd[bi_]->get_cell(ii_,j)->uf = 0.0;
	}
    }
    
    if (DB0) cout << "leaving shadow_cells()" << endl;
    return SUCCESS;
} 

// ----------------------------------------------------------------------------
// Piston boundary condition implementation

int PistonFace:: apply_boundary_condition()
{
/**
\brief
Here we use flex_cell values of pressure and piston velocity to set the flux
vector. These are dirichlet boundary conditions for a moving wall. 
**/ 
    FV_Interface *ifp;

    if (type_ == WEST_FACE) {
	for (size_t j = jmin_; j <= jmax_; ++j) {
	    ifp = flex_cells_[j]->iface[EAST];
	    // This is the flux vector for a moving wall boundary condition
	    ifp->F->mass = 0.0;
	    ifp->F->momentum.x = flex_cells_[j]->fs->gas->p;
	    ifp->F->total_energy = flex_cells_[j]->fs->gas->p*vel_->x;
	}
    }
    else {
	for (size_t j = jmin_; j <= jmax_; ++j) {
	    ifp = flex_cells_[j]->iface[WEST];
	    // This is the flux vector for a moving wall boundary condition
	    ifp->F->mass = 0.0;
	    ifp->F->momentum.x = flex_cells_[j]->fs->gas->p;
	    ifp->F->total_energy = flex_cells_[j]->fs->gas->p*vel_->x;
	}
    }

    return SUCCESS;
}

// ----------------------------------------------------------------------------
// Piston implementation

// This variable is used as a default to set the
// vanishing distance.  We hope nobody ever simulates
// something of this dimension.
static const double VERY_LARGE_X = 1.0e6; // m

// default constructor
Piston::Piston() 
    : id_(0),
      ncf_(0),
      imin_(0), imax_(0),
      status_(P_ACTIVE),
      vanish_at_x_(VERY_LARGE_X),
      const_velocity_(false),
      postv_velocity_(false),
      dx_(0.0), diam_(1.0), area_(1.0),
      m_(1.0), acc_(0.0),
      I_(1.0), rog_(0.35), beta_(0.0), n_(0.0),
      pos_(), vel_(),
      pos0_(), vel0_()
{
    if (DB0) cout << "entering default piston constructor()" << endl;
    vtx_.resize(4);
    if (DB0) cout << "leaving default piston constructor()" << endl;
}

// copy constructor
Piston::Piston(const Piston &p) 
    : id_(p.id_),
      ncf_(p.ncf_),
      imin_(p.imin_), imax_(p.imin_),
      status_(p.status_), vanish_at_x_(p.vanish_at_x_),
      const_velocity_(p.const_velocity_),
      postv_velocity_(p.postv_velocity_),
      dx_(p.dx_), diam_(p.diam_), area_(p.area_),
      m_(p.m_), acc_(p.acc_),
      I_(p.I_), rog_(p.rog_), beta_(p.beta_), n_(p.n_),
      pos_(p.pos_), vel_(p.vel_),
      pos0_(p.pos0_), vel0_(p.vel0_),
      bore_resistance_f_(p.bore_resistance_f_),
      bore_resistance_x_(p.bore_resistance_x_)
{
    if (DB0) cout << "entering piston copy constructor" << endl;
    vtx_.resize(4);
    for (size_t i = 0; i < vtx_.size(); ++i) {
	vtx_[i] = p.vtx_[i];
    }
    wPistonFace.resize(ncf_);
    ePistonFace.resize(ncf_);
    size_t j;
    for (j = 0; j < wPistonFace.size(); ++j) {
	wPistonFace[j] = new PistonFace(*p.wPistonFace[j]);
    }
    for (j = 0; j < ePistonFace.size(); ++j) {
	ePistonFace[j] = new PistonFace(*p.ePistonFace[j]);
    }
    if (DB0) cout << "leaving piston copy constructor" << endl;
}

// default destructor
Piston::~Piston() 
{
    if (DB0) cout << "entering default piston destructor()" << endl;
    size_t j;
    for(j = 0; j < wPistonFace.size(); ++j) {
	delete wPistonFace[j];
    }
    for(j = 0; j < ePistonFace.size(); ++j) {
	delete ePistonFace[j];
    }
    if (DB0) cout << "leaving default piston destructor()" << endl;
}

// assignment operator
Piston& Piston::operator=(const Piston &p)
{ 
    if (DB0) cout << "entering piston assignment operator" << endl;
    if (this == &p) {
	if (DB0) cout << "leaving piston assignment operator" << endl;
	return *this;
    } 
    else {
	id_ = p.id_;
	ncf_ = p.ncf_;
	imin_ = p.imin_; imax_ = p.imin_;
	dx_ = p.dx_; diam_ = p.diam_; area_ = p.area_;
	m_ = p.m_; acc_ = p.acc_;
	I_ = p.I_; rog_ = p.rog_; beta_ = p.beta_; n_ = p.n_;
	pos_ = p.pos_; vel_ = p.vel_;
	pos0_ = p.pos0_; vel0_ = p.vel0_;

	status_ = p.status_;
	vanish_at_x_ = p.vanish_at_x_;
	const_velocity_ = p.const_velocity_;
	postv_velocity_ = p.postv_velocity_;

	bore_resistance_f_ = p.bore_resistance_f_;
	bore_resistance_x_ = p.bore_resistance_x_;

	vtx_.resize(4);
	for (size_t i = 0; i < vtx_.size(); ++i) {
	    vtx_[i] = p.vtx_[i];
	}
	
	for( size_t j = 0; j < wPistonFace.size(); ++j ) {
	    delete wPistonFace[j];
	}

	for( size_t j = 0; j < ePistonFace.size(); ++j ) {
	    delete ePistonFace[j];
	}

	wPistonFace.resize(ncf_);
	ePistonFace.resize(ncf_);

	size_t j;
	for (j = 0; j < wPistonFace.size(); ++j) {
	    wPistonFace[j] = new PistonFace(*(p.wPistonFace[j]));
	}
	for (j = 0; j < ePistonFace.size(); ++j) {
	    ePistonFace[j] = new PistonFace(*(p.ePistonFace[j]));
	}

	if (DB0) cout << "leaving piston assignment operator" << endl;
	return *this;
    }
}

// ----------------------------------------------------------------------------
// service functions

string Piston::string_repr()
{
    ostringstream ost;
    ost << "Piston(id=" << id_ << ", x=" << pos_.x << ", v=" << vel_.x << ")";
    return ost.str();
}

string Piston::vtk_string()
{
    ostringstream ost;
    for (size_t i = 0; i < vtx_.size(); ++i) {
	ost << vtx_[i].vtk_str() << endl;
    }
    return ost.str();
}

int Piston::write_state(FILE *fp, double t)
{
    double w_avg_pg = 0.0;
    double e_avg_pg = 0.0;
    for( size_t i = 0; i < ncf_; ++i ) {
	w_avg_pg += wPistonFace[i]->get_dFg();
	e_avg_pg += ePistonFace[i]->get_dFg();
    }
    w_avg_pg /= area_;
    e_avg_pg /= area_;
    fprintf(fp, "%d %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", 
	    static_cast<int>(id_), t, pos_.x, vel_.x, acc_, w_avg_pg, e_avg_pg );
    return SUCCESS;
}

int Piston::read_state(FILE *fp)
{
/// \brief Attempt to read piston state information from fp
    char line[132];
    double t, x, u, acc, wpg, epg;
    unsigned int temporary_id, nread;
    if ( fgets( line, 131, fp ) == NULL ) {
	cerr << "Piston::read_state(): failure of fgets()" << endl;
	cerr << "Quitting program." << endl;
	exit(FILE_ERROR);
    }
    nread = sscanf(line, "%u %lf %lf %lf %lf %lf %lf", &temporary_id, &t, &x, &u, &acc, &wpg, &epg);
    if (nread != 4) {
	cout << "Piston::read_state(): id=" << id_ 
	     << " failed to read all values." << endl;
    } 
    if (temporary_id != id_) {
	cout << "Piston::read_state(): id=" << id_ 
	     << " id in file does not match " << temporary_id << endl;
    } 
    pos_.x = x;
    vel_.x = u;
    update_piston_geometry_based_on_position();
    return SUCCESS;
}

int Piston::save_state()
{
    vel0_ = vel_;
    pos0_ = pos_;
    return SUCCESS;
}

int Piston::restore_state()
{
    vel_ = vel0_;
    pos_ = pos0_;
    return SUCCESS;
}

// ----------------------------------------------------------------------------
// member functions

int Piston::set_values(size_t id,
		       double D,
		       double L,
		       double m,
		       double x,
		       double u,
		       double rifling_twist,
		       double rog,
		       double vanish_at_x,
		       vector<double> bore_resistance_f,
		       vector<double> bore_resistance_x)
{
    global_data &G = *get_global_data_ptr();
    if (DB0) cout << "entering set_values()" << endl;

    id_ = id;
    n_ = rifling_twist;
    vanish_at_x_ = vanish_at_x;

    diam_ = D;
    // Once D is known, set area_
    if ( G.axisymmetric ) {
	area_ = (M_PI*D*D/4.0)/(2.0*M_PI); // area per radian
    }
    else {
	area_ = D; // area per unit metre
    }

    if (rifling_twist > 0.0)
	beta_ = atan(rifling_twist/M_PI);
    
    if (DB4) cout << "twist angle = " << beta_*180.0/M_PI << endl;

    m_ = m;
    if (rog > 0.0)
	I_ = m*rog*rog;
    
    pos_.x = x;
    vel_.x = u;

    vtx_.resize(4);
    vtx_[0] = Vector3(x-(L/2.0), 0.0, 0.0);
    vtx_[1] = Vector3(x+(L/2.0), 0.0, 0.0);

    if ( G.axisymmetric ) {
	vtx_[2] = Vector3(x+(L/2.0), D/2.0, 0.0);
	vtx_[3] = Vector3(x-(L/2.0), D/2.0, 0.0);
    }
    else {
	vtx_[2] = Vector3(x+(L/2.0), D, 0.0);
	vtx_[3] = Vector3(x-(L/2.0), D, 0.0);
    }

    bore_resistance_f_ = bore_resistance_f;
    bore_resistance_x_ = bore_resistance_x;

    if (DB0) cout << "leaving set_values()" << endl;
    return SUCCESS;
}

int Piston::initialise_from_config_file(global_data &G,
					Block *(bd[]),
					ConfigParser &cfg, 
					string section)
{
    cout << "Entering piston initialise from config...\n";
    size_t id;
    if( ! cfg.parse_size_t(section, "id", id, 0) ) {
	cout << "Error reading id in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    double D;
    if( ! cfg.parse_double(section, "D", D, 0.0) ) {
	cout << "Error reading D in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    if( D <= 0.0 ) {
	cout << "Error reading D in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "A value greater than 0.0 is expected.  Value received: " << D << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    double L;
    if( ! cfg.parse_double(section, "L", L, 0.0) ) {
	cout << "Error reading L in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    if( L <= 0.0 ) {
	cout << "Error reading L in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "A value greater than 0.0 is expected.  Value received: " << L << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    double m;
    if( ! cfg.parse_double(section, "m", m, 1.0) ) {
	cout << "Error reading m in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    if( m <= 0.0 ) {
	cout << "Error reading m in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "A value greater than 0.0 is expected.  Value received: " << m << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    double x0;
    if( ! cfg.parse_double(section, "x0", x0, 0.0) ) {
	cout << "Error reading x0 in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    double v0;
    if( ! cfg.parse_double(section, "v0", v0, 0.0) ) {
	cout << "Error reading v0 in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    double rifling_twist;
    if( ! cfg.parse_double(section, "rifling_twist", rifling_twist, 0.0) ) {
	cout << "Error reading rifling_twist in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    double rog;
    if( ! cfg.parse_double(section, "rog", rog, D/pow(8.0, 0.5)) ) {
	cout << "Error reading rog in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    double vanish_at_x;
    if( ! cfg.parse_double(section, "vanish_at_x", vanish_at_x, VERY_LARGE_X) ) {
	cout << "Error reading vanish_at_x in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    vector<double> bore_resistance_x;
    vector<double> not_found;
    not_found.push_back(0.0);
    if( ! cfg.parse_vector_of_doubles(section, "bore_resistance_x", bore_resistance_x, not_found) ) {
	cout << "Error reading bore_resistance_x in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    if (bore_resistance_x.empty()) {
	cout << "Bore_resistance_x in section: " << section << " is empty" << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
        exit(BAD_INPUT_ERROR);
    }
    
    vector<double> bore_resistance_f;
    if( ! cfg.parse_vector_of_doubles(section, "bore_resistance_f", bore_resistance_f, not_found)) {
	cout << "Error reading bore_resistance_f in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    if (bore_resistance_f.empty()) {
	cout << "bore_resistance_f in section: " << section << " is empty" << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
        exit(BAD_INPUT_ERROR);
    }

    if (bore_resistance_x.size() != bore_resistance_f.size()) {
	cout << "Input error in section: " << section << endl 
	     << "of input file: " << cfg.file_name << endl
	     << "Equal vector sizes expected.  Values received: " 
	     << bore_resistance_x.size() << " "
	     << bore_resistance_f.size() << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    string const_v_flag;
    if( ! cfg.parse_string(section, "constant_velocity_flag", const_v_flag, "False") ) {
	cout << "Error reading constant_velocity_flag in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    if( const_v_flag == "True" )
	const_velocity_ = true;
    else
	const_velocity_ = false;

    string postv_v_flag;
    if( ! cfg.parse_string(section, "positive_velocity_flag", postv_v_flag, "False") ) {
	cout << "Error reading positive_velocity_flag in section: " << section << endl
	     << "of input file: " << cfg.file_name << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    if( postv_v_flag == "True" )
	postv_velocity_ = true;
    else
	postv_velocity_ = false;

    int status;
    
    cout << "Attempting to set piston values from config file..." << endl;
    status = set_values(id, D, L, m, x0, v0, rifling_twist, rog, vanish_at_x, bore_resistance_f, bore_resistance_x);
    if (status != SUCCESS) {
	cout << "Failure!" << endl;
	return FAILURE;
    }
    else cout << "Success!" << endl;

    cout << "Attempting to set pistonface values from block data..." << endl;
    status = initialise_piston_faces(G, bd);
    if (status != SUCCESS) {
	cout << "Failure!" << endl;
	return FAILURE;
    }
    else cout << "Success!" << endl;

    return SUCCESS;
}

int Piston::initialise_piston_faces(global_data &G, Block *(bd[]))
{
    if (DB0) cout << "entering initialise_piston_faces()" << endl;
    size_t jb;

    Vector3 vtxL, vtxR;
    
    if (locate_west_block_index(jb, vtx_[0], G, bd) != SUCCESS) return FAILURE ;

#if 0
    // FIX-ME
    for (; jb != -1; jb = bd[jb]->bcp[NORTH]->neighbour_block) {
	vtxL = Vector3(vtx_[0].x, bd[jb]->get_vtx(bd[jb]->imin,bd[jb]->jmin)->pos[0].y, 0.0);
	vtxR = Vector3(vtx_[0].x, bd[jb]->get_vtx(bd[jb]->imin,bd[jb]->jmax+1)->pos[0].y, 0.0);

	if( vtxL.y >= vtx_[3].y ) {
	    // We don't need to continue making faces.
	    break;
	}

	wPistonFace.push_back(new PistonFace());
	wPistonFace.back()->set_values(-1.0, vtxL, vtxR, &vel_, WEST_FACE);
	wPistonFace.back()->set_indices_for_flex_cells(*bd[jb]);
	wPistonFace.back()->allocate_memory();
    }
#endif

    // Piston vertices must lie within block
    if (vtx_[3].y > vtxR.y) {
	cout << "Error in initialise_piston_faces()\n";
	cout << "Entire face must lie within the grid\n";
	cout << "No block found for node " << vtxR << endl;

	return FAILURE;
    }
    
    if (locate_east_block_index(jb, vtx_[1], G, bd) != SUCCESS) return FAILURE ;

#if 0
    // FIX-ME
    for (; jb != -1; jb = bd[jb]->bcp[NORTH]->neighbour_block) {
	vtxL = Vector3(vtx_[1].x, bd[jb]->get_vtx(bd[jb]->imax+1,bd[jb]->jmax+1)->pos[0].y, 0.0);
	vtxR = Vector3(vtx_[1].x, bd[jb]->get_vtx(bd[jb]->imax+1,bd[jb]->jmin)->pos[0].y, 0.0);

	if( vtxR.y >= vtx_[2].y ) {
	    // We don't need to continue making faces.
	    break;
	}

	ePistonFace.push_back(new PistonFace());
	ePistonFace.back()->set_values(+1.0, vtxL, vtxR, &vel_, EAST_FACE);
	ePistonFace.back()->set_indices_for_flex_cells(*bd[jb]);
	ePistonFace.back()->allocate_memory();
    }
#endif

    // Piston vertices must lie within block
    if (vtx_[2].y > vtxL.y) {
	cout << "Error in initialise_piston_faces()\n";
	cout << "Entire face must lie within the grid\n";
	cout << "No block found for node " << vtxL << endl;

	return FAILURE;
    }
    
    update_piston_geometry_based_on_vertices();
    
    // Piston face vectors must be of equal size
    if (wPistonFace.size() != ePistonFace.size()) {
	cout << "Error in initialise_piston_faces()\n";
	cout << "Piston" << id_ << " face vectors of unequal length.\n";
	cout << "West length = " << wPistonFace.size() << endl;
	cout << "East length = " << ePistonFace.size() << endl;
    }
    ncf_ = wPistonFace.size();
    bmin_.resize(ncf_, 0);
    bmax_.resize(ncf_, 0);

    return SUCCESS;
}

int Piston::change_of_state_due_to_motion(global_data &G, Block *(bd[]), double dt)
{
    int status;

    //cout << "======== 1. copy_flow_data_to_piston() ============\n";
    copy_data_to_piston(bd);
    if (DB5) print_flex_cell_flow_states();

    //cout << "======== 2. update_flex_cell_states() ============\n";
    update_flex_cell_states(dt);
    //if (DB5) print_flex_cell_flow_states();

    //cout << "======== 3. update_piston_state() ================\n";
    status = update_piston_state(G, bd, dt);
    //if (DB5) print_flex_cell_flow_states();

    //cout << "======== 4. copy_flow_data_to_grid() =============\n";
    copy_flow_data_to_grid(bd);
    if (DB5) print_flex_cell_flow_states();

    if (status == FAILURE) {
	cout << "returning FAILURE status" << endl;
	return PISTON_FAILURE;
    }
    else if (status == END_OF_BLOCK) {
	cout << "returning END_OF_BLOCK status" << endl;
	return END_OF_BLOCK;
    }
    else {
	return SUCCESS;
    }
}

int Piston::print_flex_cell_flow_states()
{
    FV_Cell *cp;
    cout << "WEST PistonFace" << endl;
    for (size_t jp = 0; jp < ncf_; ++jp) {
	for (size_t j = wPistonFace[jp]->get_jmin(); j <= wPistonFace[jp]->get_jmax(); ++j) {
	    cout << "index-j = " << j << endl;
	    cp = wPistonFace[jp]->get_flex_cell(j);
	    print_flow_state(cp);
	}
    }
    cout << "EAST PistonFace" << endl;
    for (size_t jp = 0; jp < ncf_; ++jp) {
	for (size_t j = ePistonFace[jp]->get_jmin(); j <= ePistonFace[jp]->get_jmax(); ++j) {
	    cout << "index-j = " << j << endl;
	    cp = ePistonFace[jp]->get_flex_cell(j);
	    print_flow_state(cp);
	}
    }
    return SUCCESS;
}

static const double allowable_movement_fraction = 0.5;
static const double GRAVITY = 9.81;

int Piston::update_piston_state(global_data &G, Block *(bd[]), double dt)
{
/** 
Attempt the piston dynamics subproblem
 1. Use the averaged face pressure values to find net force
 2. Apply the effect of resistance.
 3. Find and apply the new piston velocity and displacement

**/
    if (DB0) cout << "entering update_piston_state()" << endl;

    size_t jp;
    double dF, dR, m_eff;
    
    for (jp = 0; jp < ncf_; ++jp) {
	wPistonFace[jp]->compute_forces(bore_resistance_f_, bore_resistance_x_);
	ePistonFace[jp]->compute_forces(bore_resistance_f_, bore_resistance_x_);
    }

    dF = dR = 0.0;
    for (jp = 0; jp < ncf_; ++jp) {
	dF += wPistonFace[jp]->get_dF() - ePistonFace[jp]->get_dF();
    }

    if ( G.axisymmetric ) {
	dF *= 2.0*M_PI;
    }
    
    dR = wPistonFace[0]->get_dR()*m_*GRAVITY;

    // WARNING: the piston will oscillate rather than stopping inside the chamber
    if (fabs(dF) > dR || fabs(vel_.x) > 0.0) dF = dF - SIGN(vel_.x)*dR; 
    else dF = 0.0;

    if (beta_ == 0)
	m_eff = m_;
    else
	m_eff = (m_ + (4*M_PI*I_*tan(beta_)/(n_*diam_*diam_)));
	  
    if (!const_velocity_) {
	acc_ = dF/m_eff;
	dx_ = vel_.x*dt + 0.5*acc_*dt*dt;
	vel_.x += acc_*dt;
    }
    else {
	acc_ = 0.0;
        dx_ = vel_.x*dt;
    }
    
    if (postv_velocity_) {
	// forbid negative velocities
	if (vel_.x < 0.0) {
	    acc_ = 0.0;
	    vel_.x = 0.0;
	}
    }
    
    if (DB4) {
	cout << setprecision(6) 
	     << dF << "N over " << ncf_ << " component face(s) with " 
	     << dR << "N resistance over " << dt << "s results in:\n"
	     << "acc     = " << acc_ << "m/s**2\n"
	     << "vel.x   = " << vel_.x  << "m/s\n"
	     << "dx      = " << dx_  << "m\n";
    }

    // Check that piston did not move too far
    // "Too far" is a criterion based on allowable_movement_fraction

    double cell_width = 0.0;
    // 1. Check at west end.
    if( wPistonFace[0]->covers_two_cells() ) {
	// we're interested in ii_
	cell_width = bd[get_wbi()]->get_cell(get_wii(),wPistonFace[0]->get_jmin())->iLength;
    }
    else {
	// we're interest in io_
	cell_width = bd[get_wbo()]->get_cell(get_wio(),wPistonFace[0]->get_jmin())->iLength;
    }
    
    if( fabs(dx_) >= allowable_movement_fraction*cell_width ) {
	cout << "Error in update_piston_state()\n";
	cout << "Piston" << id_ << " movement too large: dx > " << allowable_movement_fraction << "*cell length\n";
	cout << fabs(dx_) << " >= " << allowable_movement_fraction*cell_width << endl;
	cout << "This occurred at the west piston face.\n";
	cout << "Consider limiting timestep based on piston speed.\n";
	exit(PISTON_FAILURE);
    }

    // 2. Check at east end
    if( ePistonFace[0]->covers_two_cells() ) {
	// we're interested in ii_
	cell_width = bd[get_ebi()]->get_cell(get_eii(),ePistonFace[0]->get_jmin())->iLength;
    }
    else {
	// we're interest in io_
	cell_width = bd[get_ebo()]->get_cell(get_eio(),ePistonFace[0]->get_jmin())->iLength;
    }
    
    if( fabs(dx_) >= allowable_movement_fraction*cell_width ) {
	cout << "Error in update_piston_state()\n";
	cout << "Piston" << id_ << " movement too large: dx > " << allowable_movement_fraction << "*cell length\n";
	cout << dx_ << " >= " << allowable_movement_fraction*cell_width << endl;
	cout << "This occurred at the east piston face.\n";
	cout << "Consider limiting timestep based on piston speed.\n";
	exit(PISTON_FAILURE);
    }

    // -- At this point, we've passed the allowable movement test.

    // ------------------------------------------------------------------------
    // Update piston geometry
    //
    vtx_[0].x += dx_;
    vtx_[1].x += dx_;
    vtx_[2].x += dx_;
    vtx_[3].x += dx_;

    update_piston_geometry_based_on_vertices();

    // After piston update, check whether or not the
    // piston should still be active.
    // If it is set to inactive, this will NOT affect the
    // current step --- it will continue to influence the
    // flow field as previous.  In the NEXT step the piston
    // will be ignored entirely.  Caveat, once "ignored",
    // the piston will be ignored for the remainder of the 
    // simulation.
    if( pos_.x > vanish_at_x_ ) {
	status_ = P_INACTIVE;
	cout << "===================================================================\n";
	cout << " The centroid of Projectile-" << id_ << " has exceeded x= " << vanish_at_x_ << endl;
	cout << " and as such will be ignored for the remainder of the simulation.\n";
	cout << "===================================================================\n";
	
	// Activation of cells underneath the piston will occur 
	// during the next call to configure_cells().
	// If there were multiple projectiles, issuing this 
	// command here would remove all of them.
    }

    // change west flex cell geometry by dx
    for (jp = 0; jp < ncf_; ++jp) {
	wPistonFace[jp]->update_west_flex_cell_after_motion();
    }
    // change east flex cell geometry by dx
    for (jp = 0; jp < ncf_; ++jp) {
	ePistonFace[jp]->update_east_flex_cell_after_motion();
    }

    // ------------------------------------------------------------------------

    if (DB1) {
	size_t jmin = wPistonFace[0]->get_jmin();
	FV_Cell *cell = wPistonFace[0]->get_flex_cell(jmin);
	cout << "dt, xpos, WESTx, EASTx, iLength, area, volume = " 
	     << dt << ", "
	     << cell->pos[0].x << ", " 
	     << cell->vtx[3]->pos[0].x << ", "
	     << cell->vtx[2]->pos[0].x << ", "
	     << cell->iLength << ", " 
	     << cell->area[0] << ", "
	     << cell->volume[0]
	     << endl;
    }

    return SUCCESS;
}

static const size_t MIN_CELLS = 1;

int Piston::configure_cells(global_data &G, Block *(bd[])) 
{
    if (DB0) cout << "entering configure_cells()" << endl;

    size_t jp;
    size_t nicells = 0;

    // north-south boundary conditions - to be implemented

//     // 0. Loop over all blocks re-setting the index range on
//     //    piston interface boundary conditions
//     for( size_t jb = 0; jb < int(bd.size()); ++jb ) {
// 	if( bd[jb].get_south_bc_type() == "piston_interface" ) {
// 	    // By setting the x-range in reverse, the normal exchange boundary
// 	    // condition will be used.
// 	    bd[jb].set_x_index_range_south_bc(bd[jb].get_imax(), bd[jb].get_imin());
// 	}
//     }

    for (jp = 0; jp < ncf_; ++jp) {
	if (ePistonFace[jp]->locate_shadow_indices(G, bd) != SUCCESS) return FAILURE;
	if (wPistonFace[jp]->locate_shadow_indices(G, bd) != SUCCESS) return FAILURE;
    }

    for (jp = 0; jp < ncf_; ++jp) {
	// Case A: Piston faces lie in the same block
	if (wPistonFace[jp]->get_bi() == ePistonFace[jp]->get_bi()) {
	    if (DB4) cout << "Piston" << id_ << " face component " << jp << " lies in block " 
			  << ePistonFace[jp]->get_bi() << endl;
	    
	    imin_ = wPistonFace[jp]->get_ii()+1;
	    imax_ = ePistonFace[jp]->get_ii()-1;
	    bmin_[jp] = bmax_[jp] = wPistonFace[jp]->get_bi();
	    
	    // north-south boundary conditions, yet to be implemented

// 	    if( jp == ncf_ - 1) {
// 		// When working on the north most interface.
// 		size_t btop = bmin_[jp];
// 		// If there is a block above, then we
// 		// need to play with its boundary condition.
// 		if( G.con_matrix[btop,NORTH].index != -1 ) {
// 		    size_t bnorth = G.con_matrix[btop,NORTH].index;
// 		    // Assume the south of the adjoining block 
// 		    // connects to the north when using a Piston.
// 		    // Otherwise, the cases get too complex to handle.
// 		    if( G.con_matrix[btop,NORTH].bndry != SOUTH ) {
// 			cout << "Error trying to set the Piston_interface_boundary_condition\n"
// 			     << "on an ajdoining block north of the projectile.  There is a \n"
// 			     << "hard-coded assumption that the SOUTH edge joins to the NORTH\n"
// 			     << "edge of the block containing the projectile.\n"
// 			     << "BAILING OUT\n";
// 			exit(NOT_IMPLEMENTED_ERROR);
// 		    }
		    
// 		    if( bd[bnorth].get_south_bc_type() != "piston_interface" ) {
// 			bd[bnorth].set_south_bc_to_piston_interface();
// 		    }
		    
// 		    bd[bnorth].set_x_index_range_south_bc(wPistonFace[jp]->get_io(),
// 							  ePistonFace[jp]->get_io());
// 		}
// 	    }
	}
	// Case B: Piston faces lie in different blocks
	else {
	    if (DB4) cout << "Piston" << id_ << " W and E face-component-" << jp << "\'s lie in blocks " 
		 << wPistonFace[jp]->get_bi() << " and " << ePistonFace[jp]->get_bi() << " respectively\n";
	    imin_ = wPistonFace[jp]->get_ii()+1;
	    imax_ = ePistonFace[jp]->get_ii()-1;
	    bmin_[jp] = wPistonFace[jp]->get_bi();
	    bmax_[jp] = ePistonFace[jp]->get_bi();
	    // Check imin_ is appropriate.
	    if( imin_ > bd[bmin_[jp]]->imax ) {
		// Shuffle interior index to next block
		imin_ = bd[bmax_[jp]]->imin;
		bmin_[jp] = bmax_[jp];
	    }
	    if( imax_ < bd[bmax_[jp]]->imin ) {
		imax_ = bd[bmin_[jp]]->imax;
		bmax_[jp] = bmin_[jp];
	    }
// 	    if( jp == ncf_ - 1) {
// 		// When working on the north most interface.
// 		// First consider left block.
// 		size_t btop = bmin_[jp];
// 		// If there is a block above, then we
// 		// need to play with its boundary condition.
// 		if( G.con_matrix[btop,NORTH].index != -1 ) {
// 		    size_t bnorth = G.con_matrix[btop,NORTH].index;
// 		    // Assume the south of the adjoining block 
// 		    // connects to the north when using a Piston.
// 		    // Otherwise, the cases get too complex to handle.
// 		    if( G.con_matrix[btop,NORTH].bndry != SOUTH ) {
// 			cout << "Error trying to set the Piston_interface_boundary_condition\n"
// 			     << "on an ajdoining block north of the projectile.  There is a \n"
// 			     << "hard-coded assumption that the SOUTH edge joins to the NORTH\n"
// 			     << "edge of the block containing the projectile.\n"
// 			     << "BAILING OUT\n";
// 			exit(NOT_IMPLEMENTED_ERROR);
// 		    }
		    
// 		    if( bd[bnorth].get_south_bc_type() != "piston_interface" ) {
// 			bd[bnorth].set_south_bc_to_piston_interface();
// 		    }
		    
// 		    bd[bnorth].set_x_index_range_south_bc(wPistonFace[jp]->get_io(),
// 							  bd[bnorth].get_imax());
// 		}
		
// 		// Then consider right block.
// 		btop = bmax_[jp];
// 		// If there is a block above, then we
// 		// need to play with its boundary condition.
// 		if( G.con_matrix[btop,NORTH].index != -1 ) {
// 		    size_t bnorth = G.con_matrix[btop,NORTH].index;
// 		    // Assume the south of the adjoining block 
// 		    // connects to the north when using a Piston.
// 		    // Otherwise, the cases get too complex to handle.
// 		    if( G.con_matrix[btop,NORTH].bndry != SOUTH ) {
// 			cout << "Error trying to set the Piston_interface_boundary_condition\n"
// 			     << "on an ajdoining block north of the projectile.  There is a \n"
// 			     << "hard-coded assumption that the SOUTH edge joins to the NORTH\n"
// 			     << "edge of the block containing the projectile.\n"
// 			     << "BAILING OUT\n";
// 			exit(NOT_IMPLEMENTED_ERROR);
// 		    }
		    
// 		    if( bd[bnorth].get_south_bc_type() != "piston_interface" ) {
// 			bd[bnorth].set_south_bc_to_piston_interface();
// 		    }
		    
// 		    bd[bnorth].set_x_index_range_south_bc(bd[bnorth].get_imin(),
// 							  ePistonFace[jp]->get_io());
// 		}
// 	    }
	}
    }
    
    for (jp = 0; jp < ncf_; ++jp) {
	if (bmin_[jp] == bmax_[jp]) {
	    nicells = imax_ - imin_ + 1;
	}
	else {
	    nicells = (bd[bmin_[jp]]->imax - imin_ + 1) + (imax_ - bd[bmax_[jp]]->imin + 1);
	}
    }
    
    if (nicells < MIN_CELLS) {
	cout << "The projectile location in the grid does not cover a minimum of " 
	     << MIN_CELLS << " cells." << endl
	     << "It only covers " <<  nicells << endl;
	cout << "Two suggestions:\n";
	cout << "   1. Use a longer projectile (which covers more cells)\n";
	cout << "   2. Use a finer grid (which provides more cells under the projectile)\n";
	return FAILURE;
    }
    
    // Loop over projectiles: set shadowed, then deactivate
    for (jp = 0; jp < ncf_; ++jp) {
	wPistonFace[jp]->shadow_cells(bd);
	ePistonFace[jp]->shadow_cells(bd);
    }
    deactivate_cells(bd);

    if (DB4) {
	//Print information regarding piston location
	cout << "Piston" << id_ << " exists within x-range = " 
	     << setprecision(5) << showpoint
	     << wPistonFace[0]->get_pos()->x << " : " 
	     << ePistonFace[0]->get_pos()->x << endl;
	
	cout << "Piston" << id_ << " covers cell-index ranges\n";
	for( jp = 0; jp < ncf_; ++jp ) {
	    cout << "For face component: " << jp << endl;
	    if( bmin_[jp] == bmax_[jp] ) {
		cout << " block: " << bmin_[jp] << endl;
		cout << " i-range: " << imin_ << " : " << imax_ << endl;
		cout << " j-range: " << wPistonFace[jp]->get_jmin() << " : " 
		     << wPistonFace[jp]->get_jmax() << endl;
	    }
	    else {
		cout << "block-range: " << bmin_[jp] << " : " << bmax_[jp] << endl;
		cout << " i-range: " << imin_ << " : " << bd[bmin_[jp]]->imax
		     << " (in block: " << bmin_[jp] << "),";
		cout << bd[bmax_[jp]]->imin << " : " << imax_
		     << " (in block: " << bmax_[jp] << ")" << endl;
		cout << " j-range: " << wPistonFace[jp]->get_jmin() << " : " 
		     << wPistonFace[jp]->get_jmax()
		     << " (in block: " << bmin_[jp] << ")" << endl;
	    }
	}
    }

    if (DB0) cout << "leaving configure_cells()" << endl;
    return SUCCESS;
}

int Piston::deactivate_cells(Block *(bd[])) 
{
    for (size_t jp = 0; jp < ncf_; ++jp) {
	if (bmin_[jp] == bmax_[jp]) {
	    for (size_t j = wPistonFace[jp]->get_jmin(); j <= wPistonFace[jp]->get_jmax(); ++j) {
		for (size_t i = imin_; i <= imax_; ++i) {
		    bd[bmin_[jp]]->get_cell(i,j)->status = MASKED_CELL;
		}
	    }
	}
	else {
	    for (size_t j = wPistonFace[jp]->get_jmin(); j <= wPistonFace[jp]->get_jmax(); ++j) {
		for (size_t i = imin_; i <= bd[bmin_[jp]]->imax; ++i) {
		    bd[bmin_[jp]]->get_cell(i,j)->status = MASKED_CELL;
		}
		for (size_t i = bd[bmax_[jp]]->imin; i <= imax_; ++i) {
		    bd[bmax_[jp]]->get_cell(i,j)->status = MASKED_CELL;
		}
	    }
	}
    }
    
    return SUCCESS;
}

int Piston::update_flex_cell_states(double dt)
{
    if (DB0) cout << "entering update_flex_cell_states()" << endl;
    for (size_t jp = 0; jp < ncf_; ++jp) {
	wPistonFace[jp]->apply_boundary_condition();
 	ePistonFace[jp]->apply_boundary_condition();
    }

    for (size_t jp = 0; jp < ncf_; ++jp) {
	wPistonFace[jp]->update_flex_cell_states(dt);
	ePistonFace[jp]->update_flex_cell_states(dt);
    }
    
    if (DB0) cout << "leaving update_flex_cell_states()" << endl;
    return SUCCESS;
}

int Piston::update_piston_geometry_based_on_vertices() 
{
    for (size_t jp = 0; jp < ncf_; ++jp) {
	wPistonFace[jp]->move_face(vtx_[0]);
	ePistonFace[jp]->move_face(vtx_[1]);
    }
    pos_.x =  (vtx_[0].x + vtx_[1].x)/2.0;

    return SUCCESS;
}

int Piston::update_piston_geometry_based_on_position()
{
    double length = vtx_[1].x - vtx_[0].x;
    double x_left = pos_.x - length/2.0;
    double x_right = pos_.x + length/2.0;

    vtx_[0].x = x_left;
    vtx_[1].x = x_right;
    vtx_[2].x = x_right;
    vtx_[3].x = x_left;

    for (size_t jp = 0; jp < ncf_; ++jp) {
	wPistonFace[jp]->move_face(vtx_[0]);
	ePistonFace[jp]->move_face(vtx_[1]);
    }

    return SUCCESS;
}

int Piston::copy_data_to_piston(Block *(bd[]))
{
/**
Copy grid flow data from deactivated cells to the pistonFace cells and 
interfaces.

**/
    if (DB0) cout << "entering copy_data_to_piston()" << endl;
    for (size_t jp = 0; jp < ncf_; ++jp) {
	wPistonFace[jp]->copy_geometry_data_to_flex_cells(bd);
	ePistonFace[jp]->copy_geometry_data_to_flex_cells(bd);

	wPistonFace[jp]->copy_flow_data_to_flex_cells(bd);
	ePistonFace[jp]->copy_flow_data_to_flex_cells(bd);
	
	wPistonFace[jp]->copy_flow_data_to_interfaces(bd);
	ePistonFace[jp]->copy_flow_data_to_interfaces(bd);
    }

    if (DB0) cout << "leaving copy_data_to_piston()" << endl;
    return SUCCESS;
}

int Piston::copy_flow_data_to_grid(Block *(bd[]))
{
    if (DB0) cout << "entering copy_flow_data_to_grid()" << endl;
    FV_Cell *cp;

    for (size_t jp = 0; jp < ncf_; ++jp) {
	if (bmin_[jp] == bmax_[jp]) {
	    for (size_t j = wPistonFace[jp]->get_jmin(); j <= wPistonFace[jp]->get_jmax(); ++j) {
		cp = wPistonFace[jp]->get_flex_cell(j);
		for (size_t i = imin_; i <= imax_; ++i) {
		    bd[bmin_[jp]]->get_cell(i,j)->copy_values_from(*cp, COPY_FLOW_STATE, 0);
		    update_gas_conserved(bd[bmin_[jp]]->get_cell(i,j));
		}
	    }
	}
	else {
	    for (size_t j = wPistonFace[jp]->get_jmin(); j <= wPistonFace[jp]->get_jmax(); ++j) {
		cp = wPistonFace[jp]->get_flex_cell(j);
		for (size_t i = imin_; i <= bd[bmin_[jp]]->imax; ++i) {
		    bd[bmin_[jp]]->get_cell(i,j)->copy_values_from(*cp, COPY_FLOW_STATE, 0);
		    update_gas_conserved(bd[bmin_[jp]]->get_cell(i,j));
		}
		cp = ePistonFace[jp]->get_flex_cell(j);
		for (size_t i = bd[bmax_[jp]]->imin; i <= imax_; ++i) {
		    bd[bmax_[jp]]->get_cell(i,j)->copy_values_from(*cp, COPY_FLOW_STATE, 0);
		    update_gas_conserved(bd[bmax_[jp]]->get_cell(i,j));
		}
	    }
	}
    }
    
    // copy shadow cell data to grid
    for (size_t jp = 0; jp < ncf_; ++jp) {
	wPistonFace[jp]->copy_flow_data_to_grid(bd);
	ePistonFace[jp]->copy_flow_data_to_grid(bd);
    }

    if (DB0) cout << "leaving copy_flow_data_to_grid()" << endl;
    return SUCCESS;
}


double Piston::compute_signal_frequency()
{
    size_t jmin = wPistonFace[0]->get_jmin();

    // Check all signals and find largest.
    double u_mag, a, L; 
    FV_Cell *cpw, *cpe;
    
    u_mag = vabs(vel_.x);
    cpw = wPistonFace[0]->get_flex_cell(jmin);
    cpe = ePistonFace[0]->get_flex_cell(jmin);

    if (cpw->fs->gas->a > cpe->fs->gas->a) {
	a = cpw->fs->gas->a;
	L = cpw->L_min;
    }
    else {
	a = cpe->fs->gas->a;
	L = cpe->L_min;
    }

    // ASSERT(!isnan(u_mag) && !isinf(u_mag));
    // ASSERT(!isnan(a) && !isinf(a));

    double signal = (u_mag + a) / L;
    double signal_max = signal;

    // ASSERT(!isnan(signal_max) && !isinf(signal_max));
    
    return signal_max;
}

// ----------------------------------------------------------------------------
// Search functions

int activate_cells(global_data &G, Block *(bd[])) 
{
    string failure_message("error in activate_cells()");
    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	for ( FV_Cell *cp: G.my_blocks[jb]->active_cells ) set_cell_status(cp, NORMAL_CELL);
    }
    return SUCCESS;
}

int print_cell_status( FV_Cell *c ) 
{
    cout << "status = " << c->status << endl;
    return SUCCESS;
}

int print_all_cells_status(global_data &G, Block *(bd[]))
{
    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	for ( FV_Cell *cp: G.my_blocks[jb]->active_cells ) print_cell_status(cp);
    }
    return SUCCESS;
}

int locate_west_block_index(size_t &b_index, Vector3 &pos, global_data &G, Block *(bd[]))
{
    // Handles the output of locate_block_index 
    // with knowledge of face type
    // Returns one of two cases:
    //     Case A: corrected point lies within domain
    //     Case B: point lies outside domain

    size_t j;
    int next_block;
    int status;
    
    status = locate_block_index(b_index, pos, bd, G.nblock);
    if (status == END_OF_BLOCK) {
	for (j = 0; j < G.nblock; ++j) {
	    if (fabs(pos.x - bd[j]->get_vtx(bd[j]->imin,bd[j]->jmin)->pos[0].x) < WIDTH_TOL) {
		if (pos.y >= bd[j]->get_vtx(bd[j]->imin,bd[j]->jmin)->pos[0].y) {
		    b_index = j;
		}
	    }
	}
// 	cout << "Piston found at the west boundary of block " << b_index << endl;

	next_block = bd[b_index]->bcp[WEST]->neighbour_block;
// 	cout << "West block index = " << next_block << endl;
	if (next_block == -1) status = FAILURE;
	else b_index = static_cast<size_t>(next_block);
    }
    
    if (status == FAILURE) {
	cout << "Error in initialise_piston_faces().\n";
	cout << "Position exceeds grid.\n";
	cout << pos << endl;
	cout << "Check piston geometry, blocks or block connections.\n";
	
	return FAILURE;
    }

    return SUCCESS;
}

int locate_east_block_index(size_t &b_index, Vector3 &pos, global_data &G, Block *(bd[]))
{
    // Handles the output of locate_block_index 
    // with knowledge of face type
    // Returns one of two cases:
    //     Case A: corrected point lies within domain
    //     Case B: point lies outside domain

    size_t j;
    int next_block;
    int status;
    
    status = locate_block_index(b_index, pos, bd, G.nblock);
    if (status == END_OF_BLOCK) {
	for (j = 0; j < G.nblock; ++j) {
	    if (fabs(pos.x - bd[j]->get_vtx(bd[j]->imax+1,bd[j]->jmin)->pos[0].x) < WIDTH_TOL) {
		if (pos.y >= bd[j]->get_vtx(bd[j]->imax+1,bd[j]->jmin)->pos[0].y) {
		    b_index = j;
		}
	    }
	}
// 	cout << "Piston found at the east boundary of block " << b_index << endl;

	next_block = bd[b_index]->bcp[EAST]->neighbour_block;
// 	cout << "East block index = " << next_block << endl;
	if (next_block == -1) status = FAILURE;
	else b_index = static_cast<size_t>(next_block);
    }
    
    if (status == FAILURE) {
	cout << "Error in initialise_piston_faces().\n";
	cout << "Position exceeds grid.\n";
	cout << pos << endl;
	cout << "Check piston geometry, blocks or block connections.\n";
	
	return FAILURE;
    }

    return SUCCESS;
}

int locate_block_index(size_t &b_near, Vector3 &pos, Block *(bd[]), size_t nblock)
{
    // Searches for pos with respect to domain boundaries.
    // Returns one of three cases:
    //     Case A: point lies between domain boundaries
    //     Case B: point lies at domain boundary
    //     Case C: point lies outside domain boundary

    if (DB0) cout << "entering locate_block_index()" << endl;
#if 0
    size_t j;
    b_near = -1;  // unfound...

    // Case A: Point lies between domain boundaries
    for (j = 0; j < nblock; ++j) {

	if ((pos.x > bd[j]->get_vtx(bd[j]->imin,bd[j]->jmin)->pos[0].x) &&
	    (pos.x < bd[j]->get_vtx(bd[j]->imax+1,bd[j]->jmin)->pos[0].x)) {
	    
	    if ((pos.y >= bd[j]->get_vtx(bd[j]->imin,bd[j]->jmin)->pos[0].y) && 
		(pos.y < bd[j]->get_vtx(bd[j]->imin,bd[j]->jmax+1)->pos[0].y)) {
		
		b_near = j;
	    }
	}
    }

    // Case B: Point lies at domain boundary
    for (j = 0; j < nblock; ++j) {
	if ((fabs(pos.x - bd[j]->get_vtx(bd[j]->imin,bd[j]->jmin)->pos[0].x) < WIDTH_TOL) || 
	    (fabs(pos.x - bd[j]->get_vtx(bd[j]->imax+1,bd[j]->jmin)->pos[0].x) < WIDTH_TOL)) {

	    if (DB0) cout << "leaving locate_block_index()" << endl;
	    return END_OF_BLOCK;
	}
    }

    // Case C: Point lies outside the domain boundaries
    if (b_near == -1) {
	cout << pos << " lies outside the domain." << endl;
	b_near = 0;
	cout << "block[" << b_near << "] corner locations"<< endl;
	cout << "SW(" << bd[b_near]->get_vtx(bd[b_near]->imin,bd[b_near]->jmin)->pos[0].x << "," 
	     << bd[b_near]->get_vtx(bd[b_near]->imin,bd[b_near]->jmin)->pos[0].y << ")," 
	     << "SE(" << bd[b_near]->get_vtx(bd[b_near]->imax+1,bd[b_near]->jmin)->pos[0].x << "," 
	     << bd[b_near]->get_vtx(bd[b_near]->imax+1,bd[b_near]->jmin)->pos[0].y << ")\n"
	     << "NW(" << bd[b_near]->get_vtx(bd[b_near]->imin,bd[b_near]->jmax+1)->pos[0].x << "," 
	     << bd[b_near]->get_vtx(bd[b_near]->imin,bd[b_near]->jmax+1)->pos[0].y << "),"
	     << "NE(" << bd[b_near]->get_vtx(bd[b_near]->imax+1,bd[b_near]->jmax+1)->pos[0].x << "," 
	     << bd[b_near]->get_vtx(bd[b_near]->imax+1,bd[b_near]->jmax+1)->pos[0].y << ")\n";

	if (DB0) cout << "leaving locate_block_index()" << endl;
	return FAILURE;
    }
    else {
	if (DB0) cout << "leaving locate_block_index()" << endl;
	return SUCCESS;
    }
#else
    return SUCCESS;
#endif
}

static const double FLOP_TOL = 1.0e-12; // for testing floating-point equality

/// \brief Returns the index of first interface within piston.
///
/// This function searches the grid from the west to find the
/// the first interface which is "within" the piston.
/// Using this information we can locate the appropriate cells.
///
/// \author Brendan T. O'Flaherty
///
int linear_search_from_west(size_t &ival, double xval, Block &bd)
{
    // Searches a list of grid interface locations
    // Returns the first interface location AFTER xval

    size_t i, j;
    double x_test;
    ival = bd.imax+1;

    j = bd.jmin;
    for (i = bd.imin; i < bd.imax+1; ++i) {
	x_test = bd.get_ifi(i,j)->pos.x;
	if (xval < x_test || fabs(xval - x_test) <= FLOP_TOL) { 
	    ival = i;
	    break;
	}
    }

    return SUCCESS;
}

int linear_search_from_east(size_t &ival, double xval, Block &bd)
{
    // Searches a list of grid interface locations
    // Returns the first interface location AFTER xval

    size_t i, j;
    double x_test;
    ival = bd.imin;

    j = bd.jmin;
    for (i = bd.imax+1; i >= bd.imin; --i) {
	x_test = bd.get_ifi(i,j)->pos.x;
	if (xval > x_test || fabs(xval - x_test) <= FLOP_TOL) { 
	    ival = i;
	    break;
	}
    }

    return SUCCESS;
}

// ----------------------------------------------------------------------------
// Additional functions

int linear_eval(double xval, double &yval, vector<double> &x, vector<double> &y)
{
    int i0, i1;
    double wA, wB;

    // ASSERT(x.size() == y.size());
    
    // trivial cases
    if (x.size() == 1) {
	yval = y[0];
	return 0;
    }

    // Implementation of a 1-D lookup-table
    // with linear interpolation

    // Search for appropriate index (lower-bound)
    i0 = 0;

    for ( ; i0 < static_cast<int>(x.size()); ++i0) {
	if ( x[i0] >= xval ) {
	    i0--;
	    break;
	}
    }
    
    // Establish weight values

    // 1. Handle edges of table
    if( i0 == -1 ) {
	i0 = 0;
	i1 = i0 + 1;
	wA = 1.0;
	wB = 0.0;
    }
    else if( i0 == int(x.size())-1 ) {
	i0 = int(x.size()-2);
	i1 = i0 + 1;
	wA = 0.0;
	wB = 1.0;
    }
    // 2. and all other points (interior)
    else {
	i1 = i0 + 1;
	wA = (x[i1] - xval)/(x[i1] - x[i0]);
	wB = (xval - x[i0])/(x[i1] - x[i0]);
    }

    // Perform the interpolation once the
    // weights are established.
    yval = wA*y[i0] + wB*y[i1];
    
    return SUCCESS;
}

int determine_time_step_size(global_data &G, Piston &piston, double dt_global, double cfl_target)
{
    double dt_local, cfl_local, piston_signal;

    piston_signal = piston.compute_signal_frequency();
    cfl_local = dt_global * piston_signal;
    dt_local = cfl_target / piston_signal;
    
    if (cfl_local < G.cfl_min) G.cfl_min = cfl_local;
    if (cfl_local > G.cfl_max) G.cfl_max = cfl_local;
    if (dt_local < G.dt_allow) G.dt_allow = dt_local;
    
    if( G.cfl_min > 0.0 && G.cfl_max < 0.9 ) {
	return SUCCESS;
    }
    else {
	cout << "Determine timestep size: bad CFL number encountered.\n";
	cout << "cfl_max= " << G.cfl_max << " for piston: " << piston.get_id() << endl;
	return DT_SEARCH_FAILED;
    }

}

// ----------------------------------------------------------------------------

int set_cell_status(FV_Cell *cell, int status) 
{
    cell->status = status;
    return SUCCESS;
}

int print_flow_state(FV_Cell *cell)
{
    double f_sum;
    int nsp = get_gas_model_ptr()->get_number_of_species();

    cout << "--------------------------------"<< endl;
    cout << "---------- cell data -----------" << endl;

    cout << "cell status = " << cell->status << endl;

    cout << "pos = (" 
	 << cell->pos[0].x << ", " 
	 << cell->pos[0].y << ", " 
	 << cell->pos[0].z << ")" << endl;

    cout << "Primary variables:\n"
	 << "rho = " << cell->fs->gas->rho << "\n"
         << "p   = " << cell->fs->gas->p << "\n"
         << "T[0]= " << cell->fs->gas->T[0] << "\n"
	 << "u   = " << cell->fs->vel.x << "\n"
	 << "v   = " << cell->fs->vel.y << "\n"
	 << "w   = " << cell->fs->vel.z << "\n"
	 << "e[0]= " << cell->fs->gas->e[0] << endl;
    cout << endl;

    cout << "Conserved quantities:\n"
	 << "rho  = " << cell->U[1]->mass << "\n"
	 << "rv.x = " << cell->U[1]->momentum.x << "\n"
	 << "rv.y = " << cell->U[1]->momentum.y << "\n"
	 << "rv.z = " << cell->U[1]->momentum.z << "\n"
	 << "rE   = " << cell->U[1]->total_energy << endl;
    cout << endl;

    f_sum = 0.0;
    for ( int isp = 0; isp < nsp; ++isp ) {
	cout << "f[" << isp << "]  = " << cell->fs->gas->massf[isp] << endl;
	f_sum += cell->fs->gas->massf[isp];
    }
    cout << "f_sum = " << f_sum << endl;
    cout << "base_qdot = " << cell->base_qdot << endl;
    cout << endl;

    cout << "--------------------------------"<< endl;

    return SUCCESS;
}

int print_gas_fluxes(FV_Cell *cell)
{
    cout << "--------------------------------"<< endl;
    cout << "---------- cell data -----------" << endl;
    cout << " DrDt    = " << cell->dUdt[0]->mass << endl;
    cout << " DrvDt.x = " << cell->dUdt[0]->momentum.x;
    cout << " DrvDt.y = " << cell->dUdt[0]->momentum.y;
    cout << " DrvDt.z = " << cell->dUdt[0]->momentum.z;
    cout << " DrEDt   = " << cell->dUdt[0]->total_energy;
    cout << "--------------------------------"<< endl;

    return SUCCESS;
}

int update_cell_geometry_based_on_vertices(FV_Cell *cell)
{
    global_data &G = *get_global_data_ptr();
    Vector3 A, B, C, D;
    double xyarea;

    // Cell vertices
    D = cell->vtx[0]->pos[0];
    A = cell->vtx[1]->pos[0];
    B = cell->vtx[2]->pos[0];
    C = cell->vtx[3]->pos[0];

    // Cell area in the (x,y)-plane.
    xyarea = 0.5 * ((B.x + A.x) * (B.y - A.y) + (C.x + B.x) * (C.y - B.y) +
		    (D.x + C.x) * (D.y - C.y) + (A.x + D.x) * (A.y - D.y));
    
    cell->area[0] = xyarea;

    // Cell i,j lengths
    cell->iLength = A.x - D.x;
    cell->jLength = C.y - D.y;

    // Cell centroid
    cell->pos[0].x = 1.0 / (xyarea * 6.0) * 
	((B.y - A.y) * (A.x * A.x + A.x * B.x + B.x * B.x) + 
	 (C.y - B.y) * (B.x * B.x + B.x * C.x + C.x * C.x) +
	 (D.y - C.y) * (C.x * C.x + C.x * D.x + D.x * D.x) + 
	 (A.y - D.y) * (D.x * D.x + D.x * A.x + A.x * A.x));
    cell->pos[0].y = -1.0 / (xyarea * 6.0) * 
	((B.x - A.x) * (A.y * A.y + A.y * B.y + B.y * B.y) + 
	 (C.x - B.x) * (B.y * B.y + B.y * C.y + C.y * C.y) +
	 (D.x - C.x) * (C.y * C.y + C.y * D.y + D.y * D.y) + 
	 (A.x - D.x) * (D.y * D.y + D.y * A.y + A.y * A.y));

    if ( G.axisymmetric ) {
	cell->volume[0] = xyarea * cell->pos[0].y; // volume per unit radian
    } else {
	cell->volume[0] = xyarea * 1.0; // volume per unit depth
    }

    if (DB5) { 
	cout << "vertices: \n"
	     << "(" << D.x << ", " << D.y << ")\n"
	     << "(" << A.x << ", " << A.y << ")\n"
	     << "(" << B.x << ", " << B.y << ")\n"
	     << "(" << C.x << ", " << C.y << ")"
	     << endl;
	cout << "cell volume = " << cell->volume[0] << endl;
	cout << "cell area = " << cell->area[0] << endl;
    }

    // ------------------------------------------------------------------------
    // Cell interfaces

    update_interface_geometry_based_on_vertices(cell->iface[SOUTH], D, A);
    update_interface_geometry_based_on_vertices(cell->iface[EAST], B, A);
    update_interface_geometry_based_on_vertices(cell->iface[NORTH], C, B);
    update_interface_geometry_based_on_vertices(cell->iface[WEST], C, D);
    
    // ------------------------------------------------------------------------
    // error handling

    if ((cell->iLength) < 0.0) {
	cout << "Negative cell iLength = " << cell->iLength << endl;
	return FAILURE;
    }	
    if ((cell->jLength) < 0.0) {
	cout << "Negative cell jLength = " << cell->jLength << endl;
	return FAILURE;
    }
    if (xyarea < 0.0) {
	cout << "Negative flex_cell volume = " << xyarea << endl;
	cout << "pos = " 
	     << "(" << cell->pos[0].x << ", " << cell->pos[0].y << ")"
	     << endl;
	cout << "vertices: \n"
	     << "(" << D.x << ", " << D.y << ")\n"
	     << "(" << A.x << ", " << A.y << ")\n"
	     << "(" << B.x << ", " << B.y << ")\n"
	     << "(" << C.x << ", " << C.y << ")"
	     << endl;
	return FAILURE;
    }
    
    return SUCCESS;
}

int update_interface_geometry_based_on_vertices(FV_Interface *iface, Vector3 vtxL, Vector3 vtxR)
{
    global_data &G = *get_global_data_ptr();
    if (DB0) cout << "entering update_interface_geometry_based_on_vertices()" << endl;
    double L;
    
    L = sqrt((vtxR.x - vtxL.x)*(vtxR.x - vtxL.x) + (vtxR.y - vtxL.y)*(vtxR.y - vtxL.y));
    iface->n.x = (vtxR.y - vtxL.y)/L;
    iface->n.y = -(vtxR.x - vtxL.x)/L;
    iface->n.z = 0.0;
    
    iface->length = L;
    iface->Ybar = 0.5*(vtxR.y + vtxL.y);
    if ( G.axisymmetric ) {
	// Area per unit radian
	iface->area[0] = L * iface->Ybar;
    } else {
	// Area[0] per unit depth
        iface->area[0] = L;
    }
    
    iface->pos = 0.5*(vtxR + vtxL);
    
    if (L < 1.0e-9) {
	cout << "Zero interface length " << L << endl;
	return FAILURE;
    }
    
    if (DB0) cout << "leaving update_interface_geometry_based_on_vertices()" << endl;
    return SUCCESS;
}

int update_gas_extensive_conserved(FV_Cell *cell, double volume)
{
    // after this call all conserved quantites will be extensive values 
    // note that density (incl. that of individual species) is converted to mass

    if (DB0) cout << "entering update_gas_extensive_conserved()" << endl;

    size_t nsp = get_gas_model_ptr()->get_number_of_species();
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();

    // ASSERT(volume > 0.0);
    update_gas_conserved(cell);

    cell->U[1]->mass *= volume;  
    cell->U[1]->momentum *= volume;
    cell->U[1]->total_energy *= volume;
    for ( size_t isp = 0; isp < nsp; ++isp) {
	cell->U[1]->massf[isp] *= volume;             // mass of species
    }
    for ( size_t imode = 0; imode < nmodes; ++imode) {
	cell->U[1]->energies[imode] *= volume;         // independent energy  
    }
    if (DB0) cout << "leaving update_gas_extensive_conserved()" << endl;

    return SUCCESS;
}

int update_gas_conserved(FV_Cell *cell) 
{
    cell->encode_conserved(0, 0, 0.0, false);
    return SUCCESS;
}

int update_gas_primaries_from_extensive(FV_Cell *cell, double volume)
{
    if (DB0) cout << "entering update_gas_primaries_from_extensive()" << endl;
    
    // ASSERT(volume > 0.0);
    size_t nsp = get_gas_model_ptr()->get_number_of_species();
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();
    
    cell->U[1]->mass /= volume;  
    cell->U[1]->momentum /= volume;
    cell->U[1]->total_energy /= volume;
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	cell->U[1]->massf[isp] /= volume;             // mass of species
    }
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	cell->U[1]->energies[imode] /= volume;         // independent energies  
    }
    update_gas_primaries(cell);

    // Leave conserved quantities as we found them.

    // n.b. rho is both a primary and conserved property. it is
    // more important that conserved data be correct at this point.
    // this is a fundamental difference to casbar, which stores 
    // conserved quantities in a separate vector.
    // remember to convert to density before copying back to the grid.

    cell->U[1]->mass *= volume; 
    cell->U[1]->momentum *= volume;
    cell->U[1]->total_energy *= volume;
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	cell->U[1]->massf[isp] *= volume;             // mass of species
    }
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	cell->U[1]->energies[imode] *= volume;         // independent energies  
    }
    if (DB0) cout << "leaving update_gas_primaries_from_extensive()" << endl;

    return SUCCESS;
}

int update_gas_primaries(FV_Cell *cell)
{
    cell->decode_conserved(0, 0, 0.0, false);
    return SUCCESS;
}

int update_extensive_gas_time_derivatives_from_fluxes(FV_Cell *cell)
{
    FV_Interface *IFn, *IFs, *IFe, *IFw;
    double integral;

    size_t nsp = get_gas_model_ptr()->get_number_of_species();

    IFn = cell->iface[NORTH];
    IFe = cell->iface[EAST];
    IFs = cell->iface[SOUTH];
    IFw = cell->iface[WEST];

    // Time-derivative for Mass.
    integral = 
	- IFe->F->mass * IFe->area[0] 
	- IFn->F->mass * IFn->area[0] 
	+ IFw->F->mass * IFw->area[0] 
	+ IFs->F->mass * IFs->area[0];
    cell->dUdt[0]->mass = integral + cell->Q->mass * cell->volume[0];

    // Time-derivative for X-Momentum.
    integral = 
	- IFe->F->momentum.x * IFe->area[0] 
	- IFn->F->momentum.x * IFn->area[0]
	+ IFw->F->momentum.x * IFw->area[0] 
	+ IFs->F->momentum.x * IFs->area[0];
    cell->dUdt[0]->momentum.x = integral + cell->Q->momentum.x * cell->volume[0];

    // Time-derivative for Y-Momentum.
    integral = 
	- IFe->F->momentum.y * IFe->area[0] 
	- IFn->F->momentum.y * IFn->area[0]
	+ IFw->F->momentum.y * IFw->area[0] 
	+ IFs->F->momentum.y * IFs->area[0];
    cell->dUdt[0]->momentum.y = integral + cell->Q->momentum.y * cell->volume[0];
    
    // Time-derivative for Z-Momentum.
    cell->dUdt[0]->momentum.z = 0.0;
    
    // Time-derivative for Total Energy.
    integral = 
	- IFe->F->total_energy * IFe->area[0] 
	- IFn->F->total_energy * IFn->area[0]
	+ IFw->F->total_energy * IFw->area[0] 
	+ IFs->F->total_energy * IFs->area[0];
    cell->dUdt[0]->total_energy = integral + cell->Q->total_energy * cell->volume[0];

    // Time-derivative for individual species.
    // Units of DfDt are 1/sec.
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	integral =
	    - IFe->F->massf[isp] * IFe->area[0]
	    - IFn->F->massf[isp] * IFn->area[0]
	    + IFw->F->massf[isp] * IFw->area[0]
	    + IFs->F->massf[isp] * IFs->area[0];
	cell->dUdt[0]->massf[isp] = integral;
    }

    return SUCCESS;
}

int update_gas_conserved_from_time_derivatives(FV_Cell *cell, double dt)
{
    size_t nsp = get_gas_model_ptr()->get_number_of_species();
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();
    cell->U[1]->mass = cell->U[0]->mass + dt * cell->dUdt[0]->mass;
    cell->U[1]->momentum.x = cell->U[0]->momentum.x + dt * cell->dUdt[0]->momentum.x;
    cell->U[1]->momentum.y = cell->U[0]->momentum.y + dt * cell->dUdt[0]->momentum.y;
    cell->U[1]->momentum.z = cell->U[0]->momentum.z + dt * cell->dUdt[0]->momentum.z;
    cell->U[1]->total_energy = cell->U[0]->total_energy + dt * cell->dUdt[0]->total_energy;
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	cell->U[1]->massf[isp] = cell->U[0]->massf[isp] + dt * cell->dUdt[0]->massf[isp];
    }
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	cell->U[1]->energies[imode] = cell->U[0]->energies[imode] + dt * cell->dUdt[0]->energies[imode];
    }
    return SUCCESS;
}

int average_two_flow_states(FV_Cell *dst,
			    FV_Cell *src0, double vol0,
			    FV_Cell *src1, double vol1)
{ 
    double total_vol = vol0 + vol1;
    size_t nsp = get_gas_model_ptr()->get_number_of_species();
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();
    dst->U[1]->mass = (src0->U[1]->mass*vol0 + src1->U[1]->mass*vol1)/total_vol;
    dst->U[1]->momentum.x = (src0->U[1]->momentum.x*vol0 + src1->U[1]->momentum.x*vol1)/total_vol;
    dst->U[1]->momentum.y = (src0->U[1]->momentum.y*vol0 + src1->U[1]->momentum.y*vol1)/total_vol;
    dst->U[1]->momentum.z = (src0->U[1]->momentum.z*vol0 + src1->U[1]->momentum.z*vol1)/total_vol;
    dst->U[1]->total_energy = (src0->U[1]->total_energy*vol0 + src1->U[1]->total_energy*vol1)/total_vol;
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	dst->U[1]->massf[isp] = (src0->U[1]->massf[isp]*vol0 + src1->U[1]->massf[isp]*vol1)/total_vol;
    }
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	dst->U[1]->energies[imode] = (src0->U[1]->energies[imode]*vol0 + src1->U[1]->energies[imode]*vol1)/total_vol;
    }
    dst->decode_conserved(0, 0, 0.0, false);
    return SUCCESS;
}

