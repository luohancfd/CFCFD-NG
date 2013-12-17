/// \file piston.hh
/// \ingroup mbcns2
/// \brief Class for handling pistons.
/// \author Brendan O'Flaherty and Rowan J. Gollan
///
/// \version 16-Mar-2007 -- Began copying from mbcns2 piston
/// \version 31-Mar-2007 -- Both codes consistant, friction added along with two test cases. 
/// \version 02-Apr-2007 -- Piston faces now composed of multiple component-faces
/// \version 07-Jan-2008 -- Ported from two_phase code
/// 
/// ---------------------------------------------------------

#ifndef PISTON_HH
#define PISTON_HH

// #include <stdio.h>
// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <iomanip>

// #include <string>
// #include <vector>

#include <string>
#include <vector>
#include <fstream>

#include "../../../lib/util/source/config_parser.hh"
#include "../../../lib/util/source/dbc_assert.hh"
#include "../../../lib/geometry2/source/geom.hh"
#include "cell.hh"
#include "block.hh"
struct global_data;

#define SIGN(a) (((a) == 0.0) ? 0.0 : (a)/(fabs(a)))

// Piston exit codes
const int END_OF_BLOCK = 2;

enum PistonFaceType { EAST_FACE, WEST_FACE, UNSPECIFIED };

// ----------------------------------------------------------------------------
// \brief PistonFace declaration
//

class PistonFace {
public:
    /// \brief Default constructor
    PistonFace();

    /// \brief Copy constructor
    PistonFace(const PistonFace &pf);

    /// \brief Default destructor
    ~PistonFace();
    PistonFace& operator=(const PistonFace &pf);
    PistonFace* clone() { return new PistonFace(*this); }
    
    // ------------------------------------------------------------------------
    // service functions

    size_t get_io() { return io_; }
    size_t get_bo() { return bo_; }
    size_t get_ii() { return ii_; }
    size_t get_bi() { return bi_; }

    size_t get_jmin() { return jmin_; }
    size_t get_jmax() { return jmax_; }

    void set_jmin(size_t jmin) { jmin_ = jmin; }
    void set_jmax(size_t jmax) { jmax_ = jmax; }

    Vector3* get_pos() { return &vtx_C_; }
    void set_pos(Vector3 pos) { vtx_C_ = pos; }
    void move_face(Vector3 &pos);
    
    double get_velocity_x() { return vel_->x; }

    double get_dF() { return dFg_; }
    double get_dFg() { return dFg_; }
    double get_dR() { return dR_; }

    bool covers_two_cells()
    { return covers_two_cells_; }

    /// \brief Return a pointer to the jth flex cell
    FV_Cell* get_flex_cell(size_t j)
    { 
	// ASSERT(j < nni_*njcells_);
	return flex_cells_[j]; 
    }
    
    /// \brief Return a pointer to the ifi interface at location (i, j).
    FV_Interface* get_ifi(size_t i, size_t j)
    { 
	// ASSERT((j*nii_+i) < (nvi_*nvj_));
	return ifi_[j*nii_+i]; 
    }

    /// \brief Return a pointer to the ifj interface at location (j).
    FV_Interface* get_ifj(size_t j)
    { 
	// ASSERT(j < (nii_*nij_));
	return ifj_[j]; 
    }

    /// \brief Return a pointer to the vertex at location (i, j).
    FV_Vertex* get_vertex(size_t i, size_t j)
    { 
	// ASSERT((j*nvi_+i) < (nii_*nij_));
	return vtx_[j*nvi_+i];
    }

    // ------------------------------------------------------------------------
    // member functions

    int set_values(double n_x,
		   Vector3 vtxL,
		   Vector3 vtxR,
		   Vector3* vel,
		   PistonFaceType type);

    int set_indices_for_flex_cells(Block &bd);
    int locate_shadow_indices(global_data &G,
			      Block *(bd[]));
    int allocate_memory();
    int bind_interfaces_to_flex_cells();

    int update_west_flex_cell_after_motion();
    int update_east_flex_cell_after_motion();
    int update_flex_cell_states(double dt);

    int compute_forces(std::vector<double> bore_resistance_f,
		       std::vector<double> bore_resistance_x);

    int copy_geometry_data_to_flex_cells(Block *(bd[]));
    int copy_flow_data_to_flex_cells(Block *(bd[]));
    int copy_flow_data_to_interfaces(Block *(bd[]));
    int copy_flow_data_to_grid(Block *(bd[]));

    int shadow_cells(Block *(bd[]));
    int apply_boundary_condition();

private:

    size_t id_;   ///< \brief Id for pistonface.

    size_t io_;   ///< \brief index for "outer" shadow cell
    size_t bo_;   ///< \brief block index in which io_ is located
    size_t ii_;   ///< \brief index for "inner" shadow cell
    size_t bi_;   ///< \brief block index in which ii_ is located

    bool covers_two_cells_; ///< \brief flag indicating if flex_cell covers two cells

    size_t jmin_; ///< \brief minimum j index of flex_cells
    size_t jmax_; ///< \brief maximum j index of flex_cells

    size_t nni_; 
    size_t nnj_;
    size_t nvi_;
    size_t nvj_;
    size_t nii_;
    size_t nij_;
    size_t njcells_;
    
    double dFg_;          ///< \brief gas component force
    double dR_;

    Vector3 n_;           ///< \brief unit normal
    Vector3 vtx_C_;       ///< \brief face center Vector3
    Vector3 vtx_L_;       ///< \brief face left Vector3
    Vector3 vtx_R_;       ///< \brief face right Vector3

    Vector3 *vel_;        ///< \brief piston velocity

    PistonFaceType type_; ///< \brief piston face type

    std::vector<FV_Cell*> flex_cells_;
    std::vector<FV_Interface*> ifi_;
    std::vector<FV_Interface*> ifj_;
    std::vector<FV_Vertex*> vtx_;
};

// ----------------------------------------------------------------------------
// \brief Piston declaration
//

enum PistonStatus { P_ACTIVE,
                    P_INACTIVE
};

class Piston {
public:
    /// \brief Default constructor
    Piston();
    Piston(const Piston &p);

    /// \brief Default destructor
    ~Piston();
    Piston* clone() { return new Piston(*this); }
    Piston& operator=(const Piston &p);
    
    // ------------------------------------------------------------------------
    // service functions

    size_t get_id() { return id_; }
    std::string string_repr();
    std::string vtk_string();
    int write_state(FILE *fp, double t);
    int read_state(FILE *fp);
    int save_state();
    int restore_state();
    
    void set_const_v_flag(bool const_v_flag) { const_velocity_ = const_v_flag; }
    void set_postv_v_flag(bool postv_v_flag) { postv_velocity_ = postv_v_flag; }
    
    /// \brief Return a pointer to NW vertex
    FV_Vertex* get_NW_vertex() 
	{ return wPistonFace.back()->get_vertex(0,wPistonFace.back()->get_jmax()+1); }
    
    /// \brief Return a pointer to NE vertex
    FV_Vertex* get_NE_vertex() 
	{ return ePistonFace.back()->get_vertex(1,ePistonFace.back()->get_jmax()+1); }
    
    /// \brief Return a pointer to SE vertex
    FV_Vertex* get_SE_vertex() 
	{ return ePistonFace.front()->get_vertex(1,ePistonFace.front()->get_jmin()); }
    
    /// \brief Return a pointer to SW vertex
    FV_Vertex* get_SW_vertex() 
	{ return wPistonFace.front()->get_vertex(0,wPistonFace.front()->get_jmin()); }
    
    // ------------------------------------------------------------------------
    // member functions

    int set_values(size_t id,
		   double D,
		   double L,
		   double m,
		   double x,
		   double u,
		   double rifling_twist,
		   double rog,
		   double vanish_at_x,
		   std::vector<double> bore_resistance_f,
		   std::vector<double> bore_resistance_x);

    int initialise_from_config_file(global_data &G,
				    Block *(bd[]),
				    ConfigParser &cfg, 
				    string section);
    int initialise_piston_faces(global_data &G,
				Block *(bd[]));
    
    int change_of_state_due_to_motion(global_data &G, 
				      Block *(bd[]), double dt);
    int print_flex_cell_flow_states();
    
    int update_piston_state(global_data &G, 
			    Block *(bd[]), double dt);
    int configure_cells(global_data &G, 
			Block *(bd[]));

    int deactivate_cells(Block *(bd[]));
    
    int update_flex_cell_states(double dt);
    int update_piston_geometry_based_on_vertices();
    int update_piston_geometry_based_on_position(); 

    int copy_data_to_piston(Block *(bd[]));
    int copy_flow_data_to_grid(Block *(bd[]));

    /// \brief Return the signal frequency (related to the CFL check).
    double compute_signal_frequency();

    size_t get_wio()
    { return wPistonFace[0]->get_io(); }

    size_t get_wbo()
    { return wPistonFace[0]->get_bo(); }

    size_t get_wii()
    { return wPistonFace[0]->get_io(); }

    size_t get_wbi()
    { return wPistonFace[0]->get_bo(); }

    size_t get_eio()
    { return ePistonFace[0]->get_ii(); }

    size_t get_ebo()
    { return ePistonFace[0]->get_bi(); }

    size_t get_eii()
    { return ePistonFace[0]->get_ii(); }

    size_t get_ebi()
    { return ePistonFace[0]->get_bi(); }

    /// \brief Return a boolean about whether or not a piston is active
    bool is_active()
    { return (status_ == P_ACTIVE); }

    /// \brief Return a boolean about whether or not a piston is inactive
    bool is_inactive()
    { return (status_ == P_INACTIVE); }

private:
    size_t id_;        // piston id
    size_t ncf_;       // number of component faces per face
    size_t imin_;
    size_t imax_;

    PistonStatus status_; 
    double vanish_at_x_;

    bool const_velocity_;   ///< \brief flag that forces constant velocity
    bool postv_velocity_;   ///< \brief flag that allows positive velocity only

    double dx_;     // displacement (m) 
    double diam_;   // diameter (m) 
    double area_;   // discretised area (m**2)
    double m_;      // mass (kg)
    double acc_;    // acceleration (m/s**2)
    double I_;      // rotational inertia (kg.m**2)
    double rog_;    // radius of gyration (m)
    double beta_;   // rifling angle (rads)
    double n_;      // number of turns per calibre

    Vector3 pos_;   // position of center (m)
    Vector3 vel_;   // velocity (m/s)

    Vector3 pos0_;
    Vector3 vel0_;

    std::vector<int> bmin_;
    std::vector<int> bmax_;

    std::vector<double> bore_resistance_f_;
    std::vector<double> bore_resistance_x_;
    
    std::vector<Vector3> vtx_;
    std::vector<PistonFace*> wPistonFace;
    std::vector<PistonFace*> ePistonFace;
};

// ----------------------------------------------------------------------------
// Search functions
//
int activate_cells(global_data &G, Block *(bd[]));
int print_cell_status(global_data &G, Block *(bd[]));
int locate_west_block_index(size_t &b_index, Vector3 &pos, global_data &G, Block *(bd[]));
int locate_east_block_index(size_t &b_index, Vector3 &pos, global_data &G, Block *(bd[]));
int locate_block_index(size_t &b_near, Vector3 &pos, Block *(bd[]), size_t nblock);

///\brief returns index of "nearest" interface
int linear_search_from_west(size_t &ival, double xval, Block &bd);
int linear_search_from_east(size_t &ival, double xval, Block &bd);

///\brief returns the linear interpolation of xval in y
int linear_eval(double xval, double &yval, vector<double> &x, vector<double> &y);
int determine_time_step_size(global_data &G, Piston &piston, 
			     double dt_global, double cfl_target);

// ----------------------------------------------------------------------------
// Additional functions
//

int set_cell_status(FV_Cell *cell, int status);
int print_flow_state(FV_Cell *cell);
int print_gas_fluxes(FV_Cell *cell);
int update_cell_geometry_based_on_vertices(FV_Cell *cell);
int update_interface_geometry_based_on_vertices(FV_Interface *iface, Vector3 vtxL, Vector3 vtxR);
int update_gas_conserved(FV_Cell *cell);
int update_gas_extensive_conserved(FV_Cell *cell, double volume);
int update_gas_primaries(FV_Cell *cell);
int update_gas_primaries_from_extensive(FV_Cell *cell, double volume);
int update_extensive_gas_time_derivatives_from_fluxes(FV_Cell *cell);
int update_gas_conserved_from_time_derivatives(FV_Cell *cell, double dt);
int average_two_flow_states(FV_Cell *dst,
			    FV_Cell *src0, double vol0,
			    FV_Cell *src1, double vol1);

#endif // PISTON_HH
