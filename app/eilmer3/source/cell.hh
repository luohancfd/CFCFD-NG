/** \file cell.hh
 * \ingroup eilmer3
 * \brief Header file for the cell centre, interface and vertex data.
 * 
 * \author PJ
 * \version 05-Aug-04 extracted from mb_cns.h
 */

#ifndef CELL_HH
#define CELL_HH

#include <string>
#include <vector>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/geometry2/source/geom.hh"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"
#include "c-flow-condition.hh"

/*--------------------------------------------------------------*/

/** These values are used by the cell data copying functions
 * to specify which parts of the cell data to copy.
 */
const int COPY_ALL_CELL_DATA = 0;
const int COPY_FLOW_STATE = 1;
const int COPY_CELL_LENGTHS = 2;
const int COPY_INTERFACE_DATA = 3;

/** Symbolic labels and indices for the cell's faces. 
 *
 * The names of the faces of the structured-grid blocks will be the same.
 */
const int NORTH = 0;
const int EAST = 1;
const int SOUTH = 2;
const int WEST = 3;
const int TOP = 4;
const int BOTTOM = 5;

/// States that a cell-center may be in with respect to pistons.
const int NORMAL_CELL = 0;
const int MASKED_CELL = 1;
const int SHADOWED_CELL = 2;
/// States that cell-interfaces may be in with respect to pistons.
const int NORMAL_IFACE = 0;
const int MASKED_IFACE = 1;

// Update scheme
enum update_scheme_t { EULER_UPDATE, PC_UPDATE, MIDPOINT_UPDATE,
		       CLASSIC_RK3_UPDATE, TVD_RK3_UPDATE, DENMAN_RK3_UPDATE };
update_scheme_t set_gasdynamic_update_scheme(update_scheme_t my_scheme);
update_scheme_t get_gasdynamic_update_scheme();
std::string get_name_of_gasdynamic_update_scheme(update_scheme_t my_scheme);
size_t number_of_stages_for_update_scheme(update_scheme_t my_scheme);
/// N_LEVEL is used to size the time-derivative vectors.
const size_t N_LEVEL = 4;

const size_t N_INTERFACE = 6; // Number of interfaces per cell
const size_t N_VERTEX = 8; // Number of vertices per cell

/// Minimum values for turbulent kinetic energy (m^2/s^2) and frequency (1/s)
/// for applying limiters in the k-omega model.
const double SMALL_TKE = 0.1;
const double SMALL_OMEGA = 1.0;

int get_face_index(const std::string name);
std::string  get_face_name(int face_index);

/** \brief Flow state describes the local flow but not the cell or interface geometry. */
class FlowState {
public:
    Gas_data *gas;       ///< \brief gas thermo state
    Vector3 vel;         ///< \brief gas velocity vector, m/s
    Vector3 B;           ///< \brief magnetic field, Tesla 
    int S;               ///< \brief flag to indicate shock-point
    double tke;          ///< \brief turbulence kinetic energy per unit mass
    double omega;        ///< \brief turbulence frequency or pseudo vorticity
    double mu_t;         ///< \brief turbulence viscosity
    double k_t;          ///< \brief turbulence conductivity
    std::vector<double> G; ///< \brief velocity dist. partial densities, kg/m**3
    std::vector<double> H; ///> \brief velocity dist. partial densities, (kg*s**2)/(m**5)
    //
    FlowState(Gas_model *gm);
    FlowState();
    FlowState(const FlowState &fs);
    FlowState & operator=(const FlowState &fs);
    ~FlowState();
    int print() const;
    int copy_values_from(const FlowState &src);
    int copy_values_from(const CFlowCondition &src);
    int average_values_from(const FlowState &src0, const FlowState &src1,
			    bool with_diff_coeff);
    double * copy_values_to_buffer(double *buf) const;
    double * copy_values_from_buffer(double *buf);
    int BGK_equilibrium(void);
};

/// \brief Conserved quantities "vector"
///
/// In a cell, the conserved quantities are per unit volume.
/// It is also used as the class for the flux-vector and the time-derivatives
/// of the conserved quantities.
class ConservedQuantities {
public:
    double mass;                  ///< \brief density, kg/m**3
    Vector3 momentum;             ///< \brief momentum/unit volume
    Vector3 B;                    ///< \brief magnetic field, Tesla
    double total_energy;          ///< \brief total energy
    std::vector<double> massf;    ///< \brief mass fractions of species
    std::vector<double> energies; ///< \brief modal energies (mode 0 is usually transrotational)
    double tke;                   ///< \brief turbulent kinetic energy
    double omega;                 ///< \brief omega from k-omega turbulence model
    std::vector<double> G; ///< \brief velocity dist. partial densities, kg/m**3
    std::vector<double> H; ///> \brief velocity dist. partial densities, (kg*s**2)/(m**5)
    //
    ConservedQuantities(Gas_model *gm);
    ConservedQuantities();
    ConservedQuantities(const ConservedQuantities &Q);
    ConservedQuantities & operator=(const ConservedQuantities &Q);
    ~ConservedQuantities();
    int print() const;
    int copy_values_from(const ConservedQuantities &src);
    int clear_values();
};

/// \brief Inviscid and viscous fluxes are transported (between cells) across cell interfaces.
class FV_Interface {
public:
    size_t id;  // allows us to work out where, in the block, the interface is
    int status; // state of the interface with respect to pistons
    // Geometry
    Vector3 pos;       ///< \brief position of the (approx) midpoint
    Vector3 vel;       ///< \brief gas velocity, m/s
    double Ybar;       ///< \brief Y-coordinate of the mid-point
    double length;     ///< \brief Interface length in the x,y-plane
    std::vector<double> area; ///< \brief  Area m**2 for each time-level.
                       ///< \brief Area per radian in axisymmetric geometry
    Vector3 n;         ///< \brief Direction cosines for unit normal
    Vector3 t1;        ///< \brief tangent vector 1 (aka p)
    Vector3 t2;        ///< \brief tangent vector 2 (aka q)
    // Flow
    FlowState *fs;     ///< \brief Flow properties
    ConservedQuantities *F; ///< \brief Flux conserved quantity per unit area
#   if WITH_IMPLICIT == 1
    // Point-implicit variables
    double Lambda[6][2], R[6][6], R_inv[6][6], sl[6][2], sl_min[6][2];
    double A[6][6], J[6][6], hl[6][2], gl[6][2];
#   endif
    //
    FV_Interface(Gas_model *gm);
    FV_Interface();
    FV_Interface(const FV_Interface &iface);
    FV_Interface & operator=(const FV_Interface &iface);
    ~FV_Interface();
    int print() const;
    int copy_values_from(const FV_Interface &src, int type_of_copy);
    int copy_grid_level_to_level(size_t from_level, size_t to_level);
}; // end of class FV_Interface


/// \brief Cell vertex holds data for the calculation of gradients for the viscous terms.
class FV_Vertex {
public:
    size_t id;  // allows us to work out where, in the block, the vertex is
    // Geometry
    std::vector<Vector3> pos;  ///< \brief x,y,z-Coordinates for time-levels, m
    std::vector<Vector3> vel;  ///< \brief velocity for time-levels, m/s
    double area;               ///< \brief x,y-plane area of secondary cells (for spatial derivatives)
    double volume;             ///< \brief volume of 3D secondary cells (for spatial derivatives)
    // Derivatives of primary-cell variables.
    double dudx, dudy, dudz;   ///< \brief velocity derivatives
    double dvdx, dvdy, dvdz;   ///< \brief velocity derivatives
    double dwdx, dwdy, dwdz;   ///< \brief velocity derivatives
    std::vector<double> dTdx, dTdy, dTdz; ///< \brief Temperature derivatives
    double dtkedx, dtkedy, dtkedz;        ///< \brief turbulence kinetic energy
    double domegadx, domegady, domegadz;  ///< \brief pseudo vorticity for k-omega turbulence
    std::vector<double> dfdx, dfdy, dfdz; ///< \brief mass fraction derivatives
    double dpedx, dpedy, dpedz;           ///< \brief electron pressure derivatives
    //
    FV_Vertex(Gas_model *gm);
    FV_Vertex();
    FV_Vertex(const FV_Vertex &vtx);
    FV_Vertex & operator=(const FV_Vertex &vtx);
    ~FV_Vertex();
    int copy_values_from(const FV_Vertex &src);
    int copy_grid_level_to_level(size_t from_level, size_t to_level);
};


/// \brief Cell-averaged data forms the core of the flow field data.
class FV_Cell {
public:
    size_t id;  ///> \brief allows us to work out where, in the block, the cell is
    int status; ///> \brief state of the cell with respect to pistons
    bool fr_reactions_allowed; ///> \brief if true, will call chemical_increment (also thermal_increment)
    double dt_chem; ///> \brief acceptable time step for finite-rate chemistry
    double dt_therm; ///> \brief acceptable time step for thermal relaxation
    bool in_turbulent_zone; ///> \brief if true, we will keep the turbulence viscosity
    double base_qdot; ///> \brief base-level of heat addition to cell, W/m**3
    // Geometry
    std::vector<Vector3> pos; ///> \brief Centre x,y,z-coordinates for time-levels, m,m,m
    std::vector<double> volume; ///> \brief Cell volume for time-levels (per unit depth or radian in 2D), m**3
    std::vector<double> area; ///> \brief (x,y)-plane area for time-levels, m**2
    double uf; ///> \brief uncovered fraction */
    double iLength; ///> \brief length in the i-index direction
    double jLength; ///> \brief length in the j-index direction
    double kLength; ///> \brief length in the k-index direction
    double L_min;   ///> \brief minimum length scale for cell
    double distance_to_nearest_wall; ///> \brief for turbulence model correction.
    double half_cell_width_at_wall;  ///> \brief ditto
    FV_Cell *cell_at_nearest_wall;   ///> \brief ditto
    // Connections
    std::vector<FV_Interface *> iface;  ///> \brief pointers to defining interfaces of cell
    std::vector<FV_Vertex *> vtx;  ///> \brief pointers to vertices for quad (2D) and hexahedral (3D) cells
    // Flow
    FlowState *fs; ///> \brief Flow properties
    std::vector<ConservedQuantities *> U;    ///> \brief Conserved flow quantities for the update stages.
    std::vector<ConservedQuantities *> dUdt; ///> \brief Time derivatives for the update stages.
    ConservedQuantities *Q; ///> \brief source (or production) terms
    // Terms for loose-coupling of radiation.
    double Q_rad_org;
    double f_rad_org;
    double Q_rE_rad; ///> \brief Rate of energy addition to cell via radiation.
    double Q_rE_rad_save; // Presently, the radiation source term is calculated
                          // at the first update stage.  We need to retain that
                          // value for all of the update stages.
    // Data for computing residuals.
    double rho_at_start_of_step, rE_at_start_of_step;
#   if WITH_IMPLICIT == 1
    /* POINT-IMPLICIT VARIABLES */
    double piM[6][6], pir[6][2], qL[6][2], g_Ll[6][2], gjvec[6][2], gjmtx[6][6];
#   endif
    //
    FV_Cell(Gas_model *gm);
    FV_Cell();
    FV_Cell(const FV_Cell &cell);
    FV_Cell & operator=(const FV_Cell &cell);
    ~FV_Cell();
    int print() const;
    int point_is_inside(Vector3 &p, int dimensions, size_t gtl) const;
    int copy_values_from(const CFlowCondition &src);
    int copy_values_from(const FV_Cell &src, int type_of_copy, size_t gtl);
    double * copy_values_to_buffer(double *buf, int type_of_copy, size_t gtl) const;
    double * copy_values_from_buffer(double *buf, int type_of_copy, size_t gtl);
    int copy_grid_level_to_level(size_t from_level, size_t to_level);
    int replace_flow_data_with_average(std::vector<FV_Cell *> src);
    int scan_values_from_string(char *bufptr);
    std::string write_values_to_string() const;
    int scan_BGK_from_string(char *bufptr);
    std::string write_BGK_to_string() const;
    int impose_chemistry_timestep(double dt);
    int impose_thermal_timestep(double dt);
    int set_fr_reactions_allowed(int flag);
    int encode_conserved(size_t gtl, size_t ftl, double omegaz, bool with_k_omega);
    int decode_conserved(size_t gtl, size_t ftl, double omegaz, bool with_k_omega);
    bool check_flow_data(void);
    int time_derivatives(size_t gtl, size_t ftl, size_t dimensions, bool with_k_omega);
    int stage_1_update_for_flow_on_fixed_grid(double dt, bool force_euler, bool with_k_omega);
    int stage_2_update_for_flow_on_fixed_grid(double dt, bool with_k_omega);
    int stage_3_update_for_flow_on_fixed_grid(double dt, bool with_k_omega);
    int stage_1_update_for_flow_on_moving_grid(double dt, bool with_k_omega);
    int stage_2_update_for_flow_on_moving_grid(double dt, bool with_k_omega);
    int chemical_increment(double dt, double T_frozen);
    int thermal_increment(double dt, double T_frozen_energy);
    double signal_frequency(size_t dimensions, bool with_k_omega);
    int turbulence_viscosity_zero();
    int turbulence_viscosity_zero_if_not_in_zone();
    int turbulence_viscosity_limit(double factor);
    int turbulence_viscosity_factor(double factor);
    int turbulence_viscosity_k_omega();
    int update_k_omega_properties(double dt);
    int k_omega_time_derivatives(double *Q_rtke, double *Q_romega, double tke, double omega);
    int clear_source_vector();
    int add_inviscid_source_vector(int gtl, double omegaz=0.0);
    int add_viscous_source_vector(bool with_k_omega);
    double calculate_wall_Reynolds_number(int which_boundary);
    int store_rad_scaling_params(void);
    int rescale_Q_rE_rad(void);
    int reset_Q_rad_to_zero(void);
    double rad_scaling_ratio(void) const;
}; // end of class FV_cell

//--------------------------------------------------------------

int number_of_values_in_cell_copy(int type_of_copy);
std::string variable_list_for_cell(void);

#endif
