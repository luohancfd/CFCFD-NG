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
// #include "flux_calc.hh"

/*--------------------------------------------------------------*/

/** These values are used by the cell data copying functions
 * to specify which parts of the cell data to copy.
 */
#define COPY_ALL_CELL_DATA 0
#define COPY_FLOW_STATE    1
#define COPY_CELL_LENGTHS  2
#define COPY_INTERFACE_DATA  3

/** Symbolic labels and indices for the cell's faces. 
 *
 * The names of the faces of the structured-grid blocks will be the same.
 */
#define NORTH  0
#define EAST   1
#define SOUTH  2
#define WEST   3
#define TOP    4
#define BOTTOM 5

/// States that a cell-center may be in with respect to pistons.
#define NORMAL_CELL   0
#define MASKED_CELL   1
#define SHADOWED_CELL 2
/// States that cell-interfaces may be in with respect to pistons.
#define NORMAL_IFACE  0
#define MASKED_IFACE  1

/** \brief Number of levels in the time-stepping procedure.
 *
 * Used below to dimension some time-derivative arrays.
 */
#define NL 4
#define NI 6
#define NV 8

/// We might update the k-omega properties in with the main predictor-corrector
/// time-stepping function or we might choose to update it separately, 
/// like the chemistry update.
#define SEPARATE_UPDATE_FOR_K_OMEGA_SOURCE 1
/// Minimum values for turbulent kinetic energy (m^2/s^2) and frequency (1/s)
/// for applying limiters in the k-omega model.
#define SMALL_TKE 0.1
#define SMALL_OMEGA 1.0

int get_face_index(const std::string name);
std::string  get_face_name(int face_index);

/** \brief Flow state describes the local flow but not the cell or interface geometry. */
class FlowState {
public:
    Gas_data *gas;       ///< \brief gas thermo state
    Vector3 vel;         ///< \brief velocity vector, m/s
    Vector3 shock_vel;   ///< \brief shock velocity vector, m/s
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
    ~FlowState();
    int print();
    int copy_values_from(FlowState &src);
    int copy_values_from(CFlowCondition &src);
    int average_values_from(FlowState &src0, FlowState &src1, bool with_diff_coeff);
    int average_values_from(FlowState &src0, double alpha0, FlowState &src1, double alpha1, bool with_diff_coeff);
    double * copy_values_to_buffer(double *buf);
    double * copy_values_from_buffer(double *buf);
    int BGK_equilibrium(void);
};

/// \brief Conserved quantities, per unit volume, are stored for each cell.
class ConservedQuantities {
public:
    double mass;                  ///< \brief density, kg/m**3
    Vector3 momentum;             ///< \brief momentum/unit volume
    Vector3 B;                    ///< \brief magnetic field, Tesla
    double total_energy;          ///< \brief total energy/unit volume
    std::vector<double> massf;    ///< \brief mass of species/unit volume
    std::vector<double> energies; ///< \brief vib. energy / unit-volume
    double tke;                   ///< \brief turbulent kinetic energy
    double omega;                 ///< \brief omega from k-omega turbulence model
    std::vector<double> G; ///< \brief velocity dist. partial densities, kg/m**3
    std::vector<double> H; ///> \brief velocity dist. partial densities, (kg*s**2)/(m**5)
    //
    ConservedQuantities(Gas_model *gm);
    ~ConservedQuantities();
    int print();
    int copy_values_from(ConservedQuantities &src);
    int clear_values();
};

/// \brief Inviscid and viscous fluxes are transported (between cells) across cell interfaces.
class FV_Interface {
public:
    int    id;
    int status; // state of the interface with respect to pistons
    // Geometry
    Vector3 pos;       /* position of the (approx) midpoint */
    Vector3 vel;       /* velocity, m/s    */
    double Ybar;       /* Y-coordinate of the mid-point     */
    double length;     /* Interface length in the x,y-plane */
    double ar[NL];     /* Area m**2 in 2D geometries for each integration time level */
    double area;       /* Area m**2 in 2D geometries        */
                       /* Area per radian in axisymmetric g */
    Vector3 n;         /* Direction cosines for unit normal */
    Vector3 t1;        /* tangent vector 1 (aka p)          */
    Vector3 t2;        /* tangent vector 2 (aka q)          */
    // Flow
    FlowState *fs;
    ConservedQuantities *F; // flux conserved quantity per unit area
#   if WITH_IMPLICIT == 1
    // Point-implicit variables
    double Lambda[6][2], R[6][6], R_inv[6][6], sl[6][2], sl_min[6][2];
    double A[6][6], J[6][6], hl[6][2], gl[6][2];
#   endif
    //
    FV_Interface(Gas_model *gm);
    ~FV_Interface();
    int print(int to_stdout);
    int copy_values_from(FV_Interface &src, int type_of_copy);
}; // end of class FV_Interface


/// \brief Cell vertex holds data for the calculation of gradients for the viscous terms.
class FV_Vertex {
public:
    int    id;
    // Geometry
    Vector3 pos;  // x,y,z-Coordinates, m
    Vector3 position[NL];  // x,y,z-Coordinates for different integration time levels, m
    Vector3 vel;  // velocity, m/s
    Vector3 velocity[NL];  // velocity for different integration time levels, m/s
    double area;  // x,y-plane area of secondary cells
    double volume;  // volume of 3D secondary cells
    // Derivatives of primary-cell variables.
    double dudx, dudy, dudz; // velocity derivatives
    double dvdx, dvdy, dvdz; // velocity derivatives
    double dwdx, dwdy, dwdz; // velocity derivatives
    std::vector<double> dTdx, dTdy, dTdz; // Temperature derivatives
    double dtkedx, dtkedy, dtkedz; // turbulence kinetic energy
    double domegadx, domegady, domegadz; // pseudo vorticity for k-omega turbulence
    std::vector<double> dfdx, dfdy, dfdz; // mass fraction derivatives
    //
    FV_Vertex(Gas_model *gm);
    ~FV_Vertex();
    int copy_values_from(FV_Vertex &src);
};


/// \brief Cell-averaged data forms the core of the flow field data.
class FV_Cell {
public:
    int id;  // may be used as a global identity in VTK files
    int status; // state of the cell with respect to pistons
    int fr_reactions_allowed; // ==1, will call chemical_increment (also thermal_increment)
    double dt_chem; // acceptable time step for finite-rate chemistry
    double dt_therm; // acceptable time step for thermal relaxation
    int in_turbulent_zone; // ==1, we will keep the turbulence viscosity
    double base_qdot;       // base-level of heat addition to cell, W/m**3
    // Geometry
    Vector3 pos;            /* Centre x,y,z-coordinates, m    */
    Vector3 position[NL];   /* Centre x,y,z-coordinates for different integration time levels, m	*/
    double vol[NL];         /* Cell volume (unit depth), m**3 */
    double volume;
    double ar[NL];          /* (x,y)-plane area for different integration time levels, m**2         */
    double area;            /* (x,y)-plane area, m**2         */
    double uf;              /* uncovered fraction */
    double iLength;         /* length in the i-index direction */
    double jLength;         /* length in the j-index direction */
    double kLength;         /* length in the k-index direction */
    double L_min;           /* minimum length scale for cell  */
    double distance_to_nearest_wall; /* for turbulence model correction. */
    double half_cell_width_at_wall;  /* ditto                            */
    FV_Cell *cell_at_nearest_wall;   /* ditto                            */
    // Connections
    FV_Interface *iface[6];  // pointers to defining interfaces of cell
    FV_Vertex *vtx[8];  // pointers to vertices for quad (2D) and hexahedral (3D) cells
    // Flow
    FlowState *fs;
    ConservedQuantities *U, *U_old; // current conserved quantities, and at time-step start
    ConservedQuantities *dUdt[NL];  // time derivatives
    ConservedQuantities *Q;         // source (or production) terms
    // Terms for loose-coupling of radiation.
    double Q_rad_org;
    double f_rad_org;
    double Q_rE_rad; //Rate of energy addition to cell via radiation.
    // Data for computing residuals.
    double rho_at_start_of_step, rE_at_start_of_step;
#   if WITH_IMPLICIT == 1
    /* POINT-IMPLICIT VARIABLES */
    double piM[6][6], pir[6][2], qL[6][2], g_Ll[6][2], gjvec[6][2], gjmtx[6][6];
#   endif
    //
    FV_Cell(Gas_model *gm);
    ~FV_Cell();
    int print();
    int point_is_inside(Vector3 &p, int dimensions);
    int copy_values_from(CFlowCondition &src);
    int copy_values_from(FV_Cell &src, int type_of_copy);
    double * copy_values_to_buffer(double *buf, int type_of_copy);
    double * copy_values_from_buffer(double *buf, int type_of_copy);
    int replace_flow_data_with_average(FV_Cell *src[], int ncell);
    int scan_values_from_string(char *bufptr);
    std::string write_values_to_string();
    int scan_BGK_from_string(char *bufptr);
    std::string write_BGK_to_string();
    int impose_chemistry_timestep(double dt);
    int impose_thermal_timestep(double dt);
    int set_fr_reactions_allowed(int flag);
    int record_conserved(void);
    int restore_conserved(void);
    int encode_conserved(double omegaz=0.0);
    int decode_conserved(double omegaz=0.0);
    int decode_conserved(int time_level, double omegaz=0.0);
    int check_flow_data(void);
    int init_time_level_geometry( void );
    int get_current_time_level_geometry( int time_level );
    int time_derivatives(int time_level, int dimensions);
    int predictor_update(double dt);
    int corrector_update(double dt);
    int rk3_update(double dt);
    int chemical_increment(double dt);
    int thermal_increment(double dt);
    double signal_frequency(int dimensions);
    int turbulence_viscosity_zero(void);
    int turbulence_viscosity_zero_if_not_in_zone(void);
    int turbulence_viscosity_limit(double factor);
    int turbulence_viscosity_factor(double factor);
    int turbulence_viscosity_k_omega(void);
    int update_k_omega_properties(double dt);
    int k_omega_time_derivatives(double *Q_rtke, double *Q_romega, double tke, double omega);
    int inviscid_source_vector(int time_level, double omegaz=0.0);
    int viscous_source_vector(void);
    double calculate_wall_Reynolds_number(int which_boundary);
    int store_rad_scaling_params(void);
    int rescale_Q_rE_rad(void);
    int reset_Q_rad_to_zero(void);
    double rad_scaling_ratio(void);
}; // end of class FV_cell

//--------------------------------------------------------------

int number_of_values_in_cell_copy(int type_of_copy);
std::string variable_list_for_cell(void);

int one_d_interp(FV_Cell &cL1, FV_Cell &cL0, 
		 FV_Cell &cR0, FV_Cell &cR1, 
		 double cL1Length, double cL0Length, 
		 double cR0Length, double cR1Length, 
		 FlowState &Lft, FlowState &Rght);

int mach_weighted_one_d_interp(FV_Cell &cL1, FV_Cell &cL0, 
			       FV_Cell &cR0, FV_Cell &cR1, 
			       double cL1Length, double cL0Length, 
			       double cR0Length, double cR1Length, 
			       FlowState &Lft, FlowState &Rght);

int onesided_interp(FV_Cell &cL0, FV_Cell &cR0, FV_Cell &cR1,
		    double cL0Length, double cR0Length, double cR1Length,
		    FlowState &Rght );

int one_d_interior_interp(FV_Cell &cL0, FV_Cell &cR0, FV_Cell &cR1, FV_Cell &cR2,
			  double cL0Length, double cR0Length, double cR1Length, double cR2Length,
			  FlowState &Lft, FlowState &Rght );

#endif
