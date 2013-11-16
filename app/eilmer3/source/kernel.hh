/// \file kernel.hh
/// \ingroup eilmer3
/// \brief Header file for the kernel module.
///
/// Contains some configuration elements, also.

#ifndef KERNEL_HH
#define KERNEL_HH

#include <string>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/kinetics/reaction-update.hh"
#include "../../../lib/gas/kinetics/energy-exchange-update.hh"
#include "c-flow-condition.hh"
#include "flux_calc.hh"
#include "cell.hh"
#include "block.hh"
#include "piston.hh"
#include "radiation_transport.hh"

//-------------------------------------------------------------------

// 1 == Echo all file names and input parameters during initialization.
// 0 == don't be so verbose
const int ECHO_ALL = 1;

enum turbulence_model_t {TM_NONE, TM_BALDWIN_LOMAX, TM_K_OMEGA, TM_SPALART_ALLMARAS};

struct CHeatZone {
    double qdot;  // rate of heat addition in W/m**3
    double x0, y0, z0;
    double x1, y1, z1;
};

struct CReactionZone {
    double x0, y0, z0;
    double x1, y1, z1;
};

struct CTurbulentZone {
    double x0, y0, z0;
    double x1, y1, z1;
};

/** \brief Global data structure for control of the overall time-stepping. */
struct global_data
{
    size_t dimensions;      // 2 or 3 dimensions
    FILE *logfile;          // log file handle
    FILE *timestampfile;
    FILE *fstctimesfile;
    std::string base_file_name;
    std::string title;
    size_t nblock;             // number of blocks in overall simulation
    size_t npiston;            // number of pistons
    std::vector<Piston *> pistons;

    // Aug-2012 rework of the block-handling code for MPI.
    // We eventually want to have each task/process look after 
    // a "bag" of blocks that may not be sequentially numbered.
    std::vector<Block> bd;  // The array of vectors of blocks, holding arrays of cells.
    std::vector<Block *> my_blocks; // Collection that we can iterate over.
    
    bool mpi_parallel;      // ==1 if we are using MPI parallel
    int num_mpi_proc;       // count of MPI tasks participating in the simulation
    int my_mpi_rank;        // identification for MPI process
    std::vector<int> mpi_rank_for_block; // process in which each block resides

    size_t step;            /* global iteration count     */
    size_t max_step;        /* global iteration limit     */
    size_t t_level;         /* global time level within update */
    int halt_now;           /* flag for premature halt    */
    size_t print_count; // Number of steps between writing status message to console.
    size_t control_count; // Number of steps between rereading .control file.

    double sim_time;        /* present simulation time    */
    double max_time;        /* final solution time, s     */
    double dt_init;         /* initial time step          */
    double dt_global;       /* simulation time step       */
    double dt_allow;        /* allowable global time step */
    double CFL;             /* target CFL (worst case)    */
    double dt_max;          // Maximum allowable time-step, after all other considerations.
    bool fixed_time_step;   /* flag for fixed time-stepping */
    int Xorder; // Low order reconstruction (1) uses just the cell-centre data as left- and right-
                // flow properties in the flux calculation.
                // High-order reconstruction (2) adds a correction term to the cell-centre values
                // to approach something like a piecewise-quadratic interpolation between the
                // cell centres.

    bool sequence_blocks;   // if true, iterate blocks sequentially (like space-marching)
    size_t max_invalid_cells;  // the maximum number of bad cells (per block) 
                            // which will be tolerated without too much complaint.
    double dt_reduction_factor; 
    /*
     * If an attempt at a time step fails because of invalid cells,
     * the time step is re-attempted with a smaller time step.
     * This reduction factor is somewhat arbitrary and can now be set
     * by the user's imput script.
     * A factor of 0.5 would seem to be not enough but a factor of
     * 0.1 would seem too expensive.  We have settled on a default of 0.2.
     */

    /// We might update some properties in with the main convective-terms
    /// time-stepping function or we might choose to update them separately, 
    /// like the chemistry update.
    bool separate_update_for_viscous_terms;
    bool separate_update_for_k_omega_source;

    double viscous_time_delay;
    double diffusion_time_delay;
    turbulence_model_t turbulence_model;
    double turbulence_prandtl;
    double turbulence_schmidt;
    double max_mu_t_factor;
    double transient_mu_t_factor;
    double shock_fitting_speed_factor;

    double t_plot;          /* time to write next soln    */
    size_t write_at_step;   /* update step at which to write a solution, 0=don't do it */
    double t_his;           /* time to write next sample  */
    double t_fstc;          /* time to write next fluid-structure exchange data*/
    double t_shock;         /* time to next adapt grid to shock    */
    double dt_shock;        /* interval for running shock adapting algorithm  */
    double dt_plot;         /* interval for writing soln  */
    double dt_his;          /* interval for writing sample */
    double dt_fstc;         /* interval for writing next f-s exchange data*/

    double cfl_target;      /* target CFL (worst case)    */
    size_t cfl_count;          /* check CFL occasionally     */
    double cfl_min;         /* current CFL minimum        */
    double cfl_max;         /* current CFL maximum        */
    double cfl_tiny;        /* smallest cfl so far        */
    double time_tiny;       /* time at which it occurred  */

    double energy_residual; /* to be monitored for steady state */
    double mass_residual;
    Vector3 energy_residual_loc, mass_residual_loc; /* location of largest value */

    std::vector<CFlowCondition*> gas_state; /* gas,flow properties */
    size_t n_gas_state;

    // Filter application parameters.
    int    filter_flag;
    double filter_tstart;
    double filter_tend;
    double filter_dt;
    double filter_next_time;
    double filter_mu;
    size_t filter_npass;

    // variables for Andrew's time averaging routine.
    size_t nav;
    double tav_0, tav_f, dtav;
    double tav;
    int do_tavg;

    // variables for profile recording.
    int do_record;
    int block_record, x_record, n_record, step_record;
    double t_record, t_start;

    // variables for time dependent profile input.
    double t_old, ta_old;
    int n_vary;
    int do_vary_in;

    // these variables are all for special cases
    int secondary_diaphragm_ruptured, d2_ruptured, dn_ruptured;
    double diaphragm_rupture_time, diaphragm_rupture_diameter;
    int diaphragm_block, drummond_progressive;

    size_t n_heat_zone;
    double heat_time_start;
    double heat_time_stop;
    std::vector<struct CHeatZone> heat_zone;

    size_t n_reaction_zone;
    double reaction_time_start;
    double T_frozen; // temperature (in K) below which reactions are frozen
    double T_frozen_energy; // temperature (in K) below which energy exchanges are skipped
    std::vector<struct CReactionZone> reaction_zone;

    size_t n_turbulent_zone;
    std::vector<struct CTurbulentZone> turbulent_zone;

    std::string udf_file; // This file will contain user-defined procedures.
    int udf_source_vector_flag; // set to 1 to use (expensive) user-defined source terms

#   ifdef CONJUGATE_HT
    // variables related to a wall model for conjugate heat transfer
    bool conjugate_ht_active;
    Wall_model *wm;
    std::vector<double> T_wall;
    std::vector<double> q_wall;
    std::vector<int> recvounts;
    std::vector<int> displs;
#   endif
};

//---------------------------------------------------------------
// Function declarations.

std::string get_revision_string();
global_data * get_global_data_ptr(void);
Gas_model *set_gas_model_ptr(Gas_model *gmptr);
Gas_model *get_gas_model_ptr();
int set_reaction_update(std::string file_name);
Reaction_update *get_reaction_update_ptr();
int set_energy_exchange_update( std::string file_name );
Energy_exchange_update *get_energy_exchange_update_ptr();
int set_radiation_transport_model(std::string file_name);
RadiationTransportModel *get_radiation_transport_model_ptr();
Block * get_block_data_ptr(size_t i);
void eilmer_finalize( void );

int set_verbose_flag( int i );
int get_verbose_flag( void );
int set_axisymmetric_flag(int ia);
int get_axisymmetric_flag(void);
int set_shock_fitting_flag(int iw);
int get_shock_fitting_flag(void);
int set_shock_fitting_decay_flag(int iw);
int get_shock_fitting_decay_flag(void);
int set_moving_grid_flag(int iw);
int get_moving_grid_flag(void);
int set_write_vertex_velocities_flag(int iw);
int get_write_vertex_velocities_flag(void);
int set_filter_flag(int iw);
int get_filter_flag(void);
int set_viscous_flag(int iv);
int get_viscous_flag(void);
int set_viscous_upwinding_flag(int iw);
int get_viscous_upwinding_flag(void);
double set_viscous_factor( double value );
double get_viscous_factor( void );
double incr_viscous_factor( double value );
double set_viscous_factor_increment( double value );
double get_viscous_factor_increment( void );
int set_diffusion_flag(int id);
int get_diffusion_flag(void);
double set_diffusion_factor( double value );
double get_diffusion_factor( void );
double incr_diffusion_factor( double value );
double set_diffusion_factor_increment( double value );
double get_diffusion_factor_increment( void );
double set_heat_factor( double value );
double get_heat_factor( void );
double incr_heat_factor( double value );
double set_heat_factor_increment( double value );
double get_heat_factor_increment( void );
int set_reacting_flag(int iv);
int get_reacting_flag(void);
int set_energy_exchange_flag(int iv);
int get_energy_exchange_flag(void);
int set_radiation_flag(int iv);
int get_radiation_flag(void);
int set_radiation_update_frequency(int ruf);
int get_radiation_update_frequency(void);
int set_implicit_flag(int imf);
int get_implicit_flag(void);
int set_bad_cell_complain_flag(int ip);
int get_bad_cell_complain_flag(void);
int set_stringent_cfl_flag( int i );
int get_stringent_cfl_flag( void );
double set_compression_tolerance( double value );
double get_compression_tolerance( void );
double set_shear_tolerance( double value );
double get_shear_tolerance( void );

int set_mhd_flag(int imf);
int get_mhd_flag(void);
int set_BGK_flag(int i);
int get_BGK_flag( void );
size_t set_velocity_buckets(size_t i);
size_t get_velocity_buckets( void );
Vector3 get_vcoord(int i);
std::vector<Vector3> *get_vcoords_ptr(void);
double get_vweight(int i);
std::vector<double> *get_vweights_ptr(void);
std::string get_name_of_turbulence_model(turbulence_model_t my_model);

int set_electric_field_work_flag(int iefw );
int get_electric_field_work_flag();
#endif
