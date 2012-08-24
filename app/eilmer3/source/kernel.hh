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
#include "bc_defs.hh"
#include "piston.hh"
#include "thermo-interpolator.hh"
#include "radiation_transport.hh"

//-------------------------------------------------------------------

/** Macro definition to be used in signalling level of verbosity.
 *
 * 1 == Echo all file names and input parameters during initialization.
 *
 * 0 == don't be so verbose
 */
#define  ECHO_ALL  1

//-------------------------------------------------------------------

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
    int dimensions;         // 2 or 3 dimensions
    FILE *logfile;          // log file handle
    FILE *timestampfile;
    FILE *fstctimesfile;
    std::string base_file_name;
    std::string title;
    int nblock;             // number of blocks in overall simulation
    int npiston;            // number of pistons
    std::vector<Piston *> pistons;

    // Aug-2012 rework of the block-handling code for MPI.
    // We eventually want to have each task/process look after 
    // a "bag" of blocks that may not be sequentially numbered.
    std::vector<Block> bd;  // The array of vectors of blocks, holding arrays of cells.
    std::vector<Block *> my_blocks; // Collection that we can iterate over.
    
    int mpi_parallel;       // ==1 if we are using MPI parallel
    int num_mpi_proc;       // count of MPI tasks participating in the simulation
    int my_mpi_rank;        // identification for MPI process
    std::vector<int> mpi_rank_for_block; // process in which each block resides

    int step;               /* global iteration count     */
    int max_step;           /* global iteration limit     */
    int halt_now;           /* flag for premature halt    */
    int print_count; // Number of steps between writing status message to console.
    int control_count; // Number of steps between rereading .control file.

    double sim_time;        /* present simulation time    */
    double max_time;        /* final solution time, s     */
    double dt_init;         /* initial time step          */
    double dt_global;       /* simulation time step       */
    double dt_allow;        /* allowable global time step */
    double CFL;             /* target CFL (worst case)    */
    bool fixed_time_step;   /* flag for fixed time-stepping */
    bool sequence_blocks;   // if true, iterate blocks sequentially (like space-marching)
    int max_invalid_cells;  // the maximum number of bad cells (per block) 
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

    double viscous_time_delay;
    double max_mu_t_factor;
    double transient_mu_t_factor;

    double t_plot;          /* time to write next soln    */
    double t_his;           /* time to write next sample  */
    double t_fstc;          /* time to write next fluid-structure exchange data*/
    double dt_plot;         /* interval for writing soln  */
    double dt_his;          /* interval for writing sample */
    double dt_fstc;         /* interval for writing next f-s exchange data*/

    double cfl_target;      /* target CFL (worst case)    */
    int cfl_count;          /* check CFL occasionally     */
    double cfl_min;         /* current CFL minimum        */
    double cfl_max;         /* current CFL maximum        */
    double cfl_tiny;        /* smallest cfl so far        */
    double time_tiny;       /* time at which it occurred  */

    double energy_residual; /* to be monitored for steady state */
    double mass_residual;
    Vector3 energy_residual_loc, mass_residual_loc; /* location of largest value */

    std::vector<CFlowCondition*> gas_state; /* gas,flow properties */
    int n_gas_state;

    // Filter application parameters.
    int    do_filter;
    double filter_tstart;
    double filter_tend;
    double filter_dt;
    double filter_next_time;
    double filter_alpha;
    int    filter_npass;

    // variables for Andrew's time averaging routine.
    int nav;
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

    int n_heat_zone;
    double heat_time_start;
    double heat_time_stop;
    std::vector<struct CHeatZone> heat_zone;

    int n_reaction_zone;
    double reaction_time_start;
    std::vector<struct CReactionZone> reaction_zone;

    int n_turbulent_zone;
    std::vector<struct CTurbulentZone> turbulent_zone;

    std::string udf_file; // This file will contain user-defined procedures.
    int udf_source_vector_flag; // set to 1 to use (expensive) user-defined source terms
};

//---------------------------------------------------------------
// Function declarations.

global_data * get_global_data_ptr(void);
Gas_model *set_gas_model_ptr(Gas_model *gmptr);
Gas_model *get_gas_model_ptr();
int set_reaction_update(std::string file_name);
Reaction_update *get_reaction_update_ptr();
int set_energy_exchange_update( std::string file_name );
Energy_exchange_update *get_energy_exchange_update_ptr();
int set_thermo_interpolator(std::string name);
Thermo_interpolator *get_thermo_interpolator_ptr();
int set_radiation_transport_model(std::string file_name);
RadiationTransportModel *get_radiation_transport_model_ptr();
Block * get_block_data_ptr(int i);
void eilmer_finalize( void );
int set_verbose_flag( int i );
int get_verbose_flag( void );
int set_axisymmetric_flag(int ia);
int get_axisymmetric_flag(void);
int set_viscous_flag(int iv);
int get_viscous_flag(void);
double set_viscous_factor( double value );
double get_viscous_factor( void );
double incr_viscous_factor( double value );
double set_viscous_factor_increment( double value );
double get_viscous_factor_increment( void );
int set_diffusion_flag(int id);
int get_diffusion_flag(void);
int set_Xorder_flag(int ix);
int get_Xorder_flag(void);
int set_Torder_flag(int it);
int get_Torder_flag(void);
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
int set_turbulence_flag(int i);
int get_turbulence_flag(void);
int set_k_omega_flag(int ikw);
int get_k_omega_flag(void);
int set_baldwin_lomax_flag(int ibl);
int get_baldwin_lomax_flag(void);
double set_turbulence_prandtl_number(double Pr);
double get_turbulence_prandtl_number(void);
double set_turbulence_schmidt_number(double Sc);
double get_turbulence_schmidt_number(void);
int set_apply_limiter_flag(int ip);
int get_apply_limiter_flag(void);
int set_extrema_clipping_flag(int ip);
int get_extrema_clipping_flag(void);
int set_suppress_reconstruction_for_species_flag(int ip);
int get_suppress_reconstruction_for_species_flag(void);
int set_bad_cell_complain_flag(int ip);
int get_bad_cell_complain_flag(void);
int set_stringent_cfl_flag( int i );
int get_stringent_cfl_flag( void );
double set_compression_tolerance( double value );
double get_compression_tolerance( void );
double set_shear_tolerance( double value );
double get_shear_tolerance( void );

#endif
