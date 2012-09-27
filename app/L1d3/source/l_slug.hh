// l_slug.hh

#ifndef L_SLUG_HH
#define L_SLUG_HH

#include <vector>
#include <string>
#include <stdio.h>
#include "l1d.hh"
#include "l_tube.hh"
#include "l_kernel.hh"
#include "l_cell.hh"

class GasSlug {
// Many Lagrangian cells make up a gas-slug.
public:
    double dt;                  /* time step for this block   */
    double dt0;                 /* initial time step          */
    double dt_allow;            /* Allowable time step        */
    double cfl_target;          /* desired CFL number         */
    double sim_time;            /* simulation time in seconds */
    double cfl_min, cfl_max;    /* estimates of CFL number    */
    double cfl_tiny, time_tiny; /* the smallest so far...     */

    double residual;            /* residual monitor for       */
                                /* steady state               */

    int max_steps;              /* max number of time steps   */
                                /* on this grid               */

    int Xorder;                 /* spatial order 1 or 2       */
    int Torder;                 /* temporal order 1 or 2      */

    int hncell;                 /* number of sample cells     */
    std::vector<int> hxcell;   /* location of sample cell    */

    int viscous_effects;        /* Flag for including viscous */
                                /* effects:                   */
                                /* = 1: include them          */
                                /* = 0: Do not include them   */
                                /* = 2: laminar mass-loss     */
                                /* = 3: turbulent mass-loss   */
                                /*      See L_source_vector() */
                                /*      in l_slug.c           */

    int adiabatic;              /* Flag for wall temperature  */
                                /* = 1: use adiabtaic temp    */
                                /* = 0: use specified temp    */

    L_flow_state *init_str; /* initial flow properties  */

    int adaptive;
    // == 1 : the slug adapts to the flow solution by refining
    // the distribution of cells
    // == 0 : the slug retains a fixed number of cells.

    int nxdim;
    // Total number of cells for this block.
    // these will be used in the array allocation routines.

    int nxmax;
    // Maximum number of cells that will fit into the domensioned arrays
    // excluding the ghost cells. 

    int nnx;
    // Number of active cells in the X-direction

    int nghost;
    // Number of ghost cells around the boundaries.

    int ixmin, ixmax;
    // These index limits are set to allow convenient access 
    // to the arrays without having to worry how many buffer 
    // cells are present.
    // Typically, ixmin <= ix <= ixmax, ixmin = 2.

    double xbegin, xend;
    /*
     * The left and right end points for the slug at t=0.
     */

    double dxmin, dxmax;
    /* Maximum and minimum cell sizes for refining the cell distribution. */

    int cluster_to_end_1, cluster_to_end_2;
    double cluster_strength;
    // The initial cell distribution is based on Robert's 
    // stretching transformation.
    // If cluster_to_end_x == 1 then the cells will be clustered
    // toward end-x with the degree of clustering controlled by
    // cluster_strength.
    // cluster_strength == 0.0 : uniform distribution
    // cluster_strength > 1.0 : clustered with stronger nonuniformity
    //                          as the value approaches 1.0.

    int left_bc_type, right_bc_type;
    // Type of boundary condition.  Shall be one of 
    // FREE_END
    // SOLID_BOUNDARY
    // PISTON
    // SLUG
    // SLUG_DIAPHRAGM

    int left_slug_id, right_slug_id;
    int left_slug_end_id, right_slug_end_id;
    int left_piston_id, right_piston_id;
    int left_diaphragm_id, right_diaphragm_id;
    // Neighbouring gas slug, piston and diaphragm identifiers.
    // xxxx-slug_id    : number of the adjoining gas slug
    // xxxx-slug_end_id: LEFT or RIGHT end adjoins
    // xxxx-piston_id  : number (id) of the adjoining piston
    // xxxx-diaphragm  : number (id) of the adjoining diaphragm

    int set_left_end_ustar, set_right_end_ustar;
    // Type of end conditions to be used when applying
    // the Riemann solver. This flag is set internally and
    // is needed by the Riemann solver.
    // 0 : standard application of the Riemann solver.
    // 1 : specified ustar

    double left_ustar, right_ustar;
    // Specified interface velocities at the ends of the
    // gas slug (used to impose a velocity bc).

    double left_pstar, right_pstar;
    // Computed interface pressures at the ends of the
    // gas slug (used for the piston dynamics).

    std::vector<LCell> Cell;
    // Most of the data is stored in the preceding array. 
    // Cell[ix] = cell values

    GasSlug(int indx, SimulationData& SD, 
	    std::string config_file_name, int echo_input);
    GasSlug(const GasSlug& gs);
    ~GasSlug(void);
    int set_index_range(void);
    int read_state(FILE* infile);
    int write_state(FILE* outfile);
    int compute_areas(TubeModel *tube);
    double total_energy();
    int maximum_p(double *p_max, double *x_max);
    int fill_data();
    int encode_conserved();
    int decode_conserved();
    int set_chemistry_timestep(double dt);
    int set_thermal_timestep(double dt);
    int chemical_increment(double dt);
    int source_vector();
    int axial_heat_flux(double k);
    int adjust_end_cells();
    int apply_rivp();
    int time_derivatives(int time_level);
    int record_state();
    int restore_state();
    int predictor_step();
    int corrector_step();
    int check_cells(int js);
    int check_cfl();
    int interpolate_cell_data(double xloc, LCell& icell);
    double end_pressure(int which_end, double dx);
    int end_properties(int which_end, double dx,
		       double* total_mass, struct L_flow_state& Q);
}; // end class GasSlug

double min_increment( double y0, double y1, double y2 );
double limit_to_within(double v, double a, double b);

#endif
