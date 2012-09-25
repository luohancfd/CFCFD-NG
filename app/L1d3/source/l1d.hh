/** \file l1d.hh
 * \ingroup l1d3
 * \brief Header file for the flow solver code l1d.cxx.
 *
 * \author PA Jacobs
 */

#ifndef L1D_HH
#define L1D_HH

#include <valarray>
#include <vector>
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/util/source/config_parser.hh"

#ifndef DEBUG
#   define  DEBUG  0
#endif
/*
 * Debugging level...
 * 0 = no debugging
 * 1 = print message on entry to infrequently used functions
 * 2 = print message on entry to frequently used functions
 * 3 = dump data every step
 *
 */

/*-----------------------------------------------------------------*/

              /******************************/
              /* Data Structure Definitions */
              /******************************/

/* 
 * Case identifiers to activate special code.
 * Presently, there is only a little bit in l_tstep.cxx.
 */
#define DRUMMOND_WITH_M4_NOZZLE 27

/*
 * Definitions used for defining boundary conditions for slugs.
 */
#define  LEFT   0
#define  RIGHT  1

/*
 * Viscous Effects Levels.
 */
#define INVISCID                   0
#define VISCOUS_LOSS_PIPE_F        1
#define CJD_MASS_LOSS_LAMINAR      2
#define CJD_MASS_LOSS_TURBULENT    3
#define VISCOUS_LOSS_FLAT_PLATE_F  4
#define DRB_MASS_LOSS_PIPE_F       5
#define DRB_MASS_LOSS_FLAT_PLATE_F 6
#define VISCOUS_LOSS_PIPE_F_HALF_H 7

/*
 * At each end of the gas slug, the following items may be attached:
 * FREE_END  : nothing at all
 * SOLID_BDY : a wall with specified velocity
 * PISTON    : a piston
 * SLUG      : another gas slug
 * SLUG_DIA  : another gas slug with diaphragm control
 */
#define  FREE_END        0
#define  SOLID_BOUNDARY  1
#define  PISTON          2
#define  SLUG            3
#define  SLUG_DIAPHRAGM  4

/*
 * Levels of adaptivity.
 */
#define ADAPT_NONE       0
#define ADAPT_BOTH       1
#define ADAPT_FUSE_ONLY  2
#define ADAPT_SPLIT_ONLY 3

/*
 * MINIMUM_MASS is the smallest value of cell mass for which
 * Con Doolan's mass-loss model will remain active.
 * If the cell mass is below this value, the model is ignored.
 */
#define MINIMUM_MASS 1.0e-8

/*
 * --------------------
 * Workspace dimensions.
 * --------------------
 *
 * NL        : number of levels in the time-stepping procedure
 */
#define NL 2

/*
 * -----------------------------
 * Simulation-control parameters
 * -----------------------------
 */
struct simulation_data
{
    int test_case;
    int nslug;          /* number of gas slugs        */
    int npiston;        /* number of pistons          */
    int ndiaphragm;     /* number of diaphragms       */

    int max_step;       /* global iteration limit     */
    double sim_time;    /* present simulation time    */
    double max_time;    /* final solution time, s     */
    double dt_init;     /* initial time step          */
    double dt_global;   /* simulation time step       */
    double dt_allow;    /* allowable global time step */
    double CFL;         /* target CFL (worst case)    */
    int Torder;         /* order of time-stepping     */
    int Xorder;         /* order of reconstruction    */
    int fr_chem;        /* flag to activate finite-rate chemistry */
    double k;           /* "thermal conductivity" for */
                        /* damping temperature glitch */

    int n_dt_plot;      /* number changes in plot intervals */
    std::vector<double> t_change; /* times at which dt_plot changes */
    std::vector<double> dt_plot; /* interval for writing soln  */
    std::vector<double> dt_his;  /* interval for writing sample */
    int hnloc;          /* number of history locations */
    std::vector<double> hxloc;   /* history locations          */
    int hncell;         /* history cell count (all slugs) */
};  /* end struct global_data */


// Gas states at the interfaces.
// This is used for interfacing to the Riemann solver.
struct L_flow_state
{
    Gas_data *gas;
    double u;       /* normal velocity, m/s           */
};


// Data stored in each Lagrangian cell.
struct L_cell
{
    /* GEOMETRY -- INTERFACES */
    double x;           /* Interface position, m          */
                        /* This is the interface to the   */
                        /* right of the cell midpoint.    */
    double area;        /* Interface area, m**2           */
    double pface;       /* Interface pressure (Riemann)   */
    double uface;       /* Interface velocity (Riemann)   */
    double qstar;       /* axial heat flux (fudge)        */
    double volume;      /* Cell volume, m**3              */
    double xmid;        /* location of midpoint           */
    double T_Wall;      /* specified wall temperature     */
    double K_over_L;    /* pipe-fitting loss coeff.       */

    /* CELL-AVERAGE VARIABLES */
    Gas_data *gas;   /* core values of gas properties             */
    Gas_data *ref;   /* reference values for some viscous effects */
    double u;              /* normal velocity, m/s                      */
    double shear_stress;   /* Wall shear stress, N/m**2                 */
    double heat_flux;      /* Heat flux at wall, W/m**2                 */
    double entropy;        /* Entropy referenced to 1 atm and 300K      */

    /* CONSERVED VARIABLES */
    double mass;        /* cell mass                      */
    double moment;      /* X-momentum/unit volume         */
    double Energy;      /* Total Energy/unit volume       */

    /* OTHER INTEGRATED VARIABLES */
    double L_bar;       /* This length scale is the       */
                        /* Integral of local velocity     */
                        /* for use in Mirel's model of    */
                        /* the tube-wall Boundary Layer   */

    /* A record of the cell state */
    /* (for adaptive stepping) */
    double x_old, mass_old, moment_old, Energy_old, L_bar_old;

    /* TIME DERIVATIVES */
    double DxDt[NL];    /* updates for interface posn.    */
    double DmDt[NL];    /* Mass                           */
    double DmomDt[NL];  /* X-momentum                     */
    double DEDt[NL];    /* Total Energy                   */

    double DLDt[NL];    /* Length scale, L_bar            */

    /* PRODUCTION VECTOR */
    double Q_m;         /* Mass from sources or sinks     */
    double Q_mom;       /* X-Momentum from body forces or */
                        /* wall friction                  */
    double Q_E;         /* Total Energy production or     */
                        /* transfer to the wall           */

    double dt_chem;     // Need to remember the suggested timestep.
    double dt_therm;    // and for nonequilibrium thermodynamics.
};


// Many cells make up a gas-slug.
struct slug_data
{
    /*
     * This data structure should contain everything needed for
     * a single-block solution -- both geometry and flow data
     */

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
                                /*      in l_tstep.c          */

    int adiabatic;              /* Flag for wall temperature  */
                                /* = 1: use adiabtaic temp    */
                                /* = 0: use specified temp    */

    L_flow_state *init_str; /* initial flow properties  */

    int adaptive;
    /* == 1 : the slug adapts to the flow solution by refining
     * the distribution of cells
     * * == 0 : the slug retains a fixed number of cells.
     */

    int nxdim;
    /*
     * Total number of cells for this block.
     * these will be used in the array allocation routines.
     */

    int nxmax;
    /* Maximum number of cells that will fit into the domensioned arrays
     * excluding the ghost cells. 
     */

    int nnx;
    /* Number of active cells in the X-direction */

    int nghost;
    /* Number of ghost cells around the boundaries. */

    int ixmin, ixmax;
    /*
     * These index limits are set to allow convenient access 
     * to the arrays without having to worry how many buffer 
     * cells are present.
     * Typically, ixmin <= ix <= ixmax, ixmin = 2.
     */

    double xbegin, xend;
    /*
     * The left and right end points for the slug at t=0.
     */

    double dxmin, dxmax;
    /* Maximum and minimum cell sizes for refining the cell distribution. */

    int cluster_to_end_1, cluster_to_end_2;
    double cluster_strength;
    /*
     * The initial cell distribution is based on Robert's 
     * stretching transformation.
     * If cluster_to_end_x == 1 then the cells will be clustered
     * toward end-x with the degree of clustering controlled by
     * cluster_strength.
     * cluster_strength == 0.0 : uniform distribution
     * cluster_strength > 1.0 : clustered with stronger nonuniformity
     *                          as the value approaches 1.0.
     */

    int left_bc_type, right_bc_type;
    /*
     * Type of boundary condition.  Shall be one of 
     * FREE_END
     * SOLID_BOUNDARY
     * PISTON
     * SLUG
     * SLUG_DIAPHRAGM
     */

    int left_slug_id, right_slug_id;
    int left_slug_end_id, right_slug_end_id;
    int left_piston_id, right_piston_id;
    int left_diaphragm_id, right_diaphragm_id;
    /*
     * Neighbouring gas slug, piston and diaphragm identifiers.
     * xxxx-slug_id    : number of the adjoining gas slug
     * xxxx-slug_end_id: LEFT or RIGHT end adjoins
     * xxxx-piston_id  : number (id) of the adjoining piston
     * xxxx-diaphragm  : number (id) of the adjoining diaphragm
     */

    int set_left_end_ustar, set_right_end_ustar;
    /*
     * Type of end conditions to be used when applying
     * the Riemann solver. This flag is set internally and
     * is needed by the Riemann solver.
     * 0 : standard application of the Riemann solver.
     * 1 : specified ustar
     */

    double left_ustar, right_ustar;
    /*
     * Specified interface velocities at the ends of the
     * gas slug (used to impose a velocity bc).
     */

    double left_pstar, right_pstar;
    /*
     * Computed interface pressures at the ends of the
     * gas slug (used for the piston dynamics).
     */

    struct L_cell *Cell;
    /*
     * Most of the data is stored in the preceding array. 
     * Cell[ix] = cell values
     */
};  /* end struct slug_data */


/*-----------------------------------------------------------------*/

/*
 * --------------------
 * Indexing of the data.
 * --------------------
 *
 * The following figure shows cell [ix] and its associated 
 * interfaces.
 *
 *
 *             +-----------------------------+
 *   West      |         cell center         |  East 
 *   face      |             [ix]            |  face
 *   [ix-1]    x              o              x  [ix]
 *             |                             |
 *             |                             |
 *             +-----------------------------+
 *
 *
 * Thus...
 * ----
 * Active cells are indexed as LCell[ix], where
 * ixmin <= ix <= ixmax.
 *
 * Acitve vertical interfaces are indexed as X[ix], where
 * ixmin-1 <= ix <= ixmax.
 *
 * Space for ghost cells is available outside these ranges but
 * within the dimensioned range  0 <= ix <= nxdim-1.
 *
 * ----------------------------------------------------------- */

#endif
