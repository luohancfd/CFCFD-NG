// l1d.hh
// Some general definitions that are needed at various parts of the code.
// FIX-ME Maybe these definitions should move to l_kernel.hh

#ifndef L1D_HH
#define L1D_HH

#ifndef DEBUG
#   define  DEBUG  0
#endif
// Debugging level...
// 0 = no debugging
// 1 = print message on entry to infrequently used functions
// 2 = print message on entry to frequently used functions
// 3 = dump data every step

// Case identifiers to activate special code.
// Presently, there is only a little bit in l_cell.cxx/l_slug.cxx.
#define DRUMMOND_WITH_M4_NOZZLE 27

// Definitions used for defining boundary conditions for slugs.
#define  LEFT   0
#define  RIGHT  1

// Viscous Effects Levels.
#define INVISCID                   0
#define VISCOUS_LOSS_PIPE_F        1
#define CJD_MASS_LOSS_LAMINAR      2
#define CJD_MASS_LOSS_TURBULENT    3
#define VISCOUS_LOSS_FLAT_PLATE_F  4
#define DRB_MASS_LOSS_PIPE_F       5
#define DRB_MASS_LOSS_FLAT_PLATE_F 6
#define VISCOUS_LOSS_PIPE_F_HALF_H 7

// At each end of the gas slug, the following items may be attached:
// FREE_END  : nothing at all
// SOLID_BDY : a wall with specified velocity
// PISTON    : a piston
// SLUG      : another gas slug
// SLUG_DIA  : another gas slug with diaphragm control
#define  FREE_END        0
#define  SOLID_BOUNDARY  1
#define  PISTON          2
#define  SLUG            3
#define  SLUG_DIAPHRAGM  4

// Levels of adaptivity.
#define ADAPT_NONE       0
#define ADAPT_BOTH       1
#define ADAPT_FUSE_ONLY  2
#define ADAPT_SPLIT_ONLY 3

// MINIMUM_MASS is the smallest value of cell mass for which
// Con Doolan's mass-loss model will remain active.
// If the cell mass is below this value, the model is ignored.
#define MINIMUM_MASS 1.0e-8

// NL        : number of levels in the time-stepping procedure
#define NL 2


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
