// block.d
// Base class for blocks of cells, for use within Eilmer3.
// Peter J. 2014-07-18 first cut.

module block;
import geom;
import fvcell;
import bc;

enum
    nghost = 2; // Number of ghost cells surrounding the active cells.


class Block {
public:
    int id; // block identifier: assumed to be the same as the block number.
    bool active; // if true, block participates in the time integration
    double omegaz; // Angular velocity (in rad/s) of the rotating frame.
                   // There is only one component, about the z-axis.
    double dt_allow;            // Allowable time step
    double cfl_min, cfl_max;    // estimates of CFL number
    double mass_residual, energy_residual; // monitor these for steady state
    Vector3 mass_residual_loc, energy_residual_loc; // locations of worst case
    int hncell;            // number of sample cells
    int mncell;            // number of monitor cells
    double[] initial_T_value; // for monitor cells to check against
    FVCell[] active_cells; // collection of references to be used in foreach statements.
    BoundaryCondition[] bc; // collection of references to the boundary conditions

} // end class Block
