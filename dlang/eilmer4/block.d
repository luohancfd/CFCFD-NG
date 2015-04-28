// block.d
// Base class for blocks of cells, for use within Eilmer3.
// Peter J. 2014-07-18 first cut.

module block;

import std.conv;
import geom;
import fvcore;
import fvcell;
import bc;

enum
    nghost = 2; // Number of ghost cells surrounding the active cells.


class Block {
public:
    int id; // block identifier: assumed to be the same as the block number.
    string label;
    bool active; // if true, block participates in the time integration
    double omegaz; // Angular velocity (in rad/s) of the rotating frame.
                   // There is only one component, about the z-axis.
    double mass_residual, energy_residual; // monitor these for steady state
    Vector3 mass_residual_loc, energy_residual_loc; // locations of worst case
    int hncell;                 // number of sample cells
    int mncell;                 // number of monitor cells
    double[] initial_T_value; // for monitor cells to check against
    FVCell[] active_cells; // collection of references to be used in foreach statements.
    BoundaryCondition[] bc; // collection of references to the boundary conditions

    override string toString() const { return "Block(id=" ~ to!string(id) ~ ")"; }
    abstract void assemble_arrays();
    abstract void bind_faces_and_vertices_to_cells();
    abstract void identify_reaction_zones(int gtl);
    abstract void identify_turbulent_zones(int gtl);
    abstract void clear_fluxes_of_conserved_quantities();
    abstract int count_invalid_cells(int gtl);
    abstract void init_residuals();
    abstract void compute_residuals(int gtl);
    abstract double determine_time_step_size(double dt_current);
    abstract void detect_shock_points();
    abstract void compute_primary_cell_geometric_data(int gtl);
    abstract void compute_distance_to_nearest_wall_for_all_cells(int gtl);
    abstract void compute_secondary_cell_geometric_data(int gtl);
    abstract void read_grid(string filename, size_t gtl=0);
    abstract void write_grid(string filename, double sim_time, size_t gtl=0);
    abstract double read_solution(string filename);
    abstract void write_solution(string filename, double sim_time);
    abstract void write_history(string filename, double sim_time, bool write_header=false);
    abstract void set_grid_velocities(double sim_time);
    abstract void convective_flux();
    abstract void viscous_flux();
    abstract void viscous_derivatives(int gtl);
    abstract void apply_menter_boundary_correction(int ftl);
    abstract void estimate_turbulence_viscosity();
    abstract void applyPreReconAction(double t, int tLevel);
    void set_cell_dt_chem(double dt_chem)
    {
	foreach ( cell; active_cells ) cell.dt_chem = dt_chem;
    }

} // end class Block
