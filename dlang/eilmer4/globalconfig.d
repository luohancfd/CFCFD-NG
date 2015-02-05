/** globalconfig.d
 * A place to keep the configuration details for the simulation.
 *
 * Author: Peter J. and Rowan G.
 * First code: 2014-07-18
 * Resumed work: 2015-02-05
 */

module globalconfig;

import std.conv;
import std.stdio;

import geom;
import gas;
import fvcore;

enum
    tm_none = 0,
    tm_baldwin_lomax = 1,
    tm_k_omega = 2,
    tm_spalart_allmaras = 3;
string[] turbulence_model_names = ["none", "baldwin_lomax", "k_omega", "spalart_allmaras"];

class BlockZone {
    Vector3 _p0;
    Vector3 _p1;
    this(in Vector3 p0, in Vector3 p1) {
	_p0 = p0;
	_p1 = p1;
    }
    override string toString() {
	return "BlockZone(" ~ _p0.toString() ~ "," ~ _p1.toString() ~ ")";
    }
    bool is_inside(in Vector3 p, int dimensions) {
	if ( p.x >= _p0.x && p.x <= _p1.x &&
	     p.y >= _p0.y && p.y <= _p1.y ) {
	    if ( dimensions == 2 ) {
		return true;
	    } else if ( p.z >= _p0.z && p.z <= _p1.z ) {
		return true;
	    }
	}
	return false;
    } // end is_inside()
}

class IgnitionZone : BlockZone {
    double _Tig; // temperature to apply within reaction_update ensure ignition
    this(in Vector3 p0, in Vector3 p1, double Tig) {
	super(p0, p1);
	_Tig = Tig;
    }
    override string toString() {
	return "IgnitionZone(" ~ _p0.toString() ~ "," ~ _p1.toString() ~ ","
	    ~ to!string(_Tig) ~ ")";
    }
}


final class GlobalConfig {
    static string base_file_name = "job"; // Change this to suit at run time.
    static string title = "Eilmer4 simulation"; // Change this to suit at run time.
    static GasModel gmodel;

    static int nBlocks; // Number of blocks in the overall simulation.
    static int dimensions = 2; // or 3
    static bool axisymmetric = false;
    // Low order reconstruction (1) uses just the cell-centre data as left- and right-
    // flow properties in the flux calculation.
    // High-order reconstruction (2) adds a correction term to the cell-centre values
    // to approach something like a piecewise-quadratic interpolation between the
    // cell centres.
    static int Xorder = 2; 

    static bool viscous = false; 
    // If true, viscous effects are included in the gas-dynamic update.
    // A factor to scale the viscosity in order to achieve a soft start. 
    // The soft-start for viscous effects may be handy for impulsively-started flows.
    // A value of 1.0 means that the viscous effects are fully applied.
    static double viscous_factor = 1.0;
    // The amount by which to increment the viscous factor during soft-start.
    static double viscous_factor_increment = 0.01;
    static double viscous_time_delay = 0.0;

    // When the diffusion is calculated is treated as part of the viscous calculation:
    //   false for neglecting multicomponent diffusion, 
    //   true when considering the diffusion 
    static bool diffusion = false; 
    // A factor to scale the diffusion in order to achieve a soft start.
    // The soft-start for diffusion effects may be handy for impulsively-started flows.
    // Note that this is separate to viscous effects.
    static double diffusion_factor = 1.0;
    // The amount by which to increment the diffusion factor during soft-start.
    static double diffusion_factor_increment = 0.01;
    static double diffusion_time_delay = 0.0;
    // The Lewis number when using the constant Lewis number diffusion model
    static double diffusion_lewis = 1.0;
    // The Schmidt number when using the constant Schmidt number diffusion model
    static double diffusion_schmidt = 0.7;

    static int turbulence_model = tm_none;
    static double turbulence_prandtl = 0.89;
    static double turbulence_schmidt = 0.75;
    static double max_mu_t_factor = 300.0;
    static double transient_mu_t_factor = 1.0;
    static BlockZone[] turbulent_zones;

    static bool moving_grid = false;
    static bool write_vertex_velocities = false;

    // Set the tolerance in relative velocity change for the shock detector.
    // This value is expected to be a negative number (for compression)
    // and not too large in magnitude.
    // We have been using a value of -0.05 for years, based on some
    // early experiments with the sod and cone20 test cases, however,
    // the values may need to be tuned for other cases, especially where
    // viscous effects are important.
    static double compression_tolerance = -0.30;

    // Default flow-data reconstruction includes interpolation of density 
    // and internal energy.  Other options for the thermodunamic properties
    // to be interpolated are pressure+temperature, density+temperature and
    // density+pressure.
    static InterpolateOption thermo_interpolator = InterpolateOption.rhoe;
    static bool apply_limiter = true;
    static bool extrema_clipping = true;
    static bool interpolate_in_local_frame = true;

    // Default flux calculator is the adaptive mix of ausmdv and efm.
    static FluxCalculator flux_calculator = FluxCalculator.adaptive;

    // Set the tolerance to shear when applying the adaptive flux calculator.
    // We don't want EFM to be applied in situations of significant shear.
    // The shear value is computed as the tangential-velocity difference across an interface
    // normalised by the local sound speed.
    static double shear_tolerance = 0.20;

    // Reference free-stream Mach number, for use in the ausm_plus_up flux calculator.
    // Choose a value for M_inf that is good for low Mach numbers.
    // To be strictly correct, we should set this at run time
    // if an M_inf value is easily defined.
    static double M_inf = 0.01;

    // Turning on the reactions activates the chemical update function calls.
    // Chemical equilibrium simulations (via Look-Up Table) does not use this
    // chemical update function call.
    static bool reacting = false;

    // With this flag on, finite-rate evolution of the vibrational energies 
    // (and in turn the total energy) is computed.
    static bool thermal_energy_exchange = false;

    static bool radiation = false;
    static int radiation_update_frequency; // = 1 for every time-step
    static bool radiation_scaling = false;

    static bool electric_field_work;

    // For Daryl Bond and Vince Wheatley's MHD additions.
    static bool MHD;

    // A flag for turning on the BGK non-equilibrium gas solver:
    //   BGK == 0: OFF
    //   BGK == 1: ON, do not try to import velocity distribution values
    //   BGK == 2: ON, read in velocity distribution values from "flow" file
    static int BGK = 0;

    static double ignition_time_start = 0.0;
    static double ignition_time_stop = 0.0;
    static IgnitionZone[] ignition_zones;
    static bool ignition_zone_active = false;

    static double reaction_time_start = 0.0;
    static double T_frozen; // temperature (in K) below which reactions are frozen
    static double T_frozen_energy; // temperature (in K) below which energy exchanges are skipped
    static BlockZone[] reaction_zones;

    static int step;                // iteration count
    static int max_step;            // iteration limit
    static int t_level;             // time level within update
    static int halt_now = 0;        // flag for premature halt
    static bool halt_on_large_flow_change = false; 
    // Set to true to halt simulation when any
    // monitor point sees a large flow change.
    static double tolerance_in_T;   // Temperature change for the flow change.
    static int print_count; // Number of steps between writing messages to console.
    static int control_count; // Number of steps between rereading .control file.

    static int verbosity_level = 1; 
    // Messages have a hierarchy:
    // 0 : only error messages will be omitted
    // 1 : emit messages that are useful for a long-running job (default)
    // 2 : plus verbose init messages
    // 3 : plus verbose boundary condition messages
    // 4 : temporary messages for debugging

    /// When decoding the array of conserved quantities, 
    /// the temperature or the density may try to go negative.  
    /// If it does and adjust_invalid_cell_data == true, the cell data
    /// is adjusted to make it reasonable.
    static bool adjust_invalid_cell_data = false;

    static double sim_time = 0.0;  // present simulation time, tracked by code
    static double max_time;        // final solution time, s, set by user
    static double dt_init;         // initial time step set by user
    static double dt_global;       // simulation time step determined by code
    static double dt_allow;        // allowable global time step determined by code
    static double cfl_value = 0.5; // target CFL number (worst case) set by user
    static bool stringent_cfl = false; 
    // If true, assume the worst with respect to cell geometry and wave speed.
    static double dt_max; // Maximum allowable time-step, after all other considerations.
    static bool fixed_time_step = false; // set true to fix dt_allow

    static int max_invalid_cells = 0;  
    // The maximum number of bad cells (per block) 
    // which will be tolerated without complaint.
    static double dt_reduction_factor = 0.2; 
    // If an attempt at a time step fails because of invalid cells,
    // the time step is re-attempted with a smaller time step.
    // This reduction factor is somewhat arbitrary and can now be set
    // by the user's imput script.
    // A factor of 0.5 would seem to be not enough but a factor of
    // 0.1 would seem too expensive.  We have settled on a default of 0.2.

    /// We might update some properties in with the main convective-terms
    /// time-stepping function or we might choose to update them separately, 
    /// like the chemistry update.
    static bool separate_update_for_viscous_terms = false;
    static bool separate_update_for_k_omega_source = false;

    static double t_plot;          // time to write next soln
    static size_t write_at_step;   // update step at which to write a solution, 0=don't do it
    static double t_his;           // time to write next sample
    static double dt_plot;         // interval for writing soln
    static double dt_his;          // interval for writing sample

    static double cfl_target;      // target CFL (worst case)
    static size_t cfl_count;       // check CFL occasionally
    static double cfl_min;         // current CFL minimum
    static double cfl_max;         // current CFL maximum
    static double cfl_tiny;        // smallest cfl so far
    static double time_tiny;       // time at which it occurred

    static double energy_residual;      // to be monitored for steady state
    static Vector3 energy_residual_loc; // location of largest value
    static double mass_residual;
    static Vector3 mass_residual_loc;

    // Filter application parameters.
    static bool   filter_flag = false;
    static double filter_tstart = 0.0;
    static double filter_tend = 0.0;
    static double filter_dt;
    static double filter_next_time;
    static double filter_mu;
    static int filter_npass;
} // end class GlobalConfig
