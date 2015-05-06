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
import std.string;

import gas;
import geom;
import fvcore;
import solidfvcore;

// Symbolic names for turbulence models.
enum TurbulenceModel { none, baldwin_lomax, k_omega, spalart_allmaras }
string turbulence_model_name(TurbulenceModel i)
{
    final switch (i) {
    case TurbulenceModel.none: return "none";
    case TurbulenceModel.baldwin_lomax: return "baldwin_lomax";
    case TurbulenceModel.k_omega: return "k_omega";
    case TurbulenceModel.spalart_allmaras: return "spalart_allmaras";
    }
} // end turbulence_model_name()

TurbulenceModel turbulence_model_from_name(string name)
{
    switch (name) {
    case "none": return TurbulenceModel.none;
    case "baldwin_lomax": return TurbulenceModel.baldwin_lomax;
    case "k_omega": return TurbulenceModel.k_omega;
    case "spalart_allmaras": return TurbulenceModel.spalart_allmaras;
    default: return TurbulenceModel.none;
    }
} // end turbulence_model_from_name()

class BlockZone {
    shared double p0x, p0y, p0z;
    shared double p1x, p1y, p1z;
    this(in Vector3 p0, in Vector3 p1) {
	p0x = p0.x; p0y = p0.y; p0z = p0.z;
	p1x = p1.x; p1y = p1.y; p1z = p1.z;  
    }
    override string toString() {
	return format("BlockZone([%g,%g,%g], [%g,%g,%g])", p0x,p0y,p0z,p1x,p1y,p1z);
    }
    shared bool is_inside(in Vector3 p, int dimensions) {
	if ( p.x >= p0x && p.x <= p1x &&
	     p.y >= p0y && p.y <= p1y ) {
	    if ( dimensions == 2 ) {
		return true;
	    } else if ( p.z >= p0z && p.z <= p1z ) {
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
	return format("IgnitionZone([%g,%g,%g], [%g,%g,%g],%g)", p0x,p0y,p0z,p1x,p1y,p1z,_Tig);
    }
}


final class GlobalConfig {
    shared static string base_file_name = "job"; // Change this to suit at run time.
    shared static string title = "Eilmer4 simulation"; // Change this to suit at run time.
    shared static string gas_model_file = "gas-model.lua";
    // Note that the following reference to the GasModel is NOT shared.
    static GasModel gmodel_master; // to provide GasModel access to the Lua domain
    // Presumably, we won't be accessing this particular gas model from the 
    // individual block computations, so that parallel computations for the blocks
    // don't trip over each other.

    shared static int nBlocks; // Number of blocks in the overall simulation.
    shared static int nSolidBlocks; // Number of solid blocks in the overall simulation.
    shared static int dimensions = 2; // default is 2, other valid option is 3
    shared static bool axisymmetric = false;

    // Parameters controlling convective update
    //
    shared static GasdynamicUpdate gasdynamic_update_scheme = GasdynamicUpdate.pc;

    // We might update some properties in with the main convective-terms
    // time-stepping function or we might choose to update them separately, 
    // like the chemistry update.
    shared static bool separate_update_for_viscous_terms = false;
    shared static bool separate_update_for_k_omega_source = false;

    /// When decoding the array of conserved quantities, 
    /// the temperature or the density may try to go negative.  
    /// If it does and adjust_invalid_cell_data == true, the cell data
    /// is adjusted to make it reasonable.
    shared static bool adjust_invalid_cell_data = false;
    // The maximum number of bad cells (per block) 
    // which will be tolerated without complaint.
    shared static int max_invalid_cells = 0;  
    // If an attempt at a time step fails because of invalid cells,
    // the time step is re-attempted with a smaller time step.
    // This reduction factor is somewhat arbitrary and can now be set
    // by the user's imput script.
    // A factor of 0.5 would seem to be not enough but a factor of
    // 0.1 would seem too expensive.  We have settled on a default of 0.2.
    shared static double dt_reduction_factor = 0.2; 

    // Low order reconstruction (1) uses just the cell-centre data as left- and right-
    // flow properties in the flux calculation.
    // High-order reconstruction (2) adds a correction term to the cell-centre values
    // to approach something like a piecewise-quadratic interpolation between the
    // cell centres.
    shared static int interpolation_order = 2; 

    // Default flow-data reconstruction includes interpolation of density 
    // and internal energy.  Other options for the thermodunamic properties
    // to be interpolated are pressure+temperature, density+temperature and
    // density+pressure.
    shared static InterpolateOption thermo_interpolator = InterpolateOption.rhoe;
    shared static bool apply_limiter = true;
    shared static bool extrema_clipping = true;
    shared static bool interpolate_in_local_frame = true;

    // Default flux calculator is the adaptive mix of ausmdv and efm.
    shared static FluxCalculator flux_calculator = FluxCalculator.adaptive;

    // Set the tolerance to shear when applying the adaptive flux calculator.
    // We don't want EFM to be applied in situations of significant shear.
    // The shear value is computed as the tangential-velocity difference across an interface
    // normalised by the local sound speed.
    shared static double shear_tolerance = 0.20;

    // Reference free-stream Mach number, for use in the ausm_plus_up flux calculator.
    // Choose a value for M_inf that is good for low Mach numbers.
    // To be strictly correct, we should set this at run time
    // if an M_inf value is easily defined.
    shared static double M_inf = 0.01;

    // Set the tolerance in relative velocity change for the shock detector.
    // This value is expected to be a negative number (for compression)
    // and not too large in magnitude.
    // We have been using a value of -0.05 for years, based on some
    // early experiments with the sod and cone20 test cases, however,
    // the values may need to be tuned for other cases, especially where
    // viscous effects are important.
    shared static double compression_tolerance = -0.30;

    shared static bool moving_grid = false;
    shared static bool write_vertex_velocities = false;

    // Parameters controlling viscous/molecular transport
    //
    shared static bool viscous = false; 
    // If true, viscous effects are included in the gas-dynamic update.
    // A factor to scale the viscosity in order to achieve a soft start. 
    // The soft-start for viscous effects may be handy for impulsively-started flows.
    // A value of 1.0 means that the viscous effects are fully applied.
    shared static double viscous_factor = 1.0;
    // The amount by which to increment the viscous factor during soft-start.
    shared static double viscous_factor_increment = 0.01;
    shared static double viscous_delay = 0.0;

    // When the diffusion is calculated is treated as part of the viscous calculation:
    //   false for neglecting multicomponent diffusion, 
    //   true when considering the diffusion 
    shared static bool diffusion = false; 
    // A factor to scale the diffusion in order to achieve a soft start.
    // The soft-start for diffusion effects may be handy for impulsively-started flows.
    // Note that this is separate to viscous effects.
    shared static double diffusion_factor = 1.0;
    // The amount by which to increment the diffusion factor during soft-start.
    shared static double diffusion_factor_increment = 0.01;
    shared static double diffusion_time_delay = 0.0;
    // The Lewis number when using the constant Lewis number diffusion model
    shared static double diffusion_lewis = 1.0;
    // The Schmidt number when using the constant Schmidt number diffusion model
    shared static double diffusion_schmidt = 0.7;

    shared static TurbulenceModel turbulence_model = TurbulenceModel.none;
    shared static double turbulence_prandtl_number = 0.89;
    shared static double turbulence_schmidt_number = 0.75;
    shared static double max_mu_t_factor = 300.0;
    shared static double transient_mu_t_factor = 1.0;
    shared static BlockZone[] turbulent_zones;

    // Parameters controlling thermochemistry
    //
    // Turning on the reactions activates the chemical update function calls.
    // Chemical equilibrium simulations (via Look-Up Table) does not use this
    // chemical update function call.
    shared static bool reacting = false;
    shared static string reactions_file = "chemistry.lua";

    // With this flag on, finite-rate evolution of the vibrational energies 
    // (and in turn the total energy) is computed.
    shared static bool thermal_energy_exchange = false;

    shared static bool radiation = false;
    shared static int radiation_update_frequency; // = 1 for every time-step
    shared static bool radiation_scaling = false;

    shared static bool electric_field_work;

    // For Daryl Bond and Vince Wheatley's MHD additions.
    //
    shared static bool MHD;

    // A flag for turning on the BGK non-equilibrium gas solver:
    //   BGK == 0: OFF
    //   BGK == 1: ON, do not try to import velocity distribution values
    //   BGK == 2: ON, read in velocity distribution values from "flow" file
    shared static int BGK = 0;

    // Parameters controlling other simulation options
    //
    shared static double ignition_time_start = 0.0;
    shared static double ignition_time_stop = 0.0;
    shared static IgnitionZone[] ignition_zones;
    shared static bool ignition_zone_active = false;

    shared static double reaction_time_start = 0.0;
    shared static double T_frozen; // temperature (in K) below which reactions are frozen
    shared static double T_frozen_energy; // temperature (in K) below which energy exchanges are skipped
    shared static BlockZone[] reaction_zones;

    shared static int max_step = 10;       // iteration limit
    shared static int t_level;             // time level within update
    shared static int halt_now = 0;        // flag for premature halt
    shared static bool halt_on_large_flow_change = false; 
    // Set to true to halt simulation when any
    // monitor point sees a large flow change.
    shared static double tolerance_in_T;   // Temperature change for the flow change.
    shared static int print_count = 20; // Number of steps between writing messages to console.
    shared static int control_count = 10; // Number of steps between rereading .control file.

    shared static int verbosity_level = 1;
    // Messages have a hierarchy:  // [TODO] we are not really abiding by this.
    // 0 : only error messages will be omitted
    // 1 : emit messages that are useful for a long-running job (default)
    // 2 : plus verbose init messages
    // 3 : plus verbose boundary condition messages

    shared static double max_time = 1.0e-3; // final solution time, s, set by user
    shared static double dt_init = 1.0e-6; // initial time step set by user
    shared static double dt_max = 1.0e-3; // Maximum allowable time-step, after all other considerations.
    shared static double cfl_value = 0.5; // target CFL number (worst case) set by user
    shared static bool stringent_cfl = false; 
    // If true, assume the worst with respect to cell geometry and wave speed.
    shared static size_t cfl_count = 10;  // steps between checking time step size
    shared static bool fixed_time_step = false; // set true to fix dt_allow

    shared static size_t write_at_step = 0; // update step at which to write a solution, 0=don't do it
    shared static double dt_plot = 1.0e-3; // interval for writing soln
    shared static double dt_history = 1.0e-3; // interval for writing sample

    shared static double energy_residual;      // to be monitored for steady state
    shared static Vector3 energy_residual_loc; // location of largest value
    shared static double mass_residual;
    shared static Vector3 mass_residual_loc;

    // Filter application parameters.
    shared static bool   filter_flag = false;
    shared static double filter_tstart = 0.0;
    shared static double filter_tend = 0.0;
    shared static double filter_dt;
    shared static double filter_next_time;
    shared static double filter_mu;
    shared static int filter_npass;

    // Indicate presence of user-defined source terms
    shared static string udf_source_terms_file = "dummy-source-terms.txt";
    shared static bool udf_source_terms = false;

    // Parameters related to the solid domain and conjugate coupling
    shared static SolidDomainUpdate solidDomainUpdateScheme = SolidDomainUpdate.pc;
    shared static bool udfSolidSourceTerms = false;

} // end class GlobalConfig
