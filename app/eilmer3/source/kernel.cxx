/// \file kernel.cxx
/// \ingroup eilmer3
/// \brief Core data structures for the 2D code.
///
/// \author PJ
/// \version Jan 2008
/// \version Elmer3 Mar 2008

#include <stdio.h>

#include "../../../lib/util/source/useful.h"

#include "block.hh"
#include "kernel.hh"

//---------------------------------------------------------------------
// The core of the data collection...

// The global control data.
static global_data gd;

global_data * get_global_data_ptr(void) 
{
    return &gd;
}

// The managed gas model lives here.
Gas_model *gmodel;

Gas_model *set_gas_model_ptr(Gas_model *gmptr)
{
    return gmodel = gmptr;
}

Gas_model *get_gas_model_ptr()
{
    return gmodel;
}

// The managed reaction update model lives here.
Reaction_update *rupdate;

int set_reaction_update(std::string file_name)
{
    rupdate = create_Reaction_update(file_name, *(get_gas_model_ptr()));
    if ( rupdate != 0 )
	return SUCCESS;
    else
	return FAILURE;
}

Reaction_update *get_reaction_update_ptr()
{
    return rupdate;
}

// The managed energy exchange update model lives here.
Energy_exchange_update *eeupdate;

int set_energy_exchange_update(std::string file_name)
{
    eeupdate = create_Energy_exchange_update(file_name, *(get_gas_model_ptr()));
    if ( eeupdate != 0 )
	return SUCCESS;
    else
	return FAILURE;
}

Energy_exchange_update *get_energy_exchange_update_ptr()
{
    return eeupdate;
}

// The managed interpolation object lives here.
Thermo_interpolator *tinterp;
int set_thermo_interpolator(std::string name)
{
    tinterp = create_Thermo_interpolator(name);
    if ( tinterp == 0 )
	return FAILURE;
    else
	return SUCCESS;
}

Thermo_interpolator *get_thermo_interpolator_ptr()
{
    return tinterp;
}

// The managed radiation transport model lives here.
RadiationTransportModel *rtm;

int set_radiation_transport_model(std::string file_name)
{
    // 1. Check for an existing managed model
    if ( rtm ) {
	cout << "set_radiation_transport_model()" << endl
	     << "A managed radiation transport model already exists." << endl
	     << "Bailing out!" << endl;
	exit(BAD_INPUT_ERROR);
    }
    
    // 2. create the radiation transport model from the given file name
    rtm = create_radiation_transport_model(file_name);
    if ( rtm != 0 )
	return SUCCESS;
    else
	return FAILURE;
}

RadiationTransportModel *get_radiation_transport_model_ptr()
{
    return rtm;
}

Block * get_block_data_ptr(int i) {
    if ( i < 0 ) return NULL;
    if ( i >= (int)gd.bd.size() ) return NULL;
    return &(gd.bd[i]);
}

void eilmer_finalize( void )
{
    // Clean up the objects created earlier.
    // This will satisfy valgrind, hopefully.
    for ( int ig = 0; ig < gd.n_gas_state; ++ig ) {
	delete gd.gas_state[ig];
    }
    gd.bd.clear();
    gd.pistons.clear();
    gd.heat_zone.clear();
    gd.reaction_zone.clear();
    gd.turbulent_zone.clear();
    gd.my_blocks.clear();
    gd.mpi_rank_for_block.clear();
    delete gmodel;
    delete tinterp;
    if ( get_radiation_flag() )	delete rtm;
    return;
}


//---------------------------------------------------------------------

/// \brief Flag to indicate that we want a verbose startup.
int verbose = 0;  // quiet by default

int set_verbose_flag( int i )
{
    if ( i == 0 ) {
        printf("Brief messages at startup.\n");
	verbose = 0;
    } else {
        printf("Verbose messages at startup.\n");
	verbose = 1;
    }
    return SUCCESS;
}

int get_verbose_flag( void ) { return verbose; }

/*-----------------------------------------------------------------*/

/// \brief Viscous flag =0 for inviscid equations, =1 for viscous terms added.
///
/// Viscous effects are included in the predictor-corrector gas-dynamic update.
int viscous = 0;

/// \brief A factor to scale the viscosity in order to achieve a soft start.
/// 
/// The soft-start for viscous effects may be handy for impulsively-started flows.
double viscous_factor = 1.0;

/// \brief The amount by which to increment the viscous factor during soft-start.
double viscous_factor_increment = 0.01;

/// \brief Viscous upwinding =0 for viscous flux from upwind direction, =1 for average of both
/// directions.
int viscous_upwinding = 0;

/// \brief Diffusion flag =0 for neglecting multicomponent diffusion, =1
//         when considering the diffusion.
//  When the diffusion is calculated is treated as part of the viscous
//  calculation.
int diffusion = 0;           

/// \brief Shock fitting =0 for no shock fitting, =1 for shock fitting
int shock_fitting = 0;

/// \brief Shock fitting decay =0 for no shock fitting decay, =1 for shock fitting decay
int shock_fitting_decay = 0;

/// \brief Adaptive reconstruction =0 for no adaptive reconstruction, =1 for adaptive reconstruction
int adaptive_reconstruction = 0;

/// \brief Axisymmetric flag =0 for 2D planar equations, =1 for 2D axisymmetric. 
///
/// The flow is still 2-dimensional but may includes source terms from the
/// out-of-plane cell interfaces.
int axisymm = 0;

/// \brief Order of time-stepping: =1 for Euler, =2 for predictor-corrector. 
///
/// Euler time-stepping just uses the predictor stage of the time-step functions.
int Torder = 2;

/// \brief Order of reconstruction: =1 for low-order, =2 for high-order. 
///
/// Low order reconstruction uses just the cell-centre data as left- and right-
/// flow properties in the flux calculation.
/// High-order reconstruction adds a correction term to the cell-centre values
/// to approach something like a piecewise-quadratic interpolation between the
/// cell centres.
int Xorder = 2;

/// \brief A factor to scale the heat-addition in order to achieve a soft start.
double heat_factor = 1.0;

/// \brief The amount by which to increment the heat_factor factor during soft-start.
double heat_factor_increment = 0.01;

/// \brief Finite-rate reaction flag: =0 for non-reacting, =1 for reacting.
///
/// Turning on the reactions activates the chemical update function calls.
/// Chemical equilibrium simulations (via Look-Up Table) does not use this
/// chemical update function call.
int reacting = 0;

/// \brief Finite-rate thermal energy exchange: =0 for frozen/equilibrium, =1 for finite-rate
///
/// With this flag on, finite-rate evolution of the vibrational energies (and in turn
/// the total energy) is computed.
int thermal_energy_exchange = 0;

/// \brief Radiation flag: =0 for non-radiating, =1 for radiating.
int radiation = 0;

/// \brief Radiation update frequency: = eg, 1 for every time-step
int radiation_update_frequency = 1;

/// \brief  implicit Flag: =0 normal explicit viscous, 
///         =1 point implicit viscous treatment enabled,
///         =2 fully implicit viscous treatment enabled.
int implicit = 0;

/// \brief any turbulence model active
int turbulence_flag = 0;

/// \brief k-omega turbulence model
int k_omega = 0;

/// \brief Baldwin-Lomax turbulence model
int baldwin_lomax = 0;

/// \brief Turbulent Prandtl number
double Pr_t = 0.89;

/// \brief Turbulent Schmidt number
double Sc_t = 0.75;

/// \brief  apply limiter Flag: =0 no limiting, =1 reconstruction limiter enabled.
int apply_limiter = 1; /* default for users */

/// \brief extrema_clipping_flag: =0 suppress clipping, =1 allow extream clipping
int extrema_clipping_flag = 1; /* default is to clip */

/// \brief suppress_reconstruction_for_species: =0 use x_order reconstruction,
///                                             =1 always use cell average values directly
int suppress_reconstruction_for_species = 0; /* default is to allow reconstruction */

/// \brief A value of 1 will allow the program to complain about bad cells
///        by printing warning messages and cell data.
///        A value of zero suppresses the complaints.
int bad_cell_complain_flag = 1; /* by default, we want to hear about bad cells. */

/// \brief Flag to indicate that we want to apply the strictest CFL check.
///
/// Maybe useful to set this for viscous calculations with large aspect-ratio cells.
int stringent_cfl = 0; // default to direction-by-direction check


/// \brief Flag to indicate that temperature instead of energy is an interpolant
///
/// By default, use internal energy.
bool interpolate_on_temperature_flag = false;


/// \brief Change in normalised velocity to indicate a shock.
///
/// The original default value of -0.05 had been selected to detect the levels of
/// shock compression observed in the inviscid-flow "sod" and "cone20" test cases.
/// It may needed to be tuned for other situations, especially when
/// viscous effects are important.  Hence new higher value should avoid the EFM
/// flux calculator being turned on inappropriately in the boundary layer. 
double compression_tolerance = -0.30;

/// \brief The tolerance to shear when applying the adaptive flux calculator.
///
/// We don't want EFM to be applied in situations of significant shear.
/// The shear value is computed as the tangential-velocity difference across an interface
/// normalised by the local sound speed.
double shear_tolerance = 0.20;

/// \brief For Daryl Bond and Vince Wheatley's MHD additions.
int mhd_flag = 0;

/// \brief A flag for turning on the BGK non-equilibrium gas solver
///
/// This flag can have various settings:
/// BGK_flag == 0: OFF
/// BGK_flag == 1: ON, do not try to import velocity distribution values
/// BGK_flag == 2: ON, read in velocity distribution values from "flow" file
int BGK_flag = 0;

/// \brief The number of velocity buckets for the rarefied gas solver
int velocity_buckets = 0;
std::vector<Vector3> vcoords; // velocity coordinates for rarefied flow
std::vector<double> vweights; // weight for each velocity coordinate

/*------------------------------------------------------------------*/

int set_axisymmetric_flag(int ia)
{
    axisymm = ia;
    if (axisymm == 0) {
        if ( get_verbose_flag() ) printf("Two-dimensional planar flow\n");
    }
    else if (axisymm == 1) {
        if ( get_verbose_flag() ) printf("Axisymmetric flow\n");
    }
    else {
        printf("Invalid axisymmetric flag value: %d\n", axisymm);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_axisymmetric_flag(void)
{
    return axisymm;
}

/*------------------------------------------------------------------*/

int set_shock_fitting_flag(int iw)
{
    shock_fitting = iw;
    if (shock_fitting == 0) {
        if ( get_verbose_flag() ) printf("Turn off shock fititng\n");
    }
    else if (shock_fitting == 1) {
        if ( get_verbose_flag() ) printf("Turn on shock fitting\n");
    }
    else {
        printf("Invalid shock fitting flag value: %d\n", shock_fitting);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_shock_fitting_flag(void)
{
    return shock_fitting;
}

/*------------------------------------------------------------------*/

int set_shock_fitting_decay_flag(int iw)
{
    shock_fitting_decay = iw;
    if (shock_fitting_decay == 0) {
        if ( get_verbose_flag() ) printf("Turn off shock fitting decay\n");
    }
    else if (shock_fitting_decay == 1) {
        if ( get_verbose_flag() ) printf("Turn on shock fitting decay\n");
    }
    else {
        printf("Invalid shock fitting decay flag value: %d\n", shock_fitting_decay);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_shock_fitting_decay_flag(void)
{
    return shock_fitting_decay;
}

/*------------------------------------------------------------------*/

int set_adaptive_reconstruction_flag(int iw)
{
    adaptive_reconstruction = iw;
    if (adaptive_reconstruction == 0) {
        if ( get_verbose_flag() ) printf("Turn off adaptive reconstruction.\n");
    }
    else if (adaptive_reconstruction == 1) {
        if ( get_verbose_flag() ) printf("Turn on adaptive reconstruction.\n");
    }
    else {
        printf("Invalid adaptive reconstruction flag value: %d\n", adaptive_reconstruction);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_adaptive_reconstruction_flag(void)
{
    return adaptive_reconstruction;
}

/*------------------------------------------------------------------*/

int set_viscous_flag(int iv)
{
    viscous = iv;
    if (viscous == 0) {
        if ( get_verbose_flag() ) printf("Turn off viscous terms.\n");
    }
    else if (viscous == 1) {
        if ( get_verbose_flag() ) printf("Turn on viscous terms.\n");
    }
    else {
        printf("Invalid viscous flag value: %d\n", viscous);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_viscous_flag(void)
{
    return viscous;
}

/*------------------------------------------------------------------*/

int set_viscous_upwinding_flag(int iw)
{
    viscous_upwinding = iw;
    if (viscous_upwinding == 0) {
        if ( get_verbose_flag() ) printf("Turn off viscous upwinding.\n");
    }
    else if (viscous_upwinding == 1) {
        if ( get_verbose_flag() ) printf("Turn on viscous upwinding.\n");
    }
    else {
        printf("Invalid viscous upwinding value: %d\n", viscous_upwinding);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_viscous_upwinding_flag(void)
{
    return viscous_upwinding;
}

/*------------------------------------------------------------------*/

/// \brief Set the viscous_factor to a specified value.
double set_viscous_factor( double value )
{
    if ( value > 1.0 ) value = 1.0;
    if ( value < 0.0 ) value = 0.0;
    viscous_factor = value;
    return viscous_factor;
}

/// \brief Get the stored value of viscous_factor.
double get_viscous_factor( void )
{
    return viscous_factor;
}

/// \brief Increment the viscous_factor to a specified value.
double incr_viscous_factor( double value )
{
    viscous_factor += value;
    if ( viscous_factor > 1.0 ) viscous_factor = 1.0;
    if ( viscous_factor < 0.0 ) viscous_factor = 0.0;
    return viscous_factor;
}

/// \brief Set the viscous_factor_increment to a specified value.
double set_viscous_factor_increment( double value )
{
    if ( value > 1.0 ) value = 1.0;
    if ( value < 0.0 ) value = 0.0;
    viscous_factor_increment = value;
    return viscous_factor_increment;
}

/// \brief Set the stored value of the increment.
double get_viscous_factor_increment( void )
{
    return viscous_factor_increment;
}

/*------------------------------------------------------------------*/

int set_diffusion_flag(int id)
{
    diffusion = id;
    if (diffusion == 0) {
        if ( get_verbose_flag() ) printf("Diffusion of species ignored.\n");
    }    
    else if (diffusion == 1) {
        if ( get_verbose_flag() ) printf("Diffusion of species treated as part of viscous terms.\n");
    }
    else {
        printf("Invalid diffusion flag value: %d\n", diffusion);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_diffusion_flag(void)
{
    return diffusion;
}

/*------------------------------------------------------------------*/

int set_Xorder_flag(int ix)
{
    Xorder = ix;
    // if ( get_verbose_flag() ) printf("Xorder=%d\n", Xorder);
    if ( Xorder != 1 && Xorder != 2 ) {
	printf("set_Xorder_flag(): Invalid Xorder flag value: %d\n", Xorder);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_Xorder_flag(void)
{
    return Xorder;
}

/*------------------------------------------------------------------*/

int set_Torder_flag(int it)
{
    Torder = it;
    // if ( get_verbose_flag() ) printf("Torder=%d\n", Torder);
    if ( Torder != 1 && Torder != 2 && Torder != 3 ) {
        printf("set_Torder_flag(): Invalid Torder flag value: %d\n", Torder);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_Torder_flag(void)
{
    return Torder;
}

/*------------------------------------------------------------------*/

/// \brief Set the heat_factor to a specified value.
double set_heat_factor( double value )
{
    if ( value > 1.0 ) value = 1.0;
    if ( value < 0.0 ) value = 0.0;
    heat_factor = value;
    if ( get_verbose_flag() ) printf("set heat_factor=%g\n", heat_factor);
    return heat_factor;
}

/// \brief Get the stored value of heat_factor.
double get_heat_factor( void )
{
    return heat_factor;
}

/// \brief Increment the heat_factor by a specified value.
double incr_heat_factor( double value )
{
    heat_factor += value;
    if ( heat_factor > 1.0 ) heat_factor = 1.0;
    if ( heat_factor < 0.0 ) heat_factor = 0.0;
    return heat_factor;
}

/// \brief Set the heat_factor_increment to a specified value.
double set_heat_factor_increment( double value )
{
    if ( value > 1.0 ) value = 1.0;
    if ( value < 0.0 ) value = 0.0;
    heat_factor_increment = value;
    if ( get_verbose_flag() ) printf("set heat_factor_increment=%g\n", heat_factor_increment);
    return heat_factor_increment;
}

/// \brief Set the stored value of the increment.
double get_heat_factor_increment( void )
{
    return heat_factor_increment;
}

/*------------------------------------------------------------------*/

int set_reacting_flag(int ir)
{
    reacting = ir;
    if (reacting == 0) {
        if ( get_verbose_flag() ) printf("Flow in chemical equilibrium (or frozen)\n");
    }    
    else if (reacting == 1) {
        if ( get_verbose_flag() ) printf("Flow in chemical nonequilibrium: source terms computed\n");
    }
    else {
        printf("Invalid reacting flag value: %d\n", reacting);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_reacting_flag(void)
{
    return reacting;
}

/*------------------------------------------------------------------*/

int set_energy_exchange_flag(int ir)
{
    thermal_energy_exchange = ir;
    if (thermal_energy_exchange == 0) {
        if ( get_verbose_flag() ) printf("Flow in thermal equilibrium\n");
    }    
    else if (thermal_energy_exchange == 1) {
        if ( get_verbose_flag() ) printf("Flow in thermal nonequilibrium: source terms computed\n");
    }    
    else {
        printf("Invalid energy_exchange_flag value: %d\n", thermal_energy_exchange);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_energy_exchange_flag(void)
{
    return thermal_energy_exchange;
}

/*------------------------------------------------------------------*/

int set_radiation_flag(int ir)
{
    radiation = ir;
    if (radiation == 0) {
        if ( get_verbose_flag() ) printf("Flow without radiation\n");
    }    
    else if (radiation == 1) {
        if ( get_verbose_flag() ) printf("Flow with radiation\n");
    }    
    else {
        printf("Invalid radiation flag value: %d\n", radiation);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_radiation_flag(void)                                           
{
    return radiation;
}

int set_radiation_update_frequency(int ruf)
{
    radiation_update_frequency = ruf;
    // if ( get_verbose_flag() ) printf("radiation_update_frequency = %d\n", radiation_update_frequency);
    if ( radiation_update_frequency < 0 ) {
	printf("ERROR: radiation_update_frequency needs to be larger than or equal to 0\n");
	printf("Bailing out!\n");
	exit( BAD_INPUT_ERROR );
    }
    
    return SUCCESS;
}

int get_radiation_update_frequency(void)
{
    return radiation_update_frequency;
}

/*------------------------------------------------------------------*/

int set_implicit_flag(int imf)
{
    implicit = imf;
    // if ( get_verbose_flag() ) printf("set implicit_flag=%d\n", implicit);
    return SUCCESS;
}

int get_implicit_flag(void)
{
    return implicit;
}

/*------------------------------------------------------------------*/

int set_turbulence_flag(int i)
{
    turbulence_flag = i;
    if ( get_verbose_flag() ) printf("set turbulence_flag=%d\n", turbulence_flag);
    return SUCCESS;
}

int get_turbulence_flag(void)
{
    return turbulence_flag;
}

int set_k_omega_flag(int ikw)
{
    k_omega = ikw;
    if ( get_verbose_flag() ) printf("set k_omega_flag=%d\n", k_omega);
    return SUCCESS;
}

int get_k_omega_flag(void)
{
    return k_omega;
}

int set_baldwin_lomax_flag(int ibl)
{
    baldwin_lomax = ibl;
    if ( get_verbose_flag() ) printf("set baldwin_lomax_flag=%d\n", baldwin_lomax);
    return SUCCESS;
}

int get_baldwin_lomax_flag(void)
{
    return baldwin_lomax;
}

double set_turbulence_prandtl_number(double Pr)
{
    Pr_t = Pr;
    if ( get_verbose_flag() ) printf("set turbulence_prandtl_number=%g\n", Pr_t);
    return Pr_t;
}

double get_turbulence_prandtl_number(void)
{
    return Pr_t;
}

double set_turbulence_schmidt_number(double Sc)
{
    Sc_t = Sc;
    if ( get_verbose_flag() ) printf("set turbulence_schmidt_number=%g\n", Pr_t);
    return Sc_t;
}

double get_turbulence_schmidt_number(void)
{
    return Sc_t;
}

/*------------------------------------------------------------------*/

int set_extrema_clipping_flag(int ip)
{
    extrema_clipping_flag = ip;
    if (extrema_clipping_flag == 0) {
        if ( get_verbose_flag() ) printf("Extrema-clipping disabled\n");
    }    
    else if (extrema_clipping_flag == 1) {
        if ( get_verbose_flag() ) printf("Extrema-clipping enabled (default)\n");
    }    
    else {
        printf("Invalid extrema clipping flag value: %d\n", extrema_clipping_flag);
        exit( VALUE_ERROR );
    }
    return SUCCESS;
}

int get_extrema_clipping_flag(void)
{
    return extrema_clipping_flag;
}

/*------------------------------------------------------------------*/

int set_apply_limiter_flag(int ip)
{
    apply_limiter = ip;
    if (apply_limiter == 0) {
        if ( get_verbose_flag() ) printf("Reconstruction limiter disabled\n");
    }    
    else if (apply_limiter == 1) {
        if ( get_verbose_flag() ) printf("Reconstruction limiter applied (default)\n");
    }    
    else {
        printf("Invalid apply_limiter flag value: %d\n", apply_limiter);
        exit( VALUE_ERROR );
    }
    return SUCCESS;
}

int get_apply_limiter_flag(void)
{
    return apply_limiter;
}

/*------------------------------------------------------------------*/

int set_suppress_reconstruction_for_species_flag(int ip)
{
    suppress_reconstruction_for_species = ip;
    if (suppress_reconstruction_for_species == 1) {
        if ( get_verbose_flag() ) printf("High-order reconstruction for species suppressed\n");
    }
    else if (suppress_reconstruction_for_species == 0) {
        if ( get_verbose_flag() ) printf("High-order reconstruction for species allowed (default)\n");
    }
    else {
        printf("Invalid suppress_reconstruction_for_species flag value: %d\n", 
	       suppress_reconstruction_for_species);
        exit( VALUE_ERROR );
    }
    return SUCCESS;
}

int get_suppress_reconstruction_for_species_flag(void)
{
    return suppress_reconstruction_for_species;
}

/*------------------------------------------------------------------*/

int set_bad_cell_complain_flag(int ip)
{
    bad_cell_complain_flag = ip;
    if (bad_cell_complain_flag == 0) {
        if ( get_verbose_flag() ) printf("Will not complain about bad cells.\n");
    }
    else if ( bad_cell_complain_flag == 1 ) {
        if ( get_verbose_flag() ) printf("Will complain about bad cells.\n");
    }
    else {
        printf("Invalid bad_cell_complain_flag value: %d\n", bad_cell_complain_flag);
        exit( VALUE_ERROR );
    }
    return SUCCESS;
}

int get_bad_cell_complain_flag(void)
{
    return bad_cell_complain_flag;
}

/*------------------------------------------------------------------*/

int set_stringent_cfl_flag( int i )
{
    if ( i == 0 ) {
	stringent_cfl = 0;
    } else {
	stringent_cfl = 1;
    }
    // if ( get_verbose_flag() ) printf("set stringent_cfl_flag=%d\n", stringent_cfl);
    return SUCCESS;
}

int get_stringent_cfl_flag( void )
{
    return stringent_cfl;
}

/*------------------------------------------------------------------*/

int set_interpolate_on_temperature_flag( bool flag )
{
    if ( flag == false ) {
	if ( get_verbose_flag() ) printf( "Interpolation uses gas internal energy.\n" );
	interpolate_on_temperature_flag = false;
    } else {
	if ( get_verbose_flag() ) printf( "Interpolation uses gas temperature.\n" );
	interpolate_on_temperature_flag = true;
    }
    return SUCCESS;
}

bool interpolate_on_temperature( void )
{
    return interpolate_on_temperature_flag;
}

//--------------------------------------------------------------------

/// \brief Set the tolerance in relative velocity change for the shock detector.
///
/// This value is expected to be a negative number (for compression)
/// and not too large in magnitude.
/// We have been using a value of -0.05 for years, based on some
/// early experiments with the sod and cone20 test cases, however,
/// the values may need to be tuned for other cases, especially where
/// viscous effects are important.
double set_compression_tolerance( double value )
{
    compression_tolerance = value;
    if ( get_verbose_flag() ) printf("set compression_tolerance=%g\n", compression_tolerance);
    return compression_tolerance;
}

/// \brief Get the stored value for the shock detector.
double get_compression_tolerance( void )
{
    return compression_tolerance;
}


/// \brief Set the tolerance to shear when applying the adaptive flux calculator.
///
/// We don't want EFM to be applied in situations of significant shear.
double set_shear_tolerance( double value )
{
    shear_tolerance = value;
    if ( get_verbose_flag() ) printf("set shear_tolerance=%g\n", shear_tolerance);
    return shear_tolerance;
}

double get_shear_tolerance( void )
{
    return shear_tolerance;
}

//-----------------------------------------------------------------------------

int set_mhd_flag(int i)
{
    mhd_flag = i;
    if ( get_verbose_flag() ) printf("set mhd_flag=%d\n", mhd_flag);
    return mhd_flag;
}

int get_mhd_flag(void)
{
    return mhd_flag;
}

//-----------------------------------------------------------------------------

int set_BGK_flag(int i)
{
    BGK_flag = i;
    if ( get_verbose_flag() ) printf("set BGK_flag=%d\n", BGK_flag);
    return BGK_flag;
}

int get_BGK_flag(void)
{
    return BGK_flag;
}

int set_velocity_buckets(int i)
{
    velocity_buckets = i;

    vcoords.resize(i);
    vweights.resize(i);

    if ( get_verbose_flag() ) printf("set velocity_buckets=%d\n", velocity_buckets);
    return velocity_buckets;
}

int get_velocity_buckets( void )
{
    return velocity_buckets;
}

Vector3 get_vcoord(int i)
{
    return vcoords[i];
}

std::vector<Vector3> *get_vcoords_ptr(void)
{
    return &vcoords;
}

double get_vweight(int i)
{
    return vweights[i];
}

std::vector<double> *get_vweights_ptr(void)
{
    return &vweights;
}
