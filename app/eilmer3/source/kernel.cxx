/// \file kernel.cxx
/// \ingroup eilmer3
/// \brief Core data structures for the 2D/3D code.
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

std::string revision_string = "PUT_REVISION_STRING_HERE";
std::string get_revision_string() 
{
    return revision_string;
}

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

Block * get_block_data_ptr(size_t i) 
{
    if ( i >= gd.bd.size() ) return NULL;
    return &(gd.bd[i]);
}

void eilmer_finalize( void )
{
    // Clean up the objects created earlier.
    // This will satisfy valgrind, hopefully.
    for ( size_t ig = 0; ig < gd.n_gas_state; ++ig ) {
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
    if ( gd.radiation ) delete rtm;
    if ( gd.conjugate_ht_active ) delete gd.wm;
    return;
}


//---------------------------------------------------------------------

/// \brief Diffusion flag =0 for neglecting multicomponent diffusion, =1
//         when considering the diffusion.
//  When the diffusion is calculated is treated as part of the viscous
//  calculation.
int diffusion = 0;           

/// \brief A factor to scale the diffusion in order to achieve a soft start, separate to viscous effects.
///
/// The soft-start for diffusion effects may be handy for impulsively-started flows.
double diffusion_factor = 1.0;

/// \brief The amount by which to increment the diffusion factor during soft-start.
double diffusion_factor_increment = 0.01;

/// \brief A factor to scale the heat-addition in order to achieve a soft start.
double heat_factor = 1.0;

/// \brief The amount by which to increment the heat_factor factor during soft-start.
double heat_factor_increment = 0.01;

/// \brief  implicit Flag: =0 normal explicit viscous, 
///         =1 point implicit viscous treatment enabled,
///         =2 fully implicit viscous treatment enabled.
int implicit = 0;

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
size_t velocity_buckets = 0;
std::vector<Vector3> vcoords; // velocity coordinates for rarefied flow
std::vector<double> vweights; // weight for each velocity coordinate

/// \brief Electric field work flag: =0 to omit, =1 to include.
int electric_field_work = 0;

/*------------------------------------------------------------------*/

/// \brief Increment the viscous_factor to a specified value.
double incr_viscous_factor( double value )
{
    gd.viscous_factor += value;
    if ( gd.viscous_factor > 1.0 ) gd.viscous_factor = 1.0;
    if ( gd.viscous_factor < 0.0 ) gd.viscous_factor = 0.0;
    return gd.viscous_factor;
}

/*------------------------------------------------------------------*/

int set_diffusion_flag(int id)
{
    diffusion = id;
    if (diffusion == 0) {
        if ( gd.verbose_init_messages ) printf("Diffusion of species ignored.\n");
    }    
    else if (diffusion == 1) {
        if ( gd.verbose_init_messages ) printf("Diffusion of species treated as part of viscous terms.\n");
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

/// \brief Set the diffusion_factor to a specified value.
double set_diffusion_factor( double value )
{
    if ( value > 1.0 ) value = 1.0;
    if ( value < 0.0 ) value = 0.0;
    diffusion_factor = value;
    return diffusion_factor;
}

/// \brief Get the stored value of diffusion_factor.
double get_diffusion_factor( void )
{
    return diffusion_factor;
}

/// \brief Increment the diffusion_factor to a specified value.
double incr_diffusion_factor( double value )
{
    diffusion_factor += value;
    if ( diffusion_factor > 1.0 ) diffusion_factor = 1.0;
    if ( diffusion_factor < 0.0 ) diffusion_factor = 0.0;
    return diffusion_factor;
}

/// \brief Set the diffusion_factor_increment to a specified value.
double set_diffusion_factor_increment( double value )
{
    if ( value > 1.0 ) value = 1.0;
    if ( value < 0.0 ) value = 0.0;
    diffusion_factor_increment = value;
    return diffusion_factor_increment;
}

/// \brief Set the stored value of the increment.
double get_diffusion_factor_increment( void )
{
    return diffusion_factor_increment;
}

/*------------------------------------------------------------------*/

/// \brief Set the heat_factor to a specified value.
double set_heat_factor( double value )
{
    if ( value > 1.0 ) value = 1.0;
    if ( value < 0.0 ) value = 0.0;
    heat_factor = value;
    if ( gd.verbose_init_messages ) printf("set heat_factor=%g\n", heat_factor);
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
    if ( gd.verbose_init_messages ) printf("set heat_factor_increment=%g\n", heat_factor_increment);
    return heat_factor_increment;
}

/// \brief Set the stored value of the increment.
double get_heat_factor_increment( void )
{
    return heat_factor_increment;
}
/*------------------------------------------------------------------*/

int set_implicit_flag(int imf)
{
    implicit = imf;
    // if ( gd.verbose_init_messages ) printf("set implicit_flag=%d\n", implicit);
    return SUCCESS;
}

int get_implicit_flag(void)
{
    return implicit;
}

/*------------------------------------------------------------------*/

int set_mhd_flag(int i)
{
    mhd_flag = i;
    if ( gd.verbose_init_messages ) printf("set mhd_flag=%d\n", mhd_flag);
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
    if ( gd.verbose_init_messages ) printf("set BGK_flag=%d\n", BGK_flag);
    return BGK_flag;
}

int get_BGK_flag(void)
{
    return BGK_flag;
}

size_t set_velocity_buckets(size_t i)
{
    velocity_buckets = i;

    vcoords.resize(i);
    vweights.resize(i);

    if ( gd.verbose_init_messages )
	printf("set velocity_buckets=%d\n", static_cast<int>(velocity_buckets));
    return velocity_buckets;
}

size_t get_velocity_buckets( void )
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

//---------------------------------------------------------------------

std::string get_name_of_turbulence_model(turbulence_model_t my_model)
{
    switch ( my_model ) {
    case TM_NONE: return "none";
    case TM_BALDWIN_LOMAX: return "baldwin-lomax";
    case TM_K_OMEGA: return "k-omega";
    case TM_SPALART_ALLMARAS: return "spalart-allmaras";
    default: return "unknown";
    }
} // end get_name_of_turbulence_model()

/*------------------------------------------------------------------*/

int set_electric_field_work_flag(int iefw)
{
    electric_field_work = iefw;
    if (electric_field_work == 0) {
        if ( gd.verbose_init_messages ) printf("Flow without electric field work\n");
    }
    else if (electric_field_work == 1) {
        if ( gd.verbose_init_messages ) printf("Flow with electric field work\n");
    }
    else {
        printf("Invalid electric field work flag value: %d\n", electric_field_work);
        exit(VALUE_ERROR);
    }
    return SUCCESS;
}

int get_electric_field_work_flag(void)
{
    return electric_field_work;
}
