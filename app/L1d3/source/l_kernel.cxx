// l_kernel.cxx

#include "../../../lib/util/source/useful.h"
#include "l_kernel.hh"
#include "l1d.hh"


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

