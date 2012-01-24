// l_kernel.hh

#ifndef L_KERNEL_HH
#define L_KERNEL_HH

#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/kinetics/reaction-update.hh"
#include "../../../lib/gas/kinetics/energy-exchange-update.hh"

// Managed gas models
Gas_model *set_gas_model_ptr(Gas_model *gmptr);
Gas_model *get_gas_model_ptr();
int set_reaction_update(std::string file_name);
Reaction_update *get_reaction_update_ptr();
int set_energy_exchange_update( std::string file_name );
Energy_exchange_update *get_energy_exchange_update_ptr();

#endif
