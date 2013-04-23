/** \file init.hh
 * \ingroup eilmer3
 * \brief Prototypes for the initialisation functions.
 */

#ifndef INIT_HH
#define INIT_HH

#include "../../../lib/util/source/config_parser.hh"

int init_available_schemes_map();
int init_available_calculators_map();
int init_available_turbulence_models_map();
int read_config_parameters(const std::string pname, bool master);
int read_control_parameters(const string filename, bool master, bool first_time);
int assign_blocks_to_mpi_rank(const string filename, bool master);
CFlowCondition *read_flow_condition_from_ini_dict(ConfigParser &dict, size_t indx, bool master);
int set_block_parameters(size_t id, ConfigParser &dict, bool master);

#endif
