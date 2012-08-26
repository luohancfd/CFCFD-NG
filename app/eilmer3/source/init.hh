/** \file init.hh
 * \ingroup eilmer3
 * \brief Prototypes for the initialisation functions.
 */

#ifndef INIT_HH
#define INIT_HH

#include "../../../lib/util/source/config_parser.hh"

int read_config_parameters( const std::string pname, int master );
int read_control_parameters(const string filename, int master, int first_time);
int assign_blocks_to_mpi_rank(const string filename, int master);
CFlowCondition *read_flow_condition_from_ini_dict(ConfigParser &dict, int indx, int master);
int set_block_parameters(int id, ConfigParser &dict, int master);

#endif
