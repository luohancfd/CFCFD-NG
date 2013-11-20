/// \file conj-ht-interface.hh
/// \file ingroup eilmer3
/// \brief Header file for interface functions with wall conduction model.

#ifndef CONJ_HT_INTERFACE_HH
#define CONJ_HT_INTERFACE_HH

#include <string>

#include "kernel.hh"

int add_entries_to_wall_vectors(global_data &gd, size_t bid, int nentries);
int gather_wall_fluxes(global_data &gd);
int broadcast_wall_temperatures(global_data &gd);


#endif
