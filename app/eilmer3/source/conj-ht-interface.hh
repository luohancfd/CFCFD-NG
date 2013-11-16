/// \file conj-ht-interface.hh
/// \file ingroup eilmer3
/// \brief Header file for interface functions with wall conduction model.

#ifndef CONJ_HT_INTERFACE_HH
#define CONJ_HT_INTERFACE_HH

#include <string>

#include "kernel.hh"
#include "../../../lib/geometry2/source/geom.hh"

// Forward declaration, replace with definition
// from Jared's module.
class Wall_model;

Wall_model* initialise_wall_model(std::string fname);
int initialise_wall_node_positions(Wall_model &wm, std::vector<Vector3> &wall_vtxs);
int add_entries_to_wall_vectors(global_data &gd, int nentries);
int gather_wall_fluxes(global_data &gd);
int broadcast_wall_temperatures(global_data &gd);
int update_temperatures_from_fluxes(Wall_model &wm, double dt, const std::vector<double> &q_wall, std::vector<double> &T_wall);

#endif
