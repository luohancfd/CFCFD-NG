/// \file conj-ht-interface.hh
/// \file ingroup eilmer3
/// \brief Header file for interface functions with wall conduction model.

#ifndef CONJ_HT_INTERFACE_HH
#define CONJ_HT_INTERFACE_HH

#include <vector>

#include "../../../lib/geometry2/source/geom.hh"
#include "kernel.hh"

class Conjugate_HT_Interface {
public:
    // We intend only one object of this type with all data known
    // at initialisation time, so delete some useless constructors
    Conjugate_HT_Interface() = delete;
    Conjugate_HT_Interface(const Conjugate_HT_Interface &c) = delete;
    Conjugate_HT_Interface& operator=(const Conjugate_HT_Interface &c) = delete;

    // Our only useful constructor
    Conjugate_HT_Interface(std::vector<Vector3> &gas_cell_pos, std::vector<Vector3> &solid_cell_pos,
			   std::vector<Vector3> &iface_pos, std::vector<Vector3> &iface_normal);

    // Our service
    int compute_flux(std::vector<double> &T_gas, std::vector<double> &k_gas,
		     std::vector<double> &T_solid, std::vector<double> &k_solid, 
		     std::vector<double> &q_wall, std::vector<double> &T_wall);

private:
    std::vector<Vector3> gas_cell_pos_;
    std::vector<Vector3> solid_cell_pos_;
    std::vector<Vector3> iface_pos_;
    std::vector<Vector3> iface_normal_;
    std::vector<double> T_iface;
};

int add_entries_to_wall_vectors(global_data &gd, size_t bid, int nentries);
int gather_near_wall_gas_data(global_data &gd);
int broadcast_wall_values(global_data &gd);


#endif
