/// \file conj-ht-interface.cxx
/// \file ingroup eilmer3

#ifdef _MPI
   #include <mpi.h>
#endif
#include "bc.hh"
#include "kernel.hh"
#include "conj-ht-interface.hh"
#include <cmath>
#include <iomanip>

using namespace std;

Conjugate_HT_Interface::
Conjugate_HT_Interface(vector<Vector3> &gas_cell_pos, vector<Vector3> &solid_cell_pos,
		       vector<Vector3> &iface_pos, vector<Vector3> &iface_normal)
    : gas_cell_pos_(gas_cell_pos), solid_cell_pos_(solid_cell_pos),
      iface_pos_(iface_pos), iface_normal_(iface_normal), T_iface(gas_cell_pos.size(), 0.0) {}

int
Conjugate_HT_Interface::
compute_flux(vector<double> &T_gas, vector<double> &k_gas, 
	     vector<double> &T_solid, vector<double> &k_solid, 
	     vector<double> &q_wall, vector<double> &T_wall)
{

    // Inputs:
    // T_gas -- array of temperatures at near wall cells in gas
    // k_gas -- array of thermal conductivities at near wall cells in gas
    // T_solid -- array of temperature at near wall cells in solid
    // k_solid -- array of thermal conductivities at near wall cells in solid
    //
    // Outputs:
    // q_wall -- compute fluxes at gas/solid interfaces
    // T_wall -- update corresponding interface temperatures also
    //
    // Extra data to use (set at start of simulation):
    // gas_cell_pos_ -- cell centre positions for gas cells
    // solid_cell_pos_ -- cell centre positions for solid cells
    // iface_pos_ -- interface positions
    // iface_normal_ -- interface normals
    //
    // Working arrays:
    // T_iface -- a storage area for interface temperature
    // POSSIBLY ADD dx_gas, dy_gas, dx_solid, dy_solid

    for ( size_t i = 0; i < T_gas.size(); i++) {
	double dx_gas, dx_solid, dy_gas, dy_solid, dn_gas, dn_solid, A, B, cosa, cosb; //cosa, sina, 
	//	std::cout << "i = " << i << std::endl;
	dx_gas = iface_pos_[i].x - gas_cell_pos_[i].x;
	dy_gas = iface_pos_[i].y - gas_cell_pos_[i].y;
	dx_solid = solid_cell_pos_[i].x - iface_pos_[i].x;
	dy_solid = solid_cell_pos_[i].y - iface_pos_[i].y;
	
	cosa = iface_normal_[i].x; 				
	cosb = iface_normal_[i].y; 				
	dn_gas = abs(cosa*dx_gas + cosb*dy_gas); 
	dn_solid = abs(cosa*dx_solid + cosb*dy_solid);
	//	std::cout << "cosa = " << cosa << " cosb = " << cosb << std::endl;
	A = k_gas[i]/dn_gas; 		//isotropic
	B = k_solid[i]/dn_solid; 	//isotropic;

	T_iface[i] = (T_solid[i]*B + T_gas[i]*A)/(A+B);
	T_wall[i] = T_iface[i]; //Added by Justin Beri 02-08-2014
	q_wall[i] = - (A * (T_iface[i] - T_gas[i]));
	//	std::cout << std::setprecision(12) << "i = " << i << ",  q = " << q_wall[i] << ",  T_int = " << T_iface[i] << ", T_gas = " << T_gas[i] << ", k_gas = " << k_gas[i] << ", dn_gas = "<< dn_gas << std::endl;
    }
    return SUCCESS;
}

int add_entries_to_wall_vectors(global_data &gd, size_t bid, int nentries)
{
    gd.T_gas_near_wall.insert(gd.T_gas_near_wall.end(), nentries, 0.0);
    gd.k_gas_near_wall.insert(gd.k_gas_near_wall.end(), nentries, 0.0);
    gd.T_solid_near_wall.insert(gd.T_solid_near_wall.end(), nentries, 0.0);
    gd.k_solid_near_wall.insert(gd.k_solid_near_wall.end(), nentries, 0.0);
    gd.q_wall.insert(gd.q_wall.end(), nentries, 0.0);
    gd.T_wall.insert(gd.T_wall.end(), nentries, 0.0); //Added by Justin Beri 02-08-2014
    gd.recvcounts.push_back(nentries);
    if ( bid == 0 ) {
	gd.displs.push_back(0);
    }
    else {
	// for all other ranks
	gd.displs.push_back(gd.displs.back()+nentries);
    }
    return 0;
}

int gather_near_wall_gas_data(global_data &gd)
{
#   ifdef _MPI    
    MPI_Barrier(MPI_COMM_WORLD);
    int rank = gd.my_mpi_rank;
    Block& bdp = gd.bd[rank];
    // We only need correct values in T vector on the
    // blocks that have a conjugate heat transfer bc.
    if ( bdp.bcp[NORTH]->type_code == CONJUGATE_HT ) {
	int j = bdp.jmax;
	int T_idx = gd.displs[rank]; // starting index in T vector
	for ( size_t i = bdp.imin; i <= bdp.imax; ++i ) {
	    FV_Cell &cell = *(bdp.get_cell(i, j));
	    gd.T_gas_near_wall[T_idx] = cell.fs->gas->T[0];
	    gd.k_gas_near_wall[T_idx] = cell.fs->gas->k[0];
	    ++T_idx;
	}
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Now perform our MPI magic to gather the result on rank 0
    double *pT_local = &(gd.T_gas_near_wall[gd.displs[rank]]);
    int nT_local = gd.recvcounts[rank];
    double *pT_global = &(gd.T_gas_near_wall[0]);
    int *recvcounts = &(gd.recvcounts[0]);
    int *displs = &(gd.displs[0]);
    MPI_Gatherv(pT_local, nT_local, MPI_DOUBLE, pT_global, recvcounts, displs,
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    double *pk_local = &(gd.k_gas_near_wall[gd.displs[rank]]);
    int nk_local = gd.recvcounts[rank];
    double *pk_global = &(gd.k_gas_near_wall[0]);
    MPI_Gatherv(pk_local, nk_local, MPI_DOUBLE, pk_global, recvcounts, displs,
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#   else
    // In shared memory, loop over blocks gathering near wall temperatures.
    for ( Block *bdp : gd.my_blocks ) {
	if ( bdp->bcp[NORTH]->type_code == CONJUGATE_HT ) {
	    int j = bdp->jmax;
	    int T_idx = gd.displs[bdp->id]; // starting index in T vector
	    for ( size_t i = bdp->imin; i <= bdp->imax; ++i ) {
		FV_Cell &cell = *(bdp->get_cell(i, j));
		gd.T_gas_near_wall[T_idx] = cell.fs->gas->T[0];
		gd.k_gas_near_wall[T_idx] = cell.fs->gas->k[0];
		++T_idx;
	    }
	}
    }
#   endif    
    return SUCCESS;
}

int broadcast_wall_values(global_data &gd)
{
#   ifdef _MPI
    // Broadcast wall fluxes
    double *pq = &(gd.q_wall[0]);
    int nq = gd.q_wall.size();
    MPI_Bcast(pq, nq, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    // Broadcast wall temperatures
    pq = &(gd.T_wall[0]);
    nq = gd.q_wall.size();
    MPI_Bcast(pq, nq, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#   endif
    return SUCCESS;
}


				    

