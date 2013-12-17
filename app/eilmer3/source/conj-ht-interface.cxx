/// \file conj-ht-interface.cxx
/// \file ingroup eilmer3

#ifdef _MPI
   #include <mpi.h>
#endif
#include "bc.hh"
#include "kernel.hh"
#include "conj-ht-interface.hh"

int add_entries_to_wall_vectors(global_data &gd, size_t bid, int nentries)
{
    gd.T_wall.insert(gd.T_wall.end(), nentries, 0.0);
    gd.q_wall.insert(gd.q_wall.end(), nentries, 0.0);
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

int gather_wall_fluxes(global_data &gd)
{
#   ifdef _MPI    
    MPI_Barrier(MPI_COMM_WORLD);
    int rank = gd.my_mpi_rank;
    Block& bdp = gd.bd[rank];
    // We only need correct values in q vector on the
    // blocks that have a conjugate heat transfer bc.
    if ( bdp.bcp[NORTH]->type_code == CONJUGATE_HT ) {
	int j = bdp.jmax + 1;
	int iq = gd.displs[rank]; // starting index in q vector
	for ( size_t i = bdp.imin; i <= bdp.imax; ++i ) {
	    FV_Interface &iface = *(bdp.get_ifj(i, j));
	    gd.q_wall[iq] = iface.F->total_energy;
	    ++iq;
	}
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Now perform our MPI magic to gather the result on rank 0
    double *pq_local = &(gd.q_wall[gd.displs[rank]]);
    int nq_local = gd.recvcounts[rank];
    double *pq_global = &(gd.q_wall[0]);
    int *recvcounts = &(gd.recvcounts[0]);
    int *displs = &(gd.displs[0]);
    MPI_Gatherv(pq_local, nq_local, MPI_DOUBLE, pq_global, recvcounts, displs,
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#   else
    // In shared memory, loop over blocks gathering fluxes.
    for ( Block *bdp : gd.my_blocks ) {
	if ( bdp->bcp[NORTH]->type_code == CONJUGATE_HT ) {
	    int j = bdp->jmax + 1;
	    int iq = gd.displs[bdp->id]; // starting index in q vector
	    for ( size_t i = bdp->imin; i <= bdp->imax; ++i ) {
		FV_Interface &iface = *(bdp->get_ifj(i, j));
		gd.q_wall[iq] = iface.F->total_energy;
		++iq;
	    }
	}
    }
#   endif    
    return SUCCESS;
}

int broadcast_wall_temperatures(global_data &gd)
{
#   ifdef _MPI
    double *pT = &(gd.T_wall[0]);
    int nT = gd.T_wall.size();
    MPI_Bcast(pT, nT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
#   endif
    return SUCCESS;
}


				    

