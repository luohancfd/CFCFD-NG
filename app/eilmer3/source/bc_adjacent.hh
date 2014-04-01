// bc_adjacent.hh

#ifndef BC_ADJACENT_HH
#define BC_ADJACENT_HH

#include <vector>
#include "bc.hh"

class AdjacentBC : public BoundaryCondition {
public:
    // Note that we don't have the neighbour_* data here because they are
    // required by the exchange functions that get only pointers to the base class.
    bool reorient_vector_quantities;
    std::vector<double> Rmatrix;

public:
    AdjacentBC(Block *bdp, int which_boundary, 
	       int other_block, int other_face, int neighbour_orientation_,
	       bool reorient_vector_quantities_, vector<double>& Rmatrix_);
    AdjacentBC(const AdjacentBC &bc);
    AdjacentBC();
    AdjacentBC & operator=(const AdjacentBC &bc);
    virtual ~AdjacentBC();
    virtual void print_info(std::string lead_in);
    virtual int apply_convective(double t);
    // default apply_viscous() (does nothing)
};

// Helper functions
void apply_matrix_transform(const std::vector<double>& Rmatrix, 
			    const std::vector<double>& oldv,
			    std::vector<double> &newv);
void reorient_vector_quantities_in_cell(FV_Cell *c, const std::vector<double>& Rmatrix); 

#endif
