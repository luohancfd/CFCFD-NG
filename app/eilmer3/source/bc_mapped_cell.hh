// bc_mapped_cell.hh

#ifndef BC_MAPPED_CELL_HH
#define BC_MAPPED_CELL_HH

#include <vector>
#include "bc.hh"

class MappedCellBC : public BoundaryCondition {
public:
    // incoming- and outgoing-mapped-cell lists are in 
    // the BoundaryCondition class definition so that they're
    // accessible to the MPI exchange functions that only have
    // pointers to the base BoundaryCondition class.

    // This following data structure stores 4 integers,
    // specifying the mapped cell block and ijk-indices,
    // for each ghost cell associated with the boundary.
    std::vector<std::vector<int>> mapped_cells;
    bool reorient_vector_quantities;
    std::vector<double> Rmatrix;

public:
    MappedCellBC(Block *bdp, int which_boundary, 
		 std::string filename,
		 bool _reorient_vector_quantities, vector<double>& _Rmatrix);
    MappedCellBC(const MappedCellBC &bc);
    MappedCellBC();
    MappedCellBC & operator=(const MappedCellBC &bc);
    virtual ~MappedCellBC();
    virtual void print_info(std::string lead_in);
    virtual int apply_convective(double t);
    // default apply_viscous() (does nothing)
};

#endif
