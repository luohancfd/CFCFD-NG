#ifndef SOLID_BLOCK_HH
#define SOLID_BLOCK_HH

#include <string>
#include <vector>
#include "solid_cell.hh"

class SolidBoundaryCondition; //Required to initialise

//Block to contain FV_Cells etc
class SolidBlock {
public: //public access permission

    int id; // block number

    int nni, nnj, nnk;  //number of nodes in i,j,k directions
    int nci,ncj ,nck; //number of cells nn-1.

    double rho,k,cp, T_init;
    double k11,k12,k22;
    int anisotropic_flag;
    int time_varying_terms_flag;
    std::string time_varying_terms_dir;
    int time_varying_iteration;

    std::vector< std::vector<Solid_FV_Cell> > block_cells; // cells in block
    std::vector< std::vector<Solid_FV_Interface> > block_iface_i; //ifaces with eastward normals
    std::vector< std::vector<Solid_FV_Interface> > block_iface_j; //ifaces with northward normals
    std::vector< std::vector<Solid_FV_Vertex> > block_vertices;

    std::vector< std::vector<Solid_FV_Interface> > block_secondary_iface_i; //ifaces with eastward normals
    std::vector< std::vector<Solid_FV_Interface> > block_secondary_iface_j; //ifaces with northward normals

    std::vector<SolidBoundaryCondition *> block_bc; // [South, East, North, West]
    std::vector<int> e3connection_type_flag; //This is falg to determine flux or temp for connection.

    SolidBlock();
    void allocate_memory();
    void assign_cells_to_block();
    void assign_ifaces_to_block();
    void assign_ifaces_to_cells();
    void assign_internal_secondary_geometry_to_vertices();
    void set_BC(std::string, std::string);
    void assign_secondary_ifaces();
    void assign_properties_to_cells();
    void initialise_cells(double T_init);
    
    void update_vertex_derivatives();
    void update_interface_fluxes();
    void update_internal_secondary_interface_temperatures();
    void update_boundary_secondary_interface_temperatures();
    void OLD_update_boundary_conditions();
    void update_boundary_conditions();
    void update_cell_energy_derivative();
    void time_update(double &dt, int update_scheme);
    
};

#endif // SOLID_BLOCK_HH
