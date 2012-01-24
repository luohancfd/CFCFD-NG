/** \file flex_cell.hh
 * \ingroup basic_gas_dyn
 * \brief Header file for the flexible cell centre data.
 * 
 * \author Brendan O'Flaherty
 * \version Feb-09

 ********* FIX-ME ********
 * July 2010 -- Cut most of the code out so we don't have to reworkit immediately.
 *              Need to update to new data storage arrangement.
 */

#ifndef _FLEX_CELL_HH
#define _FLEX_CELL_HH

#include "cell.hh"

// --------------------------------------------------------------

struct flex_cell_center : FV_Cell
{
    // Some values unique to flex-cells
    Vector3* n; // pointer to surface unit vector
    Vector3 vel; // velocity vector, m/s

    // Extensive conserved variables
    double mass;  // mass
    Vector3 mv;  // momentum
    double mE;  // total energy
    std::vector<double> mf;  // mass of species
    std::vector<double> menergies;  // individual energies

    // Old extensive conserved variables
    double mass_old;  // mass
    Vector3 mv_old;  // momentum
    double mE_old;  // energy
    std::vector<double> mf_old;  // partial mass of species
    std::vector<double> menergies_old;  // species energies

    // Derivatives of extensive conserved variables
    double DmDt[NL];  // mass       
    Vector3 DmvDt[NL];  // momentum                       
    double DmEDt[NL];  // total energy                   
    std::vector<double> DmfDt[NL];  // species mass             
    std::vector<double> DmenergiesDt[NL];  // species energies             
};

// --------------------------------------------------------------

int read_solution_for_flex_cell(flex_cell_center* fcell, FILE *infile); 
int write_solution_for_flex_cell(flex_cell_center* fcell, FILE *outfile);  
int set_array_sizes_for_flex_cell(flex_cell_center &c, size_t nsp, size_t nvib);  

int record_conserved_for_flex_cell(flex_cell_center* fcell);
int restore_conserved_for_flex_cell(flex_cell_center* fcell);

int encode_intensive_to_extensive(flex_cell_center* fcell);
int decode_extensive_to_intensive(flex_cell_center* fcell); 

int copy_shadowed_cell_data_to_flex_cell(FV_Cell* cell, flex_cell_center* fcell, 
					 int type_of_copy=COPY_FLOW_STATE); 
int copy_flex_cell_data_to_shadowed_cell(flex_cell_center* fcell, FV_Cell* cell,
					 int type_of_copy=COPY_FLOW_STATE);

int add_shadowed_cell_and_flex_cell_data(FV_Cell* src, flex_cell_center* dest); 

int calculate_flex_cell_volume(flex_cell_center* fcell); 

int calculate_flex_cell_length(flex_cell_center* fcell); 

int time_derivatives_for_flex_cell(flex_cell_center *fcell, int time_level, int dimensions=2); 

int calculate_wall_flux(flex_cell_center *fcell, double u); 

int update_for_flex_cell(flex_cell_center *fcell, double u, double dt);

int inviscid_source_vector_for_flex_cell(flex_cell_center *fcell); 

int update_flex_cell_geometry_based_on_vertex_movement(flex_cell_center *src);

double signal_frequency_for_flex_cell(flex_cell_center *fcell, int dimensions=2);

int print_data_for_flex_cell(flex_cell_center *fcell); 

#endif
