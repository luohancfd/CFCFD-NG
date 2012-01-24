/*Copyright (C) Joseph Tang, 2007
                                                                                                        
    This file is part of OctVCE.
                                                                                                        
    OctVCE is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
                                                                                                        
    OctVCE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
                                                                                                        
    You should have received a copy of the GNU General Public License
    along with OctVCE.  If not, see <http://www.gnu.org/licenses/>.
*/

/**\file Function prototypes for reading/writing data to files*/

short int read_ov_param(char []);

short int read_ov_gas(char []);

short int read_IC(char []);

short int read_BC(char []);

short int get_unsteady_IC(char [], IC_heat *);

short int get_unsteady_BC(char [], BC_data *);

short int write_soln(int, int, double);

short int link_solution_files(char [], double);

short int write_oct_flow_soln(List_leaf, List_leaf, short int, int, double, char []);

short int write_history_loc(Cart_cell, double [], double, char []);

short int output_history_loc(double);

void close_history_files(void);

int get_num_boundary_cells(List_leaf, List_leaf); 

short int read_deton(char []);

void get_max_schl(void);

void write_mesh_and_soln(Cart_cell, int, double);

short int read_mesh(char *);
