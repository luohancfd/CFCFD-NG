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

/**\file Prototypes for functions related to reconstruction, including gradient computation*/

void get_grads_all_cells(List_leaf, List_leaf, short int);

double compute_spatial_grads(Cart_cell, double[][5]);

void compute_limiter(List_leaf, List_leaf, short int);

void reconstruct_refined_children(Cart_cell []);

void reconstruct_coarsened_cell(Cart_cell);

State_vector reconstruct_flow(Cart_cell, short int, short int, short int, short int, int, double *, double *, char);
