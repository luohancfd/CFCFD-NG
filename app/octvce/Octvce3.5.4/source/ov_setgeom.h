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

/**\file Prototypes for setting geometric properties of octree cartesian VCE cells*/

void compute_children_centroids(Cart_cell, double[][3]);

void compute_cell_centroid(Cart_cell, double[][3]);

void get_vtx_loc(Cart_cell, short int, double []);

short int set_cell_geom_data(Cart_cell[], short int, short int, short int, short int, short int, short int); 

short int set_parent_geom_data(Cart_cell[]);

short int set_cell_type(Cart_cell, double, double, List_bbox *, Body); 

short int set_cell_additional_refine_times(Cart_cell, List_bbox *, Body); 

short int set_face_flux_area(Cart_cell, short int, double, double, List_bbox *, Body, short int); 

short int set_IC_flag(Cart_cell[]);

void push_neighbour_flux_areas(Cart_cell, short int, double);

short int check_if_in_deton_spheres(double, double []);

void comp_new_deton_cells(double, List_leaf, List_leaf);
