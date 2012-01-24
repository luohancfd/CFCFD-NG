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

/**\file Prototypes for establishing/checking connectivities between cells - fluxes, neighbours and verticies*/

void establish_neighbour_after_refine(Cart_cell[], short int);

void establish_neighbour_after_coarsen(Cart_cell, short int);

short int check_neighbour_for_refine(Cart_cell, Cart_cell[]);

short int check_cell_for_coarsening(Cart_cell);

void assign_new_face_neighbour(Cart_cell, short int, Cart_cell);

short int assign_vtx_after_refine(Cart_cell [], short int, List_vtx_glob *, List_vtx_glob *);

short int assign_vtx_after_coarsen(Cart_cell [], short int, List_vtx_glob *, List_vtx_glob *);

void count_verticies(void); 

Cart_cell find_vtx_cell(Cart_cell [], double [], double [][3]);

short int alloc_fluxes(List_leaf, List_leaf, short int);

void reassign_fluxes_after_refine(Cart_cell [], short int);

void reassign_fluxes_after_coarsen(Cart_cell [], short int);
