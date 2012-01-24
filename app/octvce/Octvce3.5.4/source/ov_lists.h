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

/**\file Prototypes for list functions*/

short int add_to_leaf_list_front(Cart_cell, List_leaf *, List_leaf *, List_leaf *);

short int add_to_merge_list_front(Cart_cell, List_merge *);

void delete_from_leaf_list(Cart_cell, List_leaf *, List_leaf *); 

void delete_from_merge_list(Cart_cell, List_merge *);

short int add_to_Vtxlist(Vertex, List_vtx_glob *, List_vtx_glob *);

void delete_from_Vtxlist(Vertex, List_vtx_glob *, List_vtx_glob *);

void modify_lists(List_leaf[], List_leaf[]);

void break_list(List_leaf[], List_leaf[]);

void break_vtx_list(List_vtx_glob[], List_vtx_glob[]);

void merge_lists(List_leaf[], List_leaf[]);

void merge_vtx_lists(List_vtx_glob[], List_vtx_glob[]);

void delete_verticies(List_vtx_glob *, List_vtx_glob *);

void count_deletable_verticies(List_vtx_glob *, List_vtx_glob *, int, int);

void check_duplicate_cells(List_leaf, List_leaf, int, int);

void check_flux_connect(List_leaf, List_leaf, int);

void printout_adapt_cells(short int, short int, List_leaf, List_leaf, int, int);

void reset_adapt_flag(List_leaf, short int);

void printout_cells_2b_refined(List_leaf, int);

void printout_cells_2b_coarsened(List_leaf, int);
