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

/**\file Prototypes for routines associated with adaptation (except for establishing cell and flux connectivities) */

void compute_adapt_param(List_leaf, List_leaf);

short int adapt(int *, int *);

void flag_for_refine(List_leaf, List_leaf);

void flag_for_coarsen(List_leaf, List_leaf);

int par_refine(List_leaf *, List_leaf *, List_vtx_glob *, List_vtx_glob *, short int);

int par_coarsen(List_leaf *, List_leaf *, List_vtx_glob *, List_vtx_glob *, short int);

short int refine(Cart_cell, short int, short int, short int, List_leaf *);

short int coarsen(Cart_cell, Cart_cell);

Cart_cell make_root(void);

short int make_children(Cart_cell);

short int refine_specific(Cart_cell, short int, char *, Cart_cell *);
