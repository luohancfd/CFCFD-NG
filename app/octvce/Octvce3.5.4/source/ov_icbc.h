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

/**\file Prototypes for setting flow ICs and BCs for octree cartesian cells*/

State_vector set_boundary_state(State_vector *, State_vector *, double [], short int, short int, short int, double, double, Grad *, double);

void calc_trans_IC(double, IC_heat *, short int);

void calc_trans_BC(double, BC_data *, short int);

void set_cell_IC(Cart_cell);


