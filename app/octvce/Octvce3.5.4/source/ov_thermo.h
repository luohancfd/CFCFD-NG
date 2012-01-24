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

/**\file Prototypes for thermodynamic evaluations e.g. sound speed, mixture quantities etc*/

double get_mixture_soundspd(State_vector *);

double get_mixture_RhoE(State_vector *, int);

void get_mixture_Pres(State_vector *);

double get_mixture_T(double, double, double, double, double, double, double);

double JWLB_q(double, double [], double []);

double JWLB_dqdr(double, double, double [], double []);

void precompute_JRL(void);

short int treat_products_as_air(double, double);

void get_smallest_JWLB_terms(void);
