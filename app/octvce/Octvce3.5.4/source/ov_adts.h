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

/**\file Prototypes for operations on ADTs specific to cartesian cells*/

void traverse_volume_ADT(double[], double, int *, Point_tnode, List_bbox *, Body, Body *, short int *, double []);

void traverse_area_ADT(short int, double [], double [], double, double, int [], Point_tnode, List_bbox *, Body, double[], int[]);

short int build_volume_ADT(void);

short int build_area_ADT(void);

short int add_to_point_ADT(double[], short int, Point_tnode *, int); 


