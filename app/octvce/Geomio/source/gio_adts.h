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

/**\file Prototypes for constructing and using object ADTs for point inclusion interrogation*/

int build_body_ADT(void);

short int build_face_ADT(Body Solid, short int);

short int build_edge_ADT(Body2D Poly);

short int traverse_body_ADT(double [], double [], Body_tnode, List_bbox *); 

short int traverse_edge_ADT(double[], double [], double [], Edge1D_tnode, int *, double *);

short int traverse_face_ADT(double[], double [], double [], Body2D_tnode, double *);

short int add_to_body_ADT(Body, double[], double[][4], Body_tnode *, int);

short int add_to_edge_ADT(Edge3D, double[], double[][1], Edge1D_tnode *, int);

short int add_to_face_ADT(Body2D, double [], double[][2], Body2D_tnode *, int);
