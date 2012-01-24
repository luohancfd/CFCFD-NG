/***************************************************************************
#*   Copyright (C) 2009-2011 by Steve Walter, Oliver Borm, Franz Blaim     *
#*                              Lionel Gamet                               *
#*   steve.walter@mytum.de, oli.borm@web.de, franz.blaim@gmx.de,           *
#*   lionel.gamet@fluorem.com                                              *
#*                                                                         *
#*   This program is free software; you can redistribute it and/or modify  *
#*   it under the terms of the GNU General Public License as published by  *
#*   the Free Software Foundation; either version 3 of the License, or     *
#*   (at your option) any later version.                                   *
#*                                                                         *
#*   This program is distributed in the hope that it will be useful,       *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#*   GNU General Public License for more details.                          *
#*                                                                         *
#*   You should have received a copy of the GNU General Public License     *
#*   along with this program; if not, write to the                         *
#*   Free Software Foundation, Inc.,                                       *
#*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
#***************************************************************************

# Author: Steve Walter, Franz Blaim, Oliver Borm, Lionel Gamet
# Date: December 2010
*/

%module CGNS
%{
#include "cgnslib.h"
%}

/* Variadic argument treatments */
%varargs(const char *argcha = NULL, int arginta=0,
         const char *argchb = NULL, int argintb=0,
         const char *argchc = NULL, int argintc=0,
         const char *argchd = NULL, int argintd=0) cg_goto;
%varargs(const char *argcha = NULL, int arginta=0,
         const char *argchb = NULL, int argintb=0,
         const char *argchc = NULL, int argintc=0,
         const char *argchd = NULL, int argintd=0) cg_gorel;

/* SWIG include files section */
%include "carrays.i"
%include "cstring.i"
%include "cpointer.i"

/* Special treatment of functions that internally allocate char** arguments */
%cstring_output_allocate(char **connectname     , free(*$1));
%cstring_output_allocate(char **descr_text      , free(*$1));
%cstring_output_allocate(char **donorname       , free(*$1));
%cstring_output_allocate(char **filename        , free(*$1));
%cstring_output_allocate(char **geo_file        , free(*$1));
%cstring_output_allocate(char **link_path       , free(*$1));
%cstring_output_allocate(char **NormDefinitions , free(*$1));
%cstring_output_allocate(char **StateDescription, free(*$1));
%cstring_output_allocate(char **zonename        , free(*$1));

/* Special treatment of functions that return strings in char* arguments */
%cstring_bounded_output(char *ArrayName       , 1024);
%cstring_bounded_output(char *basename        , 1024);
%cstring_bounded_output(char *bitername       , 1024);
%cstring_bounded_output(char *boconame        , 1024);
%cstring_bounded_output(char *CAD_name        , 1024);
%cstring_bounded_output(char *connectname     , 1024);
%cstring_bounded_output(char *coordname       , 1024);
%cstring_bounded_output(char *descr_name      , 1024);
%cstring_bounded_output(char *discrete_name   , 1024);
%cstring_bounded_output(char *donorname       , 1024);
%cstring_bounded_output(char *fambc_name      , 1024);
%cstring_bounded_output(char *family_name     , 1024);
%cstring_bounded_output(char *fieldname       , 1024);
%cstring_bounded_output(char *geo_name        , 1024);
%cstring_bounded_output(char *gridname        , 1024);
%cstring_bounded_output(char *holename        , 1024);
%cstring_bounded_output(char *IntegralDataName, 1024);
%cstring_bounded_output(char *name            , 1024);
%cstring_bounded_output(char *node_name       , 1024);
%cstring_bounded_output(char *part_name       , 1024);
%cstring_bounded_output(char *RegionName      , 1024);
%cstring_bounded_output(char *SectionName     , 1024);
%cstring_bounded_output(char *solname         , 1024);
%cstring_bounded_output(char *user_data_name  , 1024);
%cstring_bounded_output(char *zitername       , 1024);
%cstring_bounded_output(char *zonename        , 1024);

/* Skip functions that we know to not work */
/* - Functions not working because of character strings arrays */
%ignore cg_golist;
%ignore cg_where;

/* Include file section */
%include "cgnslib.h"

/* Pointer wrapper */
%array_class(int,intArray)
%array_class(float,floatArray)
%array_class(double,doubleArray)
%pointer_class(int,intp);
%pointer_class(float,floatp);
%pointer_class(double,doublep);

%pointer_class(               ZoneType_t,                ZoneType_tp);
%pointer_class(             AngleUnits_t,              AngleUnits_tp);
%pointer_class(              MassUnits_t,               MassUnits_tp);
%pointer_class(            LengthUnits_t,             LengthUnits_tp);
%pointer_class(              TimeUnits_t,               TimeUnits_tp);
%pointer_class(       TemperatureUnits_t,        TemperatureUnits_tp);
%pointer_class(   ElectricCurrentUnits_t,    ElectricCurrentUnits_tp);
%pointer_class( LuminousIntensityUnits_t,  LuminousIntensityUnits_tp);
%pointer_class(              DataClass_t,               DataClass_tp);
%pointer_class(           GridLocation_t,            GridLocation_tp);
%pointer_class(   GridConnectivityType_t,    GridConnectivityType_tp);
%pointer_class(             BCDataType_t,              BCDataType_tp);
%pointer_class(           PointSetType_t,            PointSetType_tp);
%pointer_class( GoverningEquationsType_t,  GoverningEquationsType_tp);
%pointer_class(              ModelType_t,               ModelType_tp);
%pointer_class(                 BCType_t,                  BCType_tp);
%pointer_class(               DataType_t,                DataType_tp);
%pointer_class(            ElementType_t,             ElementType_tp);
%pointer_class(ArbitraryGridMotionType_t, ArbitraryGridMotionType_tp);
%pointer_class(         SimulationType_t,          SimulationType_tp);
%pointer_class(       WallFunctionType_t,        WallFunctionType_tp);
%pointer_class(               AreaType_t,                AreaType_tp);
%pointer_class(   AverageInterfaceType_t,    AverageInterfaceType_tp);
