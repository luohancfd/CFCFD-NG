/** \file parser.h
 *  \ingroup util
 *  \brief Header file for the parsing functions.
 *
 *  \author Rowan J Gollan (r.gollan@uq.edu.au)
 *  \date 15-Mar-2005
 **/

#include "../iniparser/src/iniparser.h"

int missing_item(char *item, char *section, char *fname, char *msg);
int parse_double(dictionary *d, char *section, char *item, char *fname, double *value);
int parse_int(dictionary *d, char *section, char *item, char *fname, int *value);
int parse_boolean(dictionary *d, char *section, char *item, char *fname, int *value);
int parse_string(dictionary *d, char *section, char *item, char *fname, char *value);

int parse_array_of_int_from_string(const char *in_string, 
				   int arry[], int n, int def_value);
int parse_array_of_double_from_string(const char *in_string, 
				      double arry[], int n, double def_value);
