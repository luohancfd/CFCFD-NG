/** \file parser.c
 *  \brief Parsing functions to add error checking to iniparser.c
 *  \ingroup util
 *
 * The iniparser library provides some convenient methods for 
 * parsing an .ini-type file BUT it lacks some error checking
 * which will be useful for our work.
 *
 * This module provides some extra functions that provide
 * the error checking and log the errors. (But is assumed
 * that the calling routine has already opened a logfile).
 *
 * \author Rowan J Gollan
 * \date 04-Mar-2005
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../iniparser/src/iniparser.h"
#include "useful.h"
#include "parser.h"

#define NOTFOUND -1
#define EMPTY "empty"
#define KDIM 64
#define MDIM 256

/**
 * \brief Create a message about the missing item in an .ini file.
 *
 * This function places a message string in *msg for use by the logging
 * function.
 *
 * \author Rowan J Gollan
 * \date 04-Mar-2005
 **/

int missing_item(char *item, char *section, char *fname, char *msg)
{
    char tmp_a[MDIM];

    sprintf(msg, "There is a problem reading \'%s\' from the [%s] section.\n", 
	    item, section);
    sprintf(tmp_a, "Check to see if it is listed in %s.\n", fname);
    msg = strcat(msg, tmp_a);

    return (0);
}

/**
 * \brief Parse a double and check for errors.
 *
 * This function parses a double from a .ini file using the
 * iniparser library.
 * The caller should test for success or failure.
 * The user passes the pointer to the variable that they want filled
 * in with the double.
 *
 * \return int SUCCESS or FAILURE
 *
 * \author Rowan J Gollan
 * \date 04-Mar-2005
 **/

int parse_double(dictionary *d, char *section, char *item, char *fname, double *value)
{
    char key[KDIM];
    char msg[MDIM];

    sprintf(key, "%s:%s", section, item);
    *value = iniparser_getdouble(d, key, NOTFOUND);
    if (*value == NOTFOUND) {
	missing_item(item, section, fname, msg);
	printf("%s", msg);
	return (FAILURE);
    }

    return (SUCCESS);
}

/**
 * \brief Parse an integer and check for errors.
 *
 * This function parses an integer from a .ini file using the
 * iniparser library.
 * The caller should test for success or failure.
 * The user passes the pointer to the variable that they want filled
 * in with the integer.
 *
 * \return int SUCCESS or FAILURE
 *
 * \note Rowan, the integer values for SUCCESS and FAILURE may clash
 *       with intended values.  Maybe the user of these functions needs
 *       to nominate what values should be used for this flag.
 *
 * \author Rowan J Gollan
 * \date 04-Mar-2005
 **/

int parse_int(dictionary *d, char *section, char *item, char *fname, int *value)
{
    char key[KDIM];
    char msg[MDIM];

    sprintf(key, "%s:%s", section, item);
    *value = iniparser_getint(d, key, NOTFOUND);
    if (*value == NOTFOUND) { 
 	missing_item(item, section, fname, msg);
	printf("%s", msg);
 	return (FAILURE); 
    } 
    return (SUCCESS);
}

/**
 * \brief Parse a boolean and check for errors.
 *
 * This function parses a boolean from a .ini file using the
 * iniparser library.
 * The caller should test for success or failure.
 * The user passes the pointer to the variable that they want filled
 * in with the integer (1 or 0 representing the boolean True or False).
 *
 * \return int SUCCESS or FAILURE
 *
 * \author Rowan J Gollan
 * \date 04-Mar-2005
 **/

int parse_boolean(dictionary *d, char *section, char *item, char *fname, int *value)
{
    char key[KDIM];
    char msg[MDIM];

    sprintf(key, "%s:%s", section, item);
    *value = iniparser_getboolean(d, key, NOTFOUND);
    if (*value == NOTFOUND) {
	missing_item(item, section, fname, msg);
	printf("%s", msg);
	return (FAILURE);
    }

    return (SUCCESS);
}

/**
 * \brief Parse a string and check for errors.
 *
 * This function parses a string from a .ini file using the
 * iniparser library.
 * The caller should test for success or failure.
 * The user passes the pointer to the variable that they want filled
 * in with the string.
 * The user should be careful that his/her pointer can contain
 * the length of string expected.
 *
 * \return int SUCCESS or FAILURE
 *
 * \author Rowan J Gollan
 * \date 04-Mar-2005
 **/

int parse_string(dictionary *d, char *section, char *item, char *fname, char *value)
{
    char key[KDIM];
    char msg[MDIM];
   
    sprintf(key, "%s:%s", section, item);
    strcpy(value, iniparser_getstring(d, key, EMPTY));
    if (strncmp(value, EMPTY, 5) == 0) {
	missing_item(item, section, fname, msg);
	printf("%s", msg);
	return (FAILURE);
    }

    return (SUCCESS);
}

/*---------------------------------------------------------------------------*/

/* Some extras for elmer... */

/** \brief Attempt to extract space-separated integer values from a string.
 *
 * \param strng     : pointer to the string we want to pull apart
 * \param arry      : array for storage of the extracted values
 * \param n         : try to get n values
 * \param def_value : value to be used if extraction from the string fails
 *
 * \author PJ
 * \version 08-Jun-2005
 */
int parse_array_of_int_from_string(const char *in_string, 
				   int arry[], int n, int def_value)
{
    int i;
    int n_success = 0;
    char local_string[256], *token;
    strncpy(local_string, in_string, sizeof(local_string)-1);
    local_string[sizeof(local_string)-1] = '\0';
    for ( i = 0; i < n; ++i ) {
	arry[i] = def_value;
    }
    token = strtok(local_string, "\t ");
    for ( i = 0; i < n; ++i ) {
	if ( token == NULL ) break;
	if ( sscanf(token, "%d", &(arry[i])) < 0 ) {
	    /* Failed conversions get the default value. */
	    arry[i] = def_value;
	} else {
	    ++n_success;
	}
	token = strtok(NULL, "\t "); /* Point to following token. */
    }
    return n_success;
}


/** \brief Attempt to extract space-separated integer values from a string.
 *
 * \param strng     : pointer to the string we want to pull apart
 * \param arry      : array for storage of the extracted values
 * \param n         : try to get n values
 * \param def_value : value to be used if extraction from the string fails
 *
 * \author PJ
 * \version 08-Jun-2005
 */
int parse_array_of_double_from_string(const char *in_string, 
				      double arry[], int n, double def_value)
{
    int i;
    int n_success = 0;
    char local_string[256], *token;
    strncpy(local_string, in_string, sizeof(local_string)-1);
    local_string[sizeof(local_string)-1] = '\0';
    for ( i = 0; i < n; ++i ) {
	arry[i] = def_value;
    }
    token = strtok(local_string, "\t ");
    for ( i = 0; i < n; ++i ) {
	if ( token == NULL ) break;
	if ( sscanf(token, "%lf", &(arry[i])) < 0 ) {
	    /* Failed conversions get the default value. */
	    arry[i] = def_value;
	} else {
	    ++n_success;
	}
	token = strtok(NULL, "\t "); /* Point to following token. */
    }
    return n_success;
}

