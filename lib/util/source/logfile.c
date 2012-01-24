/** \file logfile.c
 * \brief Message logging functions.
 * \ingroup util
 *
 * \author PA Jacobs
 */

#include<stdio.h>
#include <stdlib.h>
#include "compiler.h"

/* Static data for the message logging functions. */
static FILE *log_file;
static int  log_file_is_open = 0;

/*----------------------------------------------------------*/

/** Open a log file of the specified name.
 * \param *log_file_name  : name of the file, relative to the current directory
 */
int open_log_file (char *log_file_name)
{
    log_file = fopen (log_file_name, "w");
    if (log_file == NULL) {
	printf ("Could not open log_file: %s\n", log_file_name);
	exit (-1);
    }
    log_file_is_open = 1;
    return log_file_is_open;
}

/** Logs a message and optionally puts it on the screen.
 * \param *message  : pointer to the message string
 * \param screen    : =0  don't put the message on the screen,
 *                    =1  put the message onto the screen
 */
int log_message (char *message, int screen)
{
    if (log_file_is_open) {
	fprintf (log_file, "%s", message);
	fflush (log_file);
    }
    if (screen == 1) {
	printf ("%s", message);
	fflush (stdout);
    }
    return 0;
}

/** Cleans up.
 */
int close_log_file (void)
{
    if (log_file_is_open && log_file != NULL) fclose (log_file);
    return 0;
}

/*------------------ end of logfile.c --------------------------*/

