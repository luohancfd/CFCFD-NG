/** \file logfile.h
 * \brief Message logging routines.
 * \ingroup util
 */

int open_log_file  (char *log_file_name);
int log_message    (char *message, int screen);
int close_log_file (void);


