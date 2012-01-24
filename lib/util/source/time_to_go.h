/** \file time_to_go.h
 * \ingroup util
 * \brief Function prototype for time_to_go.c.
 */

#ifndef TIME_TO_GO_ALREADY_INCLUDED

#include <time.h>

char *time_to_go(time_t start_wall_clock,
                 time_t current_wall_clock,
                 int current_step,
                 int max_steps,
                 double current_dt, double current_t, double final_t);

#define TIME_TO_GO_ALREADY_INCLUDED
#endif
