/** \file time_to_go.c
 * \ingroup util
 * \brief Function to estimate how much longer we need to wait for the
 *        simulation to finish.
 */

#include <time.h>
#include <stdio.h>
#include "time_to_go.h"

/** \brief Return a message indicating the estimated wall clock to completion.
 *
 * There will be two estimates:
 * WCtFT == time to final simulation time;
 * WCtMS == time to maximum number of steps
 */
char *time_to_go(time_t start_wall_clock,
                 time_t current_wall_clock,
                 int current_step,
                 int max_steps,
                 double current_dt, double current_t, double final_t)
{
    double WC_elapsed, WC_per_step, WCtFT, WCtMS;
    static char time_to_go_string[132];

    WC_elapsed = (double) (current_wall_clock - start_wall_clock);
    WC_per_step = WC_elapsed / current_step;
    WCtFT = (final_t - current_t) / current_dt * WC_per_step;
    WCtMS = (max_steps - current_step) * WC_per_step;

    sprintf(time_to_go_string, "WC=%.1f WCtFT=%.1f WCtMS=%.1f",
            WC_elapsed, WCtFT, WCtMS);
    return time_to_go_string;
}   /* end time_to_go() */

/*-----------------------------------------------------------------*/
