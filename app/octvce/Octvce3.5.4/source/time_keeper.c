#include <stdio.h>
#include <stdlib.h>
#include "time_keeper.h"

static time_keeper tk;

void start_time(void) {
  tk.begin_clock = tk.save_clock = clock();
  tk.begin_time = tk.save_time = time(NULL);
}

double prn_time(void) {
  /*nt field_width, n1, n2;*/
  double clocks_per_sec = (double) CLOCKS_PER_SEC, user_time, real_time;

  user_time = (clock() - tk.save_clock)/clocks_per_sec;
  real_time = difftime(time(NULL), tk.save_time);
  tk.save_clock = clock();
  tk.save_time = time(NULL);

  /*1 = sprintf(s1, "%lf", user_time);
  n2 = sprintf(s2, "%lf", real_time);
  field_width = (n1 > n2) ? n1:n2;*/
  
  printf("User time: %g\n", user_time);
  printf("Real time: %g\n", real_time);

  return user_time;
}
