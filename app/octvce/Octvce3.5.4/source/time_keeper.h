#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXSTRING 100

typedef struct {
  clock_t begin_clock, save_clock;
  time_t begin_time, save_time;
} time_keeper;

void start_time(void);

double prn_time(void);
