#include <stdio.h>
#include <string.h>
#include "gio_kernel.h"
#include "gio_adts.h"
#include "gio_io.h"

main()
{
  int num_files = 0;
  char filename[24];
  strcpy(filename, "bodies.pvtp");
  num_files = read_body_files(filename);
  printf("Number of files read - %d\n",num_files);

  return(0);
} 
