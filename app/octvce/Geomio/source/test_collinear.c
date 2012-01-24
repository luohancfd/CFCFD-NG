#include <stdio.h>
#include "gio_kernel.h"
#include "gio_adts.h"

main()
{
  double p0[3], p1[3], p2[3];

  p0[0] = 0;
  p0[1] = 1;
  p0[2] = 0;

  p1[0] = 0;
  p1[1] = 2;
  p1[2] = 0;

  p2[0] = 0;
  p2[1] = 3;
  p2[2] = -0.00001;

  if(collinear_test(p0,p1,p2) == TRUE)
    printf("Points are collinear\n");
  else printf("Points are NOT collinear\n");
}
