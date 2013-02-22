/*
 * memread.c - write an address to it, then read back the byte at that address
 *
 *        Requires FUSD to run.
 *
 *        Zane Smith, 2004
 */




#ifndef __KERNEL__
#define __KERNEL__
#endif
/*#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <linux/kernel.h>
#include <asm/io.h>*/

#include <linux/module.h>
#include <linux/delay.h>
#include <linux/errno.h>
#include <linux/kernel.h>
#include <linux/malloc.h>
#include <linux/poll.h>
//#include <linux/segment.h>
#include <asm/io.h>
#include <linux/interrupt.h>
#include <linux/i2c.h>
#include <stdio.h>

#include "fusd.h"

#define MIN(x,y) ((x)<(y) ? (x) : (y))

int data = 255;
int adrs = 0;
int data_length = 0;

int memread_read(struct fusd_file_info *file, char *user_buffer,
	      size_t user_length, loff_t *offset)
{
  int retval = 0;
  int temp = 0;
  char buf[16];
  int fd;

  if (*offset > 0) return 0;
  
  /*temp = isa_readb(adrs);*/
  isa_memcpy_fromio(temp,adrs,1);
  sprintf(buf,"0x%x\n",temp);

  /* fd = open("/dev/mem",O_RDONLY,0);
  if ( lseek ( fd, adrs, 0) != adrs )
  {
    perror ("lseek: ");
    close (fd);
    exit (-1);
  }
  if ( read (fd, &temp, 1) != 1)
  {
    perror ("read: ");
    close (fd);
    exit (-2);
  }
  sprintf(buf,"0x%x\n",temp);
  close(fd);*/

  user_length = MIN(user_length,strlen(buf));
  memcpy(user_buffer, buf, user_length);
    *offset += user_length;
  return user_length;
}

ssize_t memread_write(struct fusd_file_info *file, const char *user_buffer,
		   size_t user_length, loff_t *offset)
{
  /* set the pointer data to the value of user_buffer */
  adrs = (int)atoi(user_buffer);
  data_length = user_length;
  return user_length;
}

int do_open_or_close(struct fusd_file_info *file)
{
  /*printf("Opened/closed\n");*/
  return 0; /* opens and closes always succeed */
}


struct fusd_file_operations echo_fops = {
  open: do_open_or_close,
  read: memread_read,
  write: memread_write,
  close: do_open_or_close
};


int main(int argc, char *argv[])
{
  char stringtest[100] = "Hello";
  sprintf(stringtest,"%s there!",stringtest);
  fprintf(stderr,stringtest);
  if (fusd_register("/dev/memread", 0666, NULL, &echo_fops) < 0) {
    perror("register of /dev/tranrec failed");
    exit(1);
  }

  fprintf(stderr, "calling fusd_run...\n");
  fusd_run();
  return 0;
}
