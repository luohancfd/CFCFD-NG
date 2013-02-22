/*
 * tr.c - Transient Recorder user-space device driver
 *
 *        Requires FUSD to run.
 *
 *        Zane Smith, 2004
 */




#ifndef __KERNEL__
#define __KERNEL__
#endif
#ifndef MODULE
#define MODULE
#endif

#include <asm/io.h>
#include <asm/system.h>
#include <asm/bitops.h>
#include <linux/module.h>
#include <linux/errno.h>
#include <linux/kernel.h>
#include <linux/version.h>
#include <linux/ioport.h>
#include <linux/slab.h>
#include <linux/poll.h>
#include <linux/mman.h>
#include <stdio.h>

#include "fusd.h"

#define MIN(x,y) ((x)<(y) ? (x) : (y))

#define ERROR_MODE 0
#define INPUT_MODE 1
#define ARM_MODE 2
#define TRIGGER_MODE 3
#define READ_TR_STATUS_MODE 4
#define READ_TB_STATUS_MODE 5
#define READ_AD_S_MODE 6
#define READ_AD_BUFFER_MODE 7
#define DIRECT_READ_MODE 8
#define DIRECT_WRITE_MODE 9
#define DIRECT_WRITE_WRITE 10
#define CARD_SELECT_MODE 11
#define GET_CARD_MODE 12

#define TR_BUFFER_WINSIZE 16384
#define TR_BUFFER_START 0xd8000
#define PORT_START 0x000
#define PORT_SIZE 0x32f

#define DEBSTR if(debug)

int data = 255;
int adrs = 0;
int mode = INPUT_MODE;
int data_length = 0;
int channelcard = 0x14;
int tempthing=0;
int debug=1;


/* Writes writebyte to the port at writeaddress */
void dooutput(int writeaddress, int writebyte)
{
  int fd;

  fd = open("/dev/port", O_RDWR);
  if (fd == -1)
    {
    perror("open /dev/port");
    exit(-1);
    }

  if ( lseek ( fd, writeaddress, 0) != writeaddress )
    {
      perror ("lseek: ");
      close (fd);
      exit (-1);
    }
  if ( write (fd, &writebyte, 1) != 1)
    {
      perror ("write: ");
      close (fd);
      exit (-2);
    }
  DEBSTR fprintf(stderr,"Set port 0x%x\n",writeaddress);
  close (fd);
}

/* Copies the byte in the port at readaddress into container */
void doinput(int readaddress, int *container)
{
  int fd;

  fd = open("/dev/port", O_RDONLY);
  if (fd == -1)
    {
    perror("open /dev/port");
    exit(-1);
    }

  if ( lseek ( fd, readaddress, 0) != readaddress )
    {
      perror ("lseek: ");
      close (fd);
      exit (-1);
    }
  if ( read (fd, container, 1) != 1)
    {
      perror ("read: ");
      close (fd);
      exit (-2);
    }
  DEBSTR fprintf(stderr,"Got port 0x%x\n",readaddress);
  close (fd);
}

/* Method to be called on a device read" */
int tranrec_read(struct fusd_file_info *file, char *user_buffer,
	      size_t user_length, loff_t *offset)
{
  int retval = 0;
  int temp = 0;
  char buf[128];
  int fd;
  unsigned char *isa_mem;

  if (mode==READ_AD_BUFFER_MODE)
    {
      if (*offset >= TR_BUFFER_WINSIZE) return 0;
    dooutput(0x31e,0x80);
    dooutput(0x320+2*(channelcard%16),channelcard);
    fd = open("/dev/mem", O_RDONLY);
    if (fd == -1)
      {
      perror("open /dev/mem");
      exit(-1);
      }
    isa_mem = mmap(0, TR_BUFFER_WINSIZE, PROT_READ, MAP_SHARED, fd, TR_BUFFER_START);
    if ((long)isa_mem == -1)
      {
      perror("memmap /dev/mem");
      exit(-1);
      }
    DEBSTR fprintf(stderr,"Memory mapped\n");
    user_length = MIN(user_length,TR_BUFFER_WINSIZE-*offset);
    memcpy(user_buffer, isa_mem+*offset, user_length);
    close(fd);
    }
  else
    {
    if (*offset > 0) return 0;

    switch(mode)
      {
      case READ_TR_STATUS_MODE:
        doinput(0x310,&temp);
        (temp<17) ? sprintf(buf,"0%x",temp) : sprintf(buf,"%x",temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        dooutput(0x310,15);
        doinput(0x311,&temp);
        (temp<17) ? sprintf(buf,"%s0%x",buf,temp) : sprintf(buf,"%s%x",buf,temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        dooutput(0x310,16);
        doinput(0x311,&temp);
        (temp<17) ? sprintf(buf,"%s0%x",buf,temp) : sprintf(buf,"%s%x",buf,temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        dooutput(0x310,17);
        doinput(0x311,&temp);
        (temp<17) ? sprintf(buf,"%s0%x",buf,temp) : sprintf(buf,"%s%x",buf,temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        mode=INPUT_MODE;
        break;
      case READ_TB_STATUS_MODE:
        doinput(0x312,&temp);
        (temp<17) ? sprintf(buf,"0%x",temp) : sprintf(buf,"%x",temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        doinput(0x314,&temp);
        (temp<17) ? sprintf(buf,"%s0%x",buf,temp) : sprintf(buf,"%s%x",buf,temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        doinput(0x316,&temp);
        (temp<17) ? sprintf(buf,"%s0%x",buf,temp) : sprintf(buf,"%s%x",buf,temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        doinput(0x313,&temp);
        (temp<17) ? sprintf(buf,"%s0%x",buf,temp) : sprintf(buf,"%s%x",buf,temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        doinput(0x315,&temp);
        (temp<17) ? sprintf(buf,"%s0%x",buf,temp) : sprintf(buf,"%s%x",buf,temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        doinput(0x317,&temp);
        (temp<17) ? sprintf(buf,"%s0%x",buf,temp) : sprintf(buf,"%s%x",buf,temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        mode=INPUT_MODE;
        break;
      case READ_AD_S_MODE:
	dooutput(0x31e,0x80);
        doinput(0x320+2*(channelcard%16),&temp);
        (temp<17) ? sprintf(buf,"0%x",temp) : sprintf(buf,"%x",temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        doinput(0x321+2*(channelcard%16),&temp);
        (temp<17) ? sprintf(buf,"%s0%x",buf,temp) : sprintf(buf,"%s%x",buf,temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        mode=INPUT_MODE;
        break;
      case DIRECT_READ_MODE:
        doinput(adrs,&temp);
        sprintf(buf,"%c\n",temp);
	DEBSTR fprintf(stderr,"0x%x\n",temp);
        mode=INPUT_MODE;
        break;
      default:
        fprintf(stderr,"Read error: No read pending\n");
        exit(-1);
      }
    user_length = MIN(user_length,strlen(buf));
    memcpy(user_buffer, buf, user_length);
    }
  

  *offset += user_length;
  return user_length;
}

/* Method to be called on device write:
 * If we are in input mode, change to the new mode.
 * Othewise, perform the required action*/
ssize_t tranrec_write(struct fusd_file_info *file, const char *user_buffer,
		   size_t user_length, loff_t *offset)
{
  int temp = 0;

  temp = (int)atoi(user_buffer);
  
  switch(mode)
  {
    case INPUT_MODE:
      mode=temp;
      switch(mode)
      {
        case CARD_SELECT_MODE:
	  DEBSTR fprintf(stderr,"Write channel&card number\n");
	  mode=GET_CARD_MODE;
	  break;
        case ARM_MODE:
	  dooutput(0x310,0x20);  //784,32
	  DEBSTR fprintf(stderr,"Trigger Unit armed\n");
	  mode=INPUT_MODE;
	  break;
        case TRIGGER_MODE:
	  dooutput(0x310,0x10);  //784,16
	  DEBSTR fprintf(stderr,"Trigger Unit triggered\n");
	  mode=INPUT_MODE;
	  break;
        case READ_TR_STATUS_MODE:
	  DEBSTR fprintf(stderr,"Read device to get Trigger Unit status registers\n");
	  break;
        case READ_TB_STATUS_MODE:
	  DEBSTR fprintf(stderr,"Read device to get Timebase Unit status registers\n");
	  break;
        case READ_AD_S_MODE:
       	  DEBSTR fprintf(stderr,"Read device to get A/D Unit S registers\n");
      	  break;
        case READ_AD_BUFFER_MODE:
	  DEBSTR fprintf(stderr,"Read device to get A/D Unit buffers\n");
	  tempthing=0;
       	  break;
        case DIRECT_READ_MODE:
	  DEBSTR fprintf(stderr,"DIRECT MODE: Write port to read from\n");
	  break;
        case DIRECT_WRITE_MODE:
	  DEBSTR fprintf(stderr,"DIRECT MODE: Write port to write to\n");
	  break;
        default:
	  DEBSTR fprintf(stderr,"Error: Unsupported action\n");
	  exit(-1);
      }
      break;
    case GET_CARD_MODE:
      channelcard=temp;
      DEBSTR fprintf(stderr,"Channel %d, card %d selected\n",channelcard/16,channelcard%16);
      mode=INPUT_MODE;
      break;
    case DIRECT_READ_MODE:
      adrs=temp;
      DEBSTR fprintf(stderr,"DIRECT MODE: Read device to get byte from port\n");
      break;
    case DIRECT_WRITE_MODE:
      adrs=temp;
      mode=DIRECT_WRITE_WRITE;
      DEBSTR fprintf(stderr,"DIRECT MODE: Write data to be written to port\n");
      break;
    case DIRECT_WRITE_WRITE:
      dooutput(adrs,temp);
      mode=INPUT_MODE;
      break;
    default:
      DEBSTR fprintf(stderr,"ERROR: Incorrect mode\n");
      exit(-1);
  }
  return user_length;
}

/* Method to be called on device open:
 * Nothing to be done, always returns success */
int do_open_or_close(struct fusd_file_info *file)
{
  DEBSTR printf("Opened/closed\n");
  return 0;
}

/* FUSD device definition */
struct fusd_file_operations echo_fops = {
  open: do_open_or_close,
  read: tranrec_read,
  write: tranrec_write,
  close: do_open_or_close
};

/* Main method creates registers device driver with fusd */
int main(int argc, char *argv[])
{
  if (fusd_register("/dev/tranrec", 0666, NULL, &echo_fops) < 0) {
    perror("register of /dev/tranrec failed");
    exit(1);
  }

  fprintf(stderr, "calling fusd_run...\n");
  fusd_run();
  return 0;
}
