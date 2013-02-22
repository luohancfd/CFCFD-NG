/** \file   dbox_driver.c
 * \brief   User-space device driver for the BCD databox.
 * \author  Zane Smith and Peter Jacobs
 * 
 * \usage   sudo ./dbox_driver [--help] \
 *                             [--buffer=0xd8000] \
 *                             [--timebase=0x310] \
 *                             [--cardbase=0x320]
 *
 * This driver exposes the control and status registers of
 * the BCD databox as files in /dev/databox/.
 * Communication with the databox registers is typically by reading 
 * and writing strings to these special files.
 * An exception is the databox data buffer which is read as 
 * one big (binary) string of 16kB characters.  
 *
 * \note Note that integer values written to the register files 
 *       are specified in hexadecimal format.
 *
 * \note This driver requires FUSD to run.
 *
 * \note It is convenient to set the owner of the executable
 *       dbox_driver as root and to then change the mode to setuid
 *       so that the driver runs with root permissions without
 *       having to have the user in the sudo list.
 *
 * \version 1.0, October 2004
 *          Rebuild Zane's prototype driver.
 * \version 1.01, December 2004
 *          Go back to using int rather than short int for 
 *          holding byte values.  The %x conversion specifier
 *          caused some grief with short ints.
 */

#include <sys/io.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "fusd.h"

#define MIN(x,y) ((x)<(y) ? (x) : (y))

#define DATA_BUFFER_SIZE   16384

/*
 * The following default values should suit most of the databoxes. 
 * The values can be overridden via the command line.
 */
static int DATA_BUFFER_ADDR = 0xd8000;
static int CARD_BASE_ADDR = 0x320;
static int TIMEBASE_BASE_ADDR = 0x310;

static int debug = 1;

static char hex_str_buf[5];
static int most_recently_selected_card = 0;

/*--------------------------------------------------------------*/
/* Convenience functions */

/** \brief Writes byte_value to the I/O-port at writeaddress. 
 */
void do_out_port(int port_address, int byte_value)
{
    int fd;
    fd = open("/dev/port", O_RDWR);
    if (fd == -1) {
	perror("open /dev/port");
	exit(-1);
    }
    if ( lseek ( fd, port_address, 0) != port_address ) {
	perror ("lseek: ");
	close (fd);
	exit (-1);
    }
    if ( write (fd, &byte_value, 1) != 1) {
	perror ("write: ");
	close (fd);
	exit (-2);
    }
    if (debug) printf("Set port 0x%x to value 0x%x\n", port_address, byte_value);
    close (fd);
}

/** \brief Copies the byte in the I/O-port into container.
 */
void do_in_port(int port_address, int *container)
{
    int fd;
    *container = 0;  /* Clear out all of the bits. */
    fd = open("/dev/port", O_RDONLY);
    if (fd == -1) {
	perror("open /dev/port");
	exit(-1);
    }
    if ( lseek ( fd, port_address, 0) != port_address ) {
	perror ("lseek: ");
	close (fd);
	exit (-1);
    }
    if ( read (fd, container, 1) != 1) {
	perror ("read: ");
	close (fd);
	exit (-2);
    }
    if (debug) printf("Read port 0x%x : value=0x%x\n", port_address, *container);
    close (fd);
}


/** \brief Always gives us two characters representing our number in hex.
 */
char *to_hex_string( short int bv )
{
    (bv <= 0x0F) ? sprintf(hex_str_buf,"0%x",bv) : sprintf(hex_str_buf,"%x",bv);
    return hex_str_buf;
}

/*-----------------------------------------------------------------*/

/* Device service functions follow... */


int trigger_status_open(struct fusd_file_info *file)
{
    if (debug) printf("Open trigger_status.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

int trigger_status_close(struct fusd_file_info *file)
{
    if (debug) printf("Close trigger_status.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

/** \brief Read the pair of trigger-unit status registers
 *         and return the values as hexadecimal (in a string). 
 */
int trigger_status_read(struct fusd_file_info *file, 
			  char *user_buffer,
			  size_t user_length, 
			  loff_t *offset)
{
    int  tA = 0;
    int  tB = 0;
    char buf[128];
    int  port_address;

    /* Subsequent reads (without closing and reopening) return EOF. */
    if (*offset > 0) return 0;

    /* Do the real work. */
    port_address = TIMEBASE_BASE_ADDR;
    do_in_port(port_address, &tA);
    port_address = TIMEBASE_BASE_ADDR + 1;
    do_in_port(port_address, &tB);
    sprintf(buf,"0x%x 0x%x",tA, tB);
    if (debug) printf("trigger_status_read: result=%s\n",buf);

    /* Now, deal with the file-like activities. */
    user_length = MIN(user_length,strlen(buf));
    memcpy(user_buffer, buf, user_length);
    *offset += user_length;
    return user_length;
} /* end trigger_status_read() */

/*----------------------------------------------------------------*/

int trigger_control_open(struct fusd_file_info *file)
{
    if (debug) printf("Open trigger_control.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

int trigger_control_close(struct fusd_file_info *file)
{
    if (debug) printf("Close trigger_control.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

/** \brief Write a single byte to the trigger_control register.
 *
 * The data is supplied as a hexadecimal string representing 
 * the integer value. 
 */
ssize_t trigger_control_write(struct fusd_file_info *file, 
			       const char *user_buffer,
			       size_t user_length, 
			       loff_t *offset)
{
    int value = 0;
    int port_address;
    sscanf(user_buffer, "%x", &value);
    port_address = TIMEBASE_BASE_ADDR;
    if (debug) printf("trigger_control_write: value=0x%x\n", value);
    do_out_port(port_address,value);
    return user_length;
} /* end trigger_control_write() */

/*----------------------------------------------------------------*/

int timebase_status_open(struct fusd_file_info *file)
{
    if (debug) printf("Open timebase_status.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

int timebase_status_close(struct fusd_file_info *file)
{
    if (debug) printf("Close timebase_status.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

/** \brief Read the three pair of timebase status registers
 *         and return the values as hexadecimal (in a string). 
 */
int timebase_status_read(struct fusd_file_info *file, 
			  char *user_buffer,
			  size_t user_length, 
			  loff_t *offset)
{
    int  u, t1A, t1B, t2A, t2B, t3A, t3B;
    char buf[128];
    int  port_address;

    /* Subsequent reads (without closing and reopening) return EOF. */
    if (*offset > 0) return 0;

    /* Do the real work. */
    u = 1;
    port_address = TIMEBASE_BASE_ADDR + 2*u;
    do_in_port(port_address, &t1A);
    port_address = TIMEBASE_BASE_ADDR + 2*u + 1;
    do_in_port(port_address, &t1B);
    u = 2;
    port_address = TIMEBASE_BASE_ADDR + 2*u;
    do_in_port(port_address, &t2A);
    port_address = TIMEBASE_BASE_ADDR + 2*u + 1;
    do_in_port(port_address, &t2B);
    u = 3;
    port_address = TIMEBASE_BASE_ADDR + 2*u;
    do_in_port(port_address, &t3A);
    port_address = TIMEBASE_BASE_ADDR + 2*u + 1;
    do_in_port(port_address, &t3B);
    sprintf(buf,"0x%x 0x%x 0x%x 0x%x 0x%x 0x%x",t1A, t1B, t2A, t2B, t3A, t3B);
    if (debug) printf("timebase_status_read: result=%s\n",buf);

    /* Now, deal with the file-like activities. */
    user_length = MIN(user_length,strlen(buf));
    memcpy(user_buffer, buf, user_length);
    *offset += user_length;
    return user_length;
} /* end timebase_status_read() */

/*----------------------------------------------------------------*/

int card_control_open(struct fusd_file_info *file)
{
    if (debug) printf("Open card_control.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

int card_control_close(struct fusd_file_info *file)
{
    if (debug) printf("Close card_control.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

/** \brief Write a single byte to the card_control register.
 *
 * The data is supplied as a string representing the hexadecimal
 * (integer) number for the card and a hexadecimal value to write
 * to the register. 
 * A card number of 0x0 indicates that we want to write to all cards 
 * simultaneously.
 */
ssize_t card_control_write(struct fusd_file_info *file, 
			   const char *user_buffer,
			   size_t user_length, 
			   loff_t *offset)
{
    int value = 0;
    int card = 0;
    int port_address;
    sscanf(user_buffer, "%x %x", &card, &value);
    if ( card > 7 ) card = 7;  /* So we don't go too far in port space. */
    if (card == 0) {
	/* Want to write to all control registers simultaneously. */
	port_address = CARD_BASE_ADDR - 1;
    } else {
	/* Write to one card register, alone. */
	port_address = CARD_BASE_ADDR + 2 * card;
    }
    do_out_port(port_address,value);
    if (debug) printf("card_control_write: port=0x%x value=0x%x\n", 
		      port_address, value);
    return user_length;
} /* end card_control_write() */

/*----------------------------------------------------------------*/

int card_status_open(struct fusd_file_info *file)
{
    if (debug) printf("Open card_status.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

int card_status_close(struct fusd_file_info *file)
{
    if (debug) printf("Close card_status.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

/** \brief Write a single byte to the *collective* card_control register
 *         to select the card for subsequent read of the status register.
 *
 * The data is expected in the form of a string of characters 
 * specifying n c and r as hexadecimal integers.
 */
ssize_t card_status_write(struct fusd_file_info *file, 
			   const char *user_buffer,
			   size_t user_length, 
			   loff_t *offset)
{
    int n = 0;  /* card number (bits 0-3)    */
    int c = 0;  /* channel number (bits 4,5) */
    int r = 0;  /* reset bit  (bit 7)        */
    int value;
    int port_address = CARD_BASE_ADDR - 1;
    sscanf(user_buffer, "%x %x %x", &n, &c, &r);
    if (debug) 
	printf("card_status_write: n=0x%x, c=0x%x, r=0x%x\n", n, c, r);
    value = (r << 7) | (c << 4) | (n);
    if (debug) printf("card_status_write: value=0x%x\n", value);
    do_out_port(port_address,value);
    most_recently_selected_card = n;
    return user_length;
} /* end card_status_write() */

/** \brief Read the pair of card status registers
 *         and return the values as hexadecimal (in a string). 
 */
int card_status_read(struct fusd_file_info *file, 
		     char *user_buffer,
		     size_t user_length, 
		     loff_t *offset)
{
    int  n, tA, tB;
    char buf[128];
    int  port_address;

    /* Subsequent reads (without closing and reopening) return EOF. */
    if (*offset > 0) return 0;

    /* Do the real work. */
    n = most_recently_selected_card;
    port_address = CARD_BASE_ADDR + 2*n;
    do_in_port(port_address, &tA);
    port_address = CARD_BASE_ADDR + 2*n + 1;
    do_in_port(port_address, &tB);
    sprintf(buf,"0x%x 0x%x",tA, tB);
    if (debug) printf("card_status_read: card_n=%d result=%s\n", n, buf);

    /* Now, deal with the file-like activities. */
    user_length = MIN(user_length,strlen(buf));
    memcpy(user_buffer, buf, user_length);
    *offset += user_length;
    return user_length;
} /* end card_status_read() */

/*----------------------------------------------------------------*/

int reset_cards_and_arm_open(struct fusd_file_info *file)
{
    if (debug) printf("Open reset_card_and_arm.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

int reset_cards_and_arm_close(struct fusd_file_info *file)
{
    if (debug) printf("Close reset_card_and_arm.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

/** \brief Reset a set of cards and arm all trigger units.
 *
 * The data is expected in the form of a string of characters 
 * specifying 7 hexadecimal values.  
 * A value of 1 for a particular card indicates that we want to
 * reset that card, otherwise it is left alone. 
 * (A card may not be present in every slot.) 
 */
ssize_t reset_cards_and_arm_write(struct fusd_file_info *file, 
				  const char *user_buffer,
				  size_t user_length, 
				  loff_t *offset)
{
    int n, card;
    int r[8];
    int port_address;
    n = sscanf(user_buffer, "%x %x %x %x %x %x %x", 
	       &r[1], &r[2], &r[3], &r[4], &r[5], &r[6], &r[7]);
    if (debug) 
	printf("reset_cards_and arm_write: received %d integers\n", n);
    for (card = 1; card <= n; ++card) {
	if (r[card] == 1) {
	    if (debug) printf("reset_cards_and_arm_write: card[%d]\n", card);
	    port_address = CARD_BASE_ADDR + 2 * card;
	    do_out_port(port_address,0x80);  /* bit-7 high */
	    do_out_port(port_address,0x00);  /* clear */
	}
    }
    /* Now that the cards are synchronised, arm the box. */
    port_address = TIMEBASE_BASE_ADDR;
    do_out_port(port_address,0x20); /* bit-5 in trigger_control register */
    return user_length;
} /* end reset_cards_and_arm_write() */

/*----------------------------------------------------------------*/

int data_buffer_open(struct fusd_file_info *file)
{
    if (debug) printf("Open data_buffer.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

int data_buffer_close(struct fusd_file_info *file)
{
    if (debug) printf("Close data_buffer.\n");
    /* Nothing to be done, always returns success. */
    return 0;
}

/** \brief Write a single byte to the *collective* card_control register
 *         to select the card for subsequent read of the status register.
 *
 * The data is expected in the form of a string of characters 
 * specifying n and c as hexadecimal integers.
 * A hexadecimal (string) value is returned.
 */
ssize_t data_buffer_write(struct fusd_file_info *file, 
			  const char *user_buffer,
			  size_t user_length, 
			  loff_t *offset)
{
    int n = 0;  /* card number (bits 0-3)    */
    int c = 0;  /* channel number (bits 4,5) */
    int r = 0;  /* reset bit  (bit 7)        */
    int value;
    int port_address = CARD_BASE_ADDR - 1;
    sscanf(user_buffer, "%x %x", &n, &c);
    if (debug) printf("data_buffer_write: n=%d, c=%d\n", n, c);
    value = (r << 7) | (c << 4) | (n);
    if (debug) printf("data_buffer_write: value=0x%x\n", value);
    do_out_port(port_address,value);
    most_recently_selected_card = n;
    return user_length;
} /* end data_buffer_write() */

/** \brief Read the data buffer of the selected card
 *         and return the values as a string of raw data bytes.
 *
 * It is assumed that a previous write to this file has selected
 * the appropriate card and channel. 
 */
int data_buffer_read(struct fusd_file_info *file, 
		     char *user_buffer,
		     size_t user_length, 
		     loff_t *offset)
{
    int  fd;
    unsigned char *isa_mem;

    if (*offset >= DATA_BUFFER_SIZE) return 0; /* EOF if we have finished */

    fd = open("/dev/mem", O_RDONLY);
    if (fd == -1) {
	perror("open /dev/mem");
	exit(-1);
    }
    isa_mem = (unsigned char *)mmap((void *)0, DATA_BUFFER_SIZE, 
				    PROT_READ, MAP_SHARED, fd, 
				    DATA_BUFFER_ADDR);
    if ((long)isa_mem == -1) {
	perror("memmap /dev/mem");
	exit(-1);
    }
    if (debug) printf("data_buffer_read: memory mapped\n");
    user_length = MIN(user_length, DATA_BUFFER_SIZE - *offset);
    memcpy(user_buffer, isa_mem + *offset, user_length);
    close(fd);
    *offset += user_length;
    return user_length;
} /* end data_buffer_read() */

/*----------------------------------------------------------------*/

/* FUSD device definitions. */

struct fusd_file_operations trigger_status_fops = {
    open:  trigger_status_open,
    read:  trigger_status_read,
    close: trigger_status_close
};

struct fusd_file_operations trigger_control_fops = {
    open:  trigger_control_open,
    write: trigger_control_write,
    close: trigger_control_close
};

struct fusd_file_operations timebase_status_fops = {
    open:  timebase_status_open,
    read:  timebase_status_read,
    close: timebase_status_close
};

struct fusd_file_operations card_control_fops = {
    open:  card_control_open,
    write: card_control_write,
    close: card_control_close
};

struct fusd_file_operations card_status_fops = {
    open:  card_status_open,
    write: card_status_write,
    read:  card_status_read,
    close: card_status_close
};

struct fusd_file_operations reset_cards_and_arm_fops = {
    open:  reset_cards_and_arm_open,
    write: reset_cards_and_arm_write,
    close: reset_cards_and_arm_close
};

struct fusd_file_operations data_buffer_fops = {
    open:  data_buffer_open,
    write: data_buffer_write,
    read:  data_buffer_read,
    close: data_buffer_close
};

/*----------------------------------------------------------------*/

int print_usage( void )
{
    printf("sudo ./dbox_driver [--help]\n");
    printf("                   [--buffer=0xd8000]\n");
    printf("                   [--timebase=0x310]\n");
    printf("                   [--cardbase=0x320]\n");
    printf("\n");
    printf("Note that sudo is not necessary if the setuid has been set.\n");
    return 0;
}

/*----------------------------------------------------------------*/

/** \brief Main function registers the call-back functions with FUSD.
 */
int main(int argc, char *argv[])
{
    int i, handle;
    int sp1, sp2;
    printf("BCD databox device driver (Zane Smith and Peter Jacobs, 2004): \n");

    /* Check for command-line arguments. */
    i = 1;
    while ( i < argc ) {
        if ( strcmp(argv[i],"--help") == 0 ) {
            print_usage();
            exit(0);
        } else if ( strncmp(argv[i],"--buffer",8) == 0 ) {
            sscanf(argv[i], "--buffer=0x%x", &DATA_BUFFER_ADDR);
            ++i;
        } else if ( strncmp(argv[i],"--timebase",10) == 0 ) {
            sscanf(argv[i], "--timebase=0x%x", &TIMEBASE_BASE_ADDR);
            ++i;
        } else if ( strncmp(argv[i],"--cardbase",10) == 0 ) {
            sscanf(argv[i], "--cardbase=0x%x", &CARD_BASE_ADDR);
            ++i;
        } else {
            printf("Unknown command-line option.\n");
            print_usage();
            exit(0);
        }
    }
    printf("    DATA_BUFFER_ADDR   = 0x%x\n", DATA_BUFFER_ADDR);
    printf("    TIMEBASE_BASE_ADDR = 0x%x\n", TIMEBASE_BASE_ADDR);
    printf("    CARD_BASE_ADDR     = 0x%x\n", CARD_BASE_ADDR);

    /*
     * Check for the presence of the databox by looking at a couple
     * of the nominated status ports.
     * We should see all high bits if the databox is not actively
     * pulling any of the bits down.
     */
    do_in_port(TIMEBASE_BASE_ADDR, &sp1);
    do_in_port(CARD_BASE_ADDR, &sp2);
    if ( sp1 == 0xFF && sp2 == 0xFF ) {
        printf("WARNING: Databox is not present at nominated address.\n");
        printf("         Will continue anyway.\n");
    }

    printf("Begin registration...\n" );
    handle = fusd_register("/dev/databox/trigger_status", 0666, NULL, 
			   &trigger_status_fops);
    if ( handle < 0) {
	perror("databox error: registration of trigger_status failed");
	exit(1);
    } else {
	printf("    created /dev/databox/trigger_status\n" );
    }

    handle = fusd_register("/dev/databox/trigger_control", 0666, NULL, 
			   &trigger_control_fops);
    if ( handle < 0) {
	perror("databox error: registration of trigger_control failed");
	exit(1);
    } else {
	printf("    created /dev/databox/trigger_control\n" );
    }

    handle = fusd_register("/dev/databox/timebase_status", 0666, NULL, 
			   &timebase_status_fops);
    if ( handle < 0) {
	perror("databox error: registration of timebase_status failed");
	exit(1);
    } else {
	printf("    created /dev/databox/timebase_status\n" );
    }

    handle = fusd_register("/dev/databox/card_control", 0666, NULL, 
			   &card_control_fops);
    if ( handle < 0) {
	perror("databox error: registration of card_control failed");
	exit(1);
    } else {
	printf("    created /dev/databox/card_control\n" );
    }

    handle = fusd_register("/dev/databox/card_status", 0666, NULL, 
			   &card_status_fops);
    if ( handle < 0) {
	perror("databox error: registration of card_status failed");
	exit(1);
    } else {
	printf("    created /dev/databox/card_status\n" );
    }

    handle = fusd_register("/dev/databox/reset_cards_and_arm", 0666, NULL, 
			   &reset_cards_and_arm_fops);
    if ( handle < 0) {
	perror("databox error: registration of reset_cards_and_arm failed");
	exit(1);
    } else {
	printf("    created /dev/databox/reset_cards_and_arm\n" );
    }

    handle = fusd_register("/dev/databox/data_buffer", 0666, NULL, 
			   &data_buffer_fops);
    if ( handle < 0) {
	perror("databox error: registration of data_buffer failed");
	exit(1);
    } else {
	printf("    created /dev/databox/data_buffer\n" );
    }

    printf("Successfully registered device driver.\n");
    printf("Now, you can interact with the device files.\n");
    printf("Calling fusd_run...(Press Ctrl-C to stop this driver.)\n");
    fusd_run();
    return 0;
} /* end main() */
