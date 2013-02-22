/*
 *
 * Copyright (c) 2003 The Regents of the University of California.  All 
 * rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Neither the name of the University nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS''
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
 

/*
 * FUSD: the Framework for User-Space Devices
 *
 * Defines the interface between the kernel module and userspace library.
 *
 */

#ifndef __KERNEL__
#include <inttypes.h>
#endif

#ifndef __FUSD_MSG_H__
#define __FUSD_MSG_H__

/* filenames */
#define DEFAULT_DEV_ROOT           "/dev/"
#define FUSD_CONTROL_FILENAME      "fusd/control"
#define FUSD_STATUS_FILENAME       "fusd/status"
#define NET_BASE_DIR               "net"

#define FUSD_NET_BASE_DIR          "fusd/net/"
#define FUSD_NET_LOG_FILENAME      FUSD_NET_BASE_DIR "log"
#define FUSD_NET_COMMAND_FILENAME  FUSD_NET_BASE_DIR "command"
#define FUSD_NET_STATUS_FILENAME   FUSD_NET_BASE_DIR "status"
#define FUSD_NET_INT_FILENAME      FUSD_NET_BASE_DIR ".int"

#define FUSD_CONTROL_DEVNAME       DEFAULT_DEV_ROOT FUSD_CONTROL_FILENAME
#define FUSD_STATUS_DEVNAME        DEFAULT_DEV_ROOT FUSD_STATUS_FILENAME
#define FUSD_NET_LOG_DEVNAME       DEFAULT_DEV_ROOT FUSD_NET_LOG_FILENAME
#define FUSD_NET_COMMAND_DEVNAME   DEFAULT_DEV_ROOT FUSD_NET_COMMAND_FILENAME
#define FUSD_NET_STATUS_DEVNAME    DEFAULT_DEV_ROOT FUSD_NET_STATUS_FILENAME
#define FUSD_NET_INT_DEVNAME       DEFAULT_DEV_ROOT FUSD_NET_INT_FILENAME

/* ioctl number to tell FUSD status device to return binary info */
#define FUSD_STATUS_USE_BINARY     _IO('F', 100)

/* ioctl number to support non-devfs FUSD daemon */
#define FUSD_STATUS_I_AM_FUSDD     _IO('F', 101)
/* ioctl number to query whether this is a DevFS system */
#define FUSD_STATUS_NO_DEVFS       _IO('F', 102)

/* constants */
#define FUSD_MAX_NAME_LENGTH       47 /* 47, to avoid expanding union size */


/* commands */
#define FUSD_REGISTER_DEVICE       0 /* device registration */
#define FUSD_UNREGISTER_DEVICE     1 /* device unregistration */

/* these two must have successive numbers */
#define FUSD_FOPS_CALL             2 /* synchronous round-trip call: request */
#define FUSD_FOPS_REPLY            (FUSD_FOPS_CALL + 1)

/* these two must have successive numbers */
#define FUSD_FOPS_NONBLOCK         4 /* call that does not block for a reply */
#define FUSD_FOPS_NONBLOCK_REPLY   (FUSD_FOPS_NONBLOCK + 1)

#define FUSD_FOPS_CALL_DROPREPLY   6 /* call that doesn't want a reply */

/* fusdnet internal message types */
#define FUSD_NET_MESSAGE           10
#define FUSD_NET_SUBOP_CONNECT     110

/* subcommands */
#define FUSD_OPEN                  100
#define FUSD_CLOSE                 101
#define FUSD_READ                  102
#define FUSD_WRITE                 103
#define FUSD_IOCTL                 104
#define FUSD_POLL_DIFF             105
#define FUSD_UNBLOCK               106

/* other constants */
#define FUSD_MSG_MAGIC      0x7a6b93cf

/* user->kernel: register a device */
typedef struct {
  char name[FUSD_MAX_NAME_LENGTH+1];
  mode_t mode;
  void *device_info;
} register_msg_t;


/* kernel->user: fops request message (common data) */
typedef struct {
  pid_t pid;
  uid_t uid;
  gid_t gid;
  unsigned int flags;		/* flags from file struct */
  void *device_info;		/* device info */
  void *private_info;		/* file info */

  /* parameters and return values for various calls.  should be a
   * union but it just makes things too complex and doesn't save all
   * that much memory anyway */
  ssize_t retval;
  size_t length;
  loff_t offset;
  unsigned int cmd;  /* ioctl cmd, poll_diff cached_state */
  unsigned long arg; /* ioctl */

  /* the following are cookies that have meaning internal to the kernel
   * but must be returned, untouched, by userspace */
  void *fusd_file;
  long transid;
  int hint;
} fops_msg_t;


/* the message struct written to FUSD control channel */
typedef struct {
  int magic;
  short int cmd;
  short int subcmd;

  char *data;  /* yes, it's slightly inefficient to push this useless
		* pointer between user and kernel space, but it makes
		* it much easier to have a pointer available in this
		* structure that both the kernel and userlib can make
		* their own use of. */
  int datalen;
  union {
    register_msg_t register_msg; /* device registration (U->K) */
    fops_msg_t fops_msg;	/* U->K and K->U fops messages */
  } parm;
} fusd_msg_t;


/* structure read from FUSD binary status device */
typedef struct {
  char name[FUSD_MAX_NAME_LENGTH+1];
  int zombie:1;
  int unmade:1;
  int full_mode:1;
  int reserved:29;
  pid_t pid;
  int num_open;
  uint16_t major;
  uint16_t minor;
} fusd_status_t;


/*  
 *  non-devfs minor numbers
 */

#define FUSD_DEV_MAJOR       242
#define FUSD_CONTROL_MINOR     0
#define FUSD_STATUS_MINOR      1
#define FUSD_DEV_MINOR         2


/*
 *  FUSDnet Port
 */

#define FUSD_NET_PORT       1212

#endif /* __FUSD_MSG_H__ */
