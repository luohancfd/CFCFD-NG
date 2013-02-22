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
 * Jeremy Elson <jelson@circlemud.org>
 * Copyright (c) Sensoria Corporation 2001
 *
 * Private header file used by the Linux Kernel Module
 *
 * $Id: kfusd.h,v 1.51 2004/04/25 03:30:23 girod Exp $
 */

#ifndef __KFUSD_H__
#define __KFUSD_H__

#include "fusd_msg.h"

/* magic numbers for structure checking; unique w.r.t
 * /usr/src/linux/Documentation/magic-number.txt */
#define FUSD_DEV_MAGIC      0x8b43a124
#define FUSD_FILE_MAGIC     0x613aa8ff

/* number of times each device can be opened simultaneously */
#define MIN_FILEARRAY_SIZE  8  /* initialize allocation */
#define MAX_FILEARRAY_SIZE  1024 /* maximum it can grow to */

/* maximum read/write size we're willing to service */
#define MAX_RW_SIZE         (1024*128)

/* maximum number of major and minor devices */
#define MAX_MAJORS 256
#define MAX_MINORS 256

/********************** Structure Definitions *******************************/

/* Container for a fusd msg */
typedef struct fusd_msgC_s_t fusd_msgC_t;

struct fusd_msgC_s_t {
  fusd_msg_t fusd_msg;		/* the message itself */
  fusd_msgC_t *next;		/* pointer to next one in the list */

  /* 1-bit flags */
  unsigned int peeked:1;	/* has the first half of this been read? */
};

/* magical forward declarations to break the circular dependency */
struct fusd_dev_t_s;
typedef struct fusd_dev_t_s fusd_dev_t;

/* state kept per opened file (i.e., an instance of a device) */
typedef struct {
  /* general state management */
  int magic;			/* magic number for sanity checking */
  fusd_dev_t *fusd_dev;		/* fusd device associated with this file */
  long fusd_dev_version;	/* version number of fusd device */
  void *private_data;		/* the user's private data (we ignore it) */
  struct file *file;		/* kernel's file pointer for this file */
  pid_t pid;			/* PID of the last user of this FD */
  int index;			/* our index in our device's file array */
  struct semaphore file_sem;	/* Semaphore for file structure */
  int cached_poll_state;	/* Latest result from a poll diff req */
  int last_poll_sent;		/* Last polldiff request we sent */
  int poll_failure;		/* errno from a failed poll_diff */
  int in_select;		/* for debug only (/dev/fusd/status) */

  /* structures used for messaging */
  wait_queue_head_t file_wait;	/* Wait on this for a user->kernel msg */
  wait_queue_head_t poll_wait;  /* Given to kernel for poll() queue */
  long transid_outstanding;	/* transid of msg we are waiting for */
  int subcmd_outstanding;	/* subcmd of msg we are waiting for */
  pid_t sys_restarting_pid;	/* PID we just returned -ERESTARTSYS to */
  int sys_restarting_subcmd;	/* subcmd of restarting syscall */
  fusd_msg_t *msg_in;		/* A reply we've just received */
} fusd_file_t;


/* state kept per device registered under fusd */
struct fusd_dev_t_s {
  int magic;			/* Magic number for sanity checking */
  long version;			/* Instance number of this device */
  int last_version;             /* Indicates version of last change */
  int zombie;			/* Is the device dead? */
  pid_t pid;			/* PID of device driver */
  char *name;			/* Name of the device under devfs (/dev) */
  void *private_data;		/* User's private data */
#ifdef CONFIG_DEVFS_FS 
#if LINUX_VERSION_CODE < KERNEL_VERSION(2,6,0)
  devfs_handle_t handle;	/* The devfs-provided handle */
#else
  /* handle not required for >= 2.6.0 */
#endif
#endif
  dev_t device;                 /* the device (major/minor) */
  int hash;                     /* the hash of the device name */

  fusd_file_t **files;		/* Array of this device's open files */
  int array_size;		/* Size of the array pointed to by 'files' */
  int num_files;		/* Number of array entries that are valid */
  int open_in_progress;		/* File is referencing this struct,
                                   but not yet part of the file array */
  /* messaging */
  fusd_msgC_t *msg_head;	/* linked list head for message queue */
  fusd_msgC_t *msg_tail;	/* linked list tail for message queue */

  /* synchronization */
  wait_queue_head_t dev_wait;	/* Wait queue for kernel->user msgs */
  struct semaphore dev_sem;	/* Sempahore for device structure */

  /* pointer to allow a dev to be placed on a dev_list */
  struct list_head devlist;

  /* Flag indicating whether MAKEDEV is done */
  uint8_t makedev_done;
};
  

/**** Function Prototypes ****/

STATIC int maybe_free_fusd_dev(fusd_dev_t *fusd_dev);

STATIC int find_fusd_file(fusd_dev_t *fusd_dev, fusd_file_t *fusd_file);
STATIC int free_fusd_file(fusd_dev_t *fusd_dev, fusd_file_t *fusd_file);

STATIC int fusd_fops_call_send(fusd_file_t *fusd_file_arg,
			       fusd_msg_t *fusd_msg);
STATIC int fusd_fops_call_wait(fusd_file_t *fusd_file_arg,
			       fusd_msg_t **reply);
STATIC void fusd_fops_call_done(fusd_file_t *fusd_file);

STATIC void fusd_forge_close(fusd_msg_t *msg, fusd_dev_t *fusd_dev);





/**** Utility functions & macros ****/

#ifdef CONFIG_FUSD_USE_WAKEUPSYNC
#define WAKE_UP_INTERRUPTIBLE_SYNC(x) wake_up_interruptible_sync(x)
#else
#define WAKE_UP_INTERRUPTIBLE_SYNC(x) wake_up_interruptible(x)
#endif /* CONFIG_FUSD_USE_WAKEUPSYNC */

#ifdef CONFIG_FUSD_DEBUG
static void rdebug_real(char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));

#define RDEBUG(message_level, args...) do { \
   if (fusd_debug_level >= message_level) rdebug_real(args); \
} while(0)
#else
#define RDEBUG(message_level, args...)
#endif /* CONFIG_FUSD_DEBUG */


#define ZOMBIE(fusd_dev)  ((fusd_dev)->zombie)


#define GET_FUSD_DEV(candidate, fusd_dev) do { \
  fusd_dev = candidate; \
  if (fusd_dev == NULL || fusd_dev->magic != FUSD_DEV_MAGIC) \
        goto invalid_dev; \
} while (0)

#define GET_FUSD_FILE_AND_DEV(candidate, fusd_file, fusd_dev) do { \
  fusd_file = candidate; \
  if (fusd_file == NULL || fusd_file->magic != FUSD_FILE_MAGIC) \
     goto invalid_file; \
  GET_FUSD_DEV(fusd_file->fusd_dev, fusd_dev); \
  if (fusd_dev->version != fusd_file->fusd_dev_version) \
    goto invalid_file; \
} while (0)


#define LOCK_FUSD_DEV(fusd_dev) \
  do { down(&fusd_dev->dev_sem); \
  if (ZOMBIE(fusd_dev)) { up(&fusd_dev->dev_sem); goto zombie_dev; } \
 } while (0)

/* rawlock does not do a zombie check */
#define RAWLOCK_FUSD_DEV(fusd_dev) \
  do { down(&fusd_dev->dev_sem); } while (0)

#define UNLOCK_FUSD_DEV(fusd_dev) \
  do { up(&fusd_dev->dev_sem); } while (0)


#define LOCK_FUSD_FILE(fusd_file) \
  do { down(&fusd_file->file_sem); \
 } while (0)

#define UNLOCK_FUSD_FILE(fusd_file) \
  do { up(&fusd_file->file_sem); } while (0)

#define FREE_FUSD_MSGC(fusd_msgc) do { \
   if ((fusd_msgc)->fusd_msg.data != NULL) KFREE(fusd_msgc->fusd_msg.data); \
   KFREE(fusd_msgc); \
} while (0)

#define NAME(fusd_dev) ((fusd_dev)->name == NULL ? \
			"<noname>" : (fusd_dev)->name)

#ifdef CONFIG_FUSD_MEMDEBUG
static int  fusd_mem_init(void);
static void fusd_mem_cleanup(void);
static void fusd_mem_add(void *ptr, int line, int size);
static void fusd_mem_del(void *ptr);
static void *fusd_kmalloc(size_t size, int type, int line);
static void fusd_kfree(void *ptr);
static void *fusd_vmalloc(size_t size, int line);
static void fusd_vfree(void *ptr);
# define KMALLOC(size, type) fusd_kmalloc(size, type, __LINE__)
# define KFREE(ptr) fusd_kfree(ptr)
# define VMALLOC(size) fusd_vmalloc(size, __LINE__)
# define VFREE(ptr) fusd_vfree(ptr)
#else /* no memory debugging */
# define KMALLOC(size, type) kmalloc(size, type)
# define KFREE(ptr) kfree(ptr)
# define VMALLOC(size) vmalloc(size)
# define VFREE(ptr) vfree(ptr)
#endif /* CONFIG_FUSD_MEMDEBUG */



/* Functions like this should be in the kernel, but they are not.  Sigh. */
#ifdef CONFIG_SMP

DECLARE_MUTEX(atomic_ops);

static __inline__ int atomic_inc_and_ret(int *i)
{
  int val;

  down(&atomic_ops);
  val = (++(*i));
  up(&atomic_ops);
  return val;
}
#else
static __inline__ int atomic_inc_and_ret(int *i)
{
  return (++(*i));
}
#endif


#endif /* __KFUSD_H__ */
