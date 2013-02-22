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
 * Public header function for user-space library.  This is the API
 * that user-space device drivers should write to.
 */

#ifndef __FUSD_H__
#define __FUSD_H__

#ifndef __KERNEL__
#include <sys/types.h>

__BEGIN_DECLS
#endif


#include "fusd_msg.h"

/* FUSD_NOREPLY is a special error code.  If a user-space driver
 * implementing a system call returns -FUSD_NOREPLY (note it's
 * negative!), the calling application will be blocked.  When
 * conditions enable a response to the system call (e.g. the read or
 * write has completed), the user-space driver must call the
 * fusd_return() function.  */
#define FUSD_NOREPLY  0x1000

/* FUSD defines several bitmasks for describing which channels of  
 * notification are being requested or signaled.  These flags are
 * used in the arguments and return value of the notify() callback. */
#define FUSD_NOTIFY_INPUT   0x1
#define FUSD_NOTIFY_OUTPUT  0x2
#define FUSD_NOTIFY_EXCEPT  0x4


struct fusd_file_info; /* forward decl */

typedef
struct fusd_file_operations {
  int (*open) (struct fusd_file_info *file);
  int (*close) (struct fusd_file_info *file);
  ssize_t (*read) (struct fusd_file_info *file, char *buffer, size_t length,
		   loff_t *offset);
  ssize_t (*write) (struct fusd_file_info *file, const char *buffer,
		    size_t length, loff_t *offset);
  int (*ioctl) (struct fusd_file_info *file, int request, void *data);
  int (*poll_diff) (struct fusd_file_info *file, unsigned int cached_state);
  int (*unblock) (struct fusd_file_info *file);    
} fusd_file_operations_t;


/* state-keeping structure passed to device driver callbacks */
typedef
struct fusd_file_info {
  void *device_info;		/* This is set by the library to
				 * whatever you passed to
				 * fusd_register.  Changing this in a
				 * file_operations callback has no
				 * effect. */

  void *private_data;		/* File-specific data you can change
				 * in a file_operations callback.
				 * e.g., you can set this in an open()
				 * callback, then get it in a
				 * corresponding read() callback. */

  unsigned int flags;		/* Kept synced with file->f_flags */
  pid_t pid;			/* PID of process making the request */
  uid_t uid;			/* UID of process making the request */
  gid_t gid;			/* GID of process making the request */

  /* other info might be added later, e.g. state needed to complete
     operations... */

  /* request message associated with this call */
  int fd;
  fusd_msg_t *fusd_msg;

  /* fusdnet specific fields: no not disturb! */
  void *remote_client;
  void *local_client;
  void **local_device_info;
} fusd_file_info_t;




/*************************** Library Functions ****************************/

/* fusd_register: create a device file and register callbacks for it
 *
 * Arguments:
 *
 *    name - the name of the device file, to be created wherever devfs
 *    is mounted (usually dev).  example: pass "mydevice" will create
 *    /dev/mydevice.
 *
 *    As a convenience, passing a string that starts with "/dev/" will
 *    automatically skip over that portion of the name.
 *
 *    mode - the file protections to be given to the device
 *
 *    device_info - you can provide arbitrary data that will later be
 *    passed back to your driver's callbacks in file->device_info.
 *    value has no effect on FUSD itself.
 *
 *    fops - a table of callbacks to be called for this device; see
 *    structure above.
 *
 * Return value:
 *    On failure, -1 is returned and errno is set to indicate the error.
 *
 *    On success, a valid file descriptor is returned which represents
 *    the control channel to your new device.  You should never read
 *    from or write to that control channel directcly, but you can
 *    select on it to see when it needs attention (see fusd_run and
 *    fusd_dispatch).
 */

int fusd_register(const char *name, mode_t mode, void *device_info,
		  struct fusd_file_operations *fops);



/* "simple" interface to fusd_register. */
#define fusd_simple_register(name, perms, arg, ops...) do { \
   struct fusd_file_operations f = { ops } ; \
   if (fusd_register(name, perms, arg, &f) < 0) \
      perror("warning: fusd unavailable"); \
} while(0)

/* fusd_unregister: unregister a previously registered device
 *
 * Arguments:
 *    fd - the file descriptor previously returned to you by fusd_register.
 *
 * Return value:
 *    0 on success.
 *   -1 on failure with errno set to indicate the failure.
 */
int fusd_unregister(int fd);


/* fusd_return: unblock a previously blocked system call
 * 
 * Arguments:
 *   file - the file info struct that was previously blocked
 *   retval - the return value that would have been returned by the
 *            returning system call
 *    
 * Return value:
 *   0 on success.
 *  -1 on failure with errno set to indicate the failure
 */
int fusd_return(struct fusd_file_info *file, ssize_t retval);


/*
 * fusd_destroy destroys all state associated with a fusd_file_info
 * pointer.  (It is implicitly called by fusd_return.)  If a driver
 * saves a fusd_file_info pointer by calling -FUSD_NOREPLY in order to
 * block a read, but gets a "close" request on the file before the
 * pointer is returned with fusd_return, it should be thrown away
 * using fusd_destroy.  
 */
void fusd_destroy(struct fusd_file_info *file);


/* fusd_dispatch: handles an event on a fusd file descriptor
 * 
 * Arguments:
 *   fd - the file descriptor of the device that received an event
 *    
 * Return value:
 *    None.
 *
 * Side effects:
 *    May (but may not) call a callback function originally passed to
 *    fusd_register.
 *
 *    Prints an error to stderr in case of a dispatching error.
 */
void fusd_dispatch(int fd);


/* 
 * fusd_run: convenience function that handles dispatch for all
 *           fusd devices
 *
 * No return value; runs forever.
 */
void fusd_run(void);


/*
 * fusd_fdset_add: given an FDSET and "max", add the currently valid
 * FUSD fds to the set and update max accordingly.
 */
void fusd_fdset_add(fd_set *set, int *max);


/*
 * fusd_dispatch_fdset: given an fd_set full of descriptors, call
 * fusd_dispatch on every descriptor in the set which is a valid FUSD
 * fd.
 */
void fusd_dispatch_fdset(fd_set *set);

/*
 *  FUSD status handling
 */

typedef int (* fusd_status_handler_cb)(fusd_status_t *fs, int count, void *private_data);
typedef int (* fusd_status_notify_cb)(void *private_data);

typedef struct fusd_status_context {
  char buf[4096];
  int len;                              /* initialize to 0! */
  int fd;
  int verbose:1;
  int full_mode:1;                      /* initialize to 0! */
  fusd_status_handler_cb new_data_cb;
  fusd_status_notify_cb begin_full_cb;
  fusd_status_notify_cb end_full_cb;
  void *private_data;
} fusd_status_context_t;

#define STATUS_OK     0
#define STATUS_DONE   1
#define STATUS_ERR    2
#define STATUS_FOUND  3

int fusd_status_process(fusd_status_context_t *ctx);
int fusd_status_wait(fusd_status_context_t *ctx);


/********************************************************************
 *
 *  Direct access API
 *
 *  This API enables a driver implementation to store state about a 
 *  blocked call more easily, extracting the call arguments directly
 *  with no need to store them separately.
 *
 ********************************************************************/

/* accessors */
static inline int fusd_get_call_type(struct fusd_file_info *file)
{ return file->fusd_msg->subcmd; }

static inline char * fusd_get_read_buffer(struct fusd_file_info *file)
{ return file->fusd_msg->data; }

static inline const char * fusd_get_write_buffer(struct fusd_file_info *file)
{ return (const char *)file->fusd_msg->data; }

static inline size_t fusd_get_length(struct fusd_file_info *file)
{ return (size_t)file->fusd_msg->datalen; }

static inline loff_t *fusd_get_offset(struct fusd_file_info *file)
{ return &(file->fusd_msg->parm.fops_msg.offset); }

static inline int fusd_get_ioctl_request(struct fusd_file_info *file)
{ return file->fusd_msg->parm.fops_msg.cmd; }

static inline int fusd_get_ioctl_arg(struct fusd_file_info *file)
{ return file->fusd_msg->parm.fops_msg.arg; }

static inline void * fusd_get_ioctl_buffer(struct fusd_file_info *file)
{ return (void *)file->fusd_msg->data; }

static inline int fusd_get_poll_diff_cached_state(struct fusd_file_info *file)
{ return file->fusd_msg->parm.fops_msg.cmd; }

/* returns static string representing the flagset (e.g. RWE) */
char *fusd_unparse_flags(int flags);


/*
 *  Low-level FUSD Library Helpers, used by the glib interface
 */

int fusd_open_control(void);
int fusd_build_reg_msg(fusd_msg_t *message, const char *name, mode_t mode, void *device_info);
int fusd_postreg_block(const char *name);
void fusd_dispatch_aux(int fd, int from_net, fusd_file_operations_t *use_fops, void **local_device_info);

#ifndef __KERNEL__
__END_DECLS
#endif

#endif /* __FUSD_H__ */
