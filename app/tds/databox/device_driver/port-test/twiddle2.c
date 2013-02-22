/* twiddle2.c 
 * Play with IO ports via /dev/port.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/io.h>

int main ()
{
    int port = 0x310;
    int byte = 0;
    int fp;

    fp = open("/dev/port", O_RDWR);
    if ( fp == -1 ) {
	perror("open /dev/port");
	exit(1);
    }

    lseek(fp, port, SEEK_SET);
    read(fp, &byte, 1);
    printf("Input from port: 0x%x\n", byte);

    printf("Attempt write to port:\n");
    fflush(stdout);
    lseek(fp, port, SEEK_SET);
    byte = 0x20;
    write(fp, &byte, 1);
    printf("done write.\n");

    lseek(fp, port, SEEK_SET);
    read(fp, &byte, 1);
    printf("Input from port: 0x%x\n", byte);

    printf("Done.\n");
    return 0;
}
