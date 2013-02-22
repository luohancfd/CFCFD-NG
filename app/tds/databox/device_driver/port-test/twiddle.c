/* twiddle.c 
 * Play with IO ports.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/io.h>

int main ()
{
    int port = 0x310;
    int byte;

    if ( ioperm(port, 1, 1) ) {
	perror("ioperm open");
	exit(1);
    }
    byte = inb(port);
    printf("Input from port: 0x%x\n", byte);

    printf("Attempt write to port:\n");
    fflush(stdout);
    outb(0x20, port);
    printf("done write.\n");

    byte = inb(port);
    printf("Input from port: 0x%x\n", byte);

    printf("Done.\n");
    return 0;
}
