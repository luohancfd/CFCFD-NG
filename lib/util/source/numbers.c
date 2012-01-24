#include <stdio.h>
#include "useful.h"

int main() {
    double a, b, c, d;
    a = 0.0; a /= a;
    b = 1.0/0.0;
    c = 0.0;
    printf("Begin numbers test.\n");
    printf("a=%g\n", a);
    printf("b=%g\n", b);
    printf("c=%g\n", c);
    if ( NUMBER_EXISTS_D(a) ) {
        printf("a exists\n");
    } else {
	printf("a is not an existant number\n");
    }
    if ( NUMBER_EXISTS_D(b) ) {
        printf("b exists\n");
    } else {
	printf("b is not an existant number\n");
    }
    if ( NUMBER_EXISTS_D(c) ) {
        printf("c exists\n");
    } else {
	printf("c is not an existant number\n");
    }

    d = MAXIMUM(1.0, 2.0);
    printf("MAXIMUM(1.0, 2.0)=%g\n", d);
    d = MINIMUM(3.0, MAXIMUM(1.0,2.0));
    printf("MINIMUM(3.0, MAXIMUM(1.0,2.0))=%g\n", d);
    d = MINIMUM(3.0, MAXIMUM(1.0,5.0));
    printf("MINIMUM(3.0, MAXIMUM(1.0,5.0))=%g\n", d);
    d = MINIMUM(3.0, MAXIMUM(1.0,-5.0));
    printf("MINIMUM(3.0, MAXIMUM(1.0,-5.0))=%g\n", d);
    d = MINIMUM(-1.0, MAXIMUM(-3.0,5.0));
    printf("MINIMUM(-1.0, MAXIMUM(-3.0,5.0))=%g\n", d);

    d = MINIMUM(3.0, MAXIMUM(1.0,a));
    printf("MINIMUM(3.0, MAXIMUM(1.0,a))=%g\n", d);
    printf("a>5.0=%d\n", (a>5.0));
    printf("a<5.0=%d\n", (a<5.0));
    printf("a==5.0=%d\n", (a==5.0));

    printf("Done.");
    return 0;
}
