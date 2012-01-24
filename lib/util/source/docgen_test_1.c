/* docgen_test_1.c */

/** These documentation comments will be extracted and written
    to a HTML file by the Awk program docgen.
    The intention is to automatically generate HTML
    documentation from the source file with a minimum of fuss.
    Of course, this preamble should be letting one know what the
    module is all about. 
 */

#include <stdio.h>

/*@function*/ double sqrit( double a ) {
    /** Accepts a double value and returns the square of that value. */
    return a * a;
} /* end function sqrit */

/* @function*/
int main ( void ) {
    /** Purpose... <BR>
        The main function, of course. 
        This program does nothing much 
        (just like 99.9% of all computer programs).
    */
    double a, b;
    printf( "Begin test program...\n" );
    a = 3.0;
    b = sqrit(a);
    printf( "%g * %g = %g\n", a, a, b );
    printf( "End test program.\n");
    return 0;
} /* end main */

/*-------------------------------------------------------*/

/* @function */ int unused_function( 
      int a, 
      double b 
   ) 
/** An unused function with its prototype having a different layout.
The function is otherwise unused by the program.
*/
{
    return 0;
} /* end unused_function */

/*-------------------------------------------------------*/

/**
    <HR>
    <H4>Postscript...</H4>
    <P>
    This is a section of documentation with its
    own HTML tags.
    </P> 
*/

