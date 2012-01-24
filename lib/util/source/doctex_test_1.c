/* doctex_test_1.c */

/** These documentation comments will be extracted and written
    to a file as fragments of LaTeX by the Awk program doctex.
    The intention is to automatically generate documentation 
    with LaTeX markup from the source file with a minimum of fuss.
    Of course, this preamble should be letting one know what the
    module is all about. 
 */

#include <stdio.h>

/****f*/ double sqrit( double a ) {
    /** Accepts a double value and returns the square of that value. */
    return a * a;
} /* end function sqrit */

//****f*
int main ( void ) {
    /** This is where it all begins...
     *
     *  This is the \textit{main} function, of course. 
     *  This program does nothing much 
     *  (just like 99.9\% of all computer programs).
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
 *  The function is otherwise unused by the program.
 */
{
    return 0;
} /* end unused_function */

/*-------------------------------------------------------*/

/** \section{Postscript...}
 *  %
 *  This is a section of documentation with its
 *  own \LaTeX markup.
 *  
 */

