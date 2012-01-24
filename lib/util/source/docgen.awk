# docgen.awk
#
# Purpose:
#     Search for documentation comments within a C source file
#     and put them to stdout.  Output is HTML.  
#
#     There are a couple of types of comment:
#     1. There are general comments of the form /** ... */ as seen
#        in Java programs.  These comments may have embedded HTML tags.
#     2. Comments introducing function prototypes are of the form
#        /* @function */
#
# Usage:
#     awk -f docgen.awk file.c > file.html
#
# Author:
#     Peter Jacobs
#     Department of Mechanical Engineering
#     The University of Queensland
#
# Revisions:
#     02-Oct-1999: general Javadoc style comments
#     30-Oct-1999: function prototypes handled specially
#----------------------------------------------------------


BEGIN {
    print "<HTML>"
    print "<!-- Don't bother editing this file; it is machine generated. -->"
    print "<HEAD>"
    print "<TITLE> Documentation extracted from ", ARGV[1], "</TITLE>"
    print "</HEAD>"
    print "<BODY>"
    print "<H2>", ARGV[1], "</H2>"
}


/\/\*\*/, /\*\// { 
    # Found a Javadoc comment of the form /** ... */
    # It may extend over several lines.

    # Print everything on each line but replace /** and */
    # with begin and end paragraph tags.

    something_of_interest_printed = 0
    for(i = 1; i <= NF; ++i) {
        if ($i == "/**") {
            printf( "<P>" )
        } else if ($i == "*/") {
            printf( "</P>\n" )
        } else {
            printf( "%s ", $i )
            something_of_interest_printed = 1
        } # end if
    } # end for i

    if ( something_of_interest_printed == 1 ) {
        printf("\n")
    } # end if
} # end of action for general documentation comment */


/\/\*( |\t)*@function( |\t)*\*\// { 
    # Start of a function prototype
    print "<!-- --------------------------------------------- --> <HR>"
}


/\/\*( |\t)*@function( |\t)*\*\//, /)/ { 
    # Found a function prototype which may be spread over several lines.
    # The beginning is marked with a comment of the form /* @function */.
    # The end is the right parenthesis closing the argument list.

    # Replace the begin and end tags with HTML tags.
    sub( /\/\*( |\t)*@function( |\t)*\*\//, "<P> <CODE> <B>" )
    sub( /)/, ") </B> </CODE> </P>" )

    # Remove the opening brace (beginning the function body)
    # if it appears on the same line as other parts of the
    # function prototype.
    sub( /{/, "" )

    # Print everything remaining items on the line.
    something_of_interest_printed = 0
    for(i = 1; i <= NF; ++i) {
        printf( "%s ", $i )
        something_of_interest_printed = 1
    } # end for i
    if ( something_of_interest_printed == 1 ) {
        printf("\n")
    } # end if
} # end of action for a function prototype


END {
    print "<HR>"
    print "<ADDRESS>"
    print "    Extracted by <B>docgen.awk</B> ", strftime()
    print "</ADDRESS>"
    print "</BODY>"
    print "</HTML>"
}
