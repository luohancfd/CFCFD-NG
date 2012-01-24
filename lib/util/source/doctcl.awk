# doctcl.awk
#
# Purpose:
#     Search for documentation comments within a Tcl source file
#     and put them to stdout.  Output is assumed to be HTML.  
#
#     There are a couple of types of comment:
#     1. There are general comments of the form 
#        #@doc
#        #... 
#        #@end
#        analogous to JavaDoc comments.
#        These comments may have embedded HTML tags.
#     2. Comments terminating procedure declarations
#        #@proc
#        This script will copy the first lines of the
#        procedure declaration, up until this comment
#        is encountered.
#
# Usage:
#     awk -f doctcl.awk file.tcl > file.html
#
# Author:
#     Peter Jacobs
#     Department of Mechanical Engineering
#     The University of Queensland
#
# Revisions:
#     25-Jan-2000: Adapted from docgen.awk
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


/\#\@doc/, /\#\@end/ { 
    # Found a Javadoc type comment.
    # It may extend over several lines.

    # Print everything except the comment characters #
    # and the continuation characters \
    # on each line but replace #@doc and #@end
    # with begin and end paragraph tags.

    something_of_interest_printed = 0
    for(i = 1; i <= NF; ++i) {
        if ($i == "#@doc") {
            printf( "<P>" )
        } else if ($i == "#@end") {
            printf( "</P>\n" )
        } else if ($i == "#") {
            # do nothing
        } else if ($i == "\\") {
            # do nothing
        } else {
            printf( "%s ", $i )
            something_of_interest_printed = 1
        } # end if
    } # end for i

    if ( something_of_interest_printed == 1 ) {
        printf("\n")
    } # end if
} # end of action for general documentation comment */


/^( |\t)*proc( |\t)*/ { 
    # Start of a function prototype
    print "<!-- --------------------------------------------- --> <HR>"
}


/^( |\t)*proc( |\t)*/, /\#@proc/ { 
    # Found a function prototype which may be spread over several lines.
    # The beginning is the keyword proc at the start of the line.
    # The end is marked with a comment of the form #@proc

    # Replace the begin and end tags with HTML tags.
    sub( /^( |\t)*proc( |\t)*/, "<P> <CODE> <B> proc " )
    sub( /\#@proc/, " </B> </CODE> </P>" )

    # Print everything remaining items on the line.
    something_of_interest_printed = 0
    for(i = 1; i <= NF; ++i) {
        printf( "%s ", $i )
        something_of_interest_printed = 1
    } # end for i
    if ( something_of_interest_printed == 1 ) {
        printf("\n")
    } # end if
} # end of action for a procedure header


END {
    print "<HR>"
    print "<ADDRESS>"
    print "    Extracted by <B>doctcl.awk</B> ", strftime()
    print "</ADDRESS>"
    print "</BODY>"
    print "</HTML>"
}
