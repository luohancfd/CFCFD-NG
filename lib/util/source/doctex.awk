# doctex.awk
#
# Purpose:
#     Search for documentation comments within a C source file
#     and put them to stdout.  Output is a LaTeX fragment.  
#
#     There are a couple of types of comment:
#     1. There are general comments of the form /** ... */ as seen
#        in Java programs.  These comments may have embedded HTML tags.
#     2. Comments introducing function prototypes are of the Robodoc form
#        /****f*/
#        //****f*
#        or of the old form used by the original docgen.awk
#        /* @function */
#
# Usage:
#     awk -f doctex.awk file.c > tmp.tex
#     awk -f doctex.awk standalone=1 file.c > tmp.tex
#     awk -f doctex.awk file1.c file2.c > tmp.tex
# Plus:
#     sed 's/&/\\&/g' tmp.tex > mydoc.tex
#
# The sed command is used to workaround a known deficiency with
# awk and the ampersand character
#
# Authors :
#     Peter Jacobs
#     Department of Mechanical Engineering
#     The University of Queensland
#
# Revisions:
#     29-Oct-2007: adapted from docgen.awk
#                  Now I remember why I seldom use regular expressions.
#----------------------------------------------------------


BEGIN {
    print "%%% Don't bother editing this file; it is machine generated.";
}

NR == 1 && standalone == 1 {
    print "\\documentclass[a4paper]{article}";
    print "\\begin{document}";
    cfilename = FILENAME;
    gsub("_", "\\_", cfilename); # LaTeX doesn't like underscores outside math mode
    print "\\section{", cfilename, "}";
}

/\/\*\* /, /\*\// { 
    # Found a Javadoc comment of the form /** ... */
    #
    # It may extend over several lines.
    # Note that a space character following the /** is now
    # required to distinguish the Javadoc comment from the 
    # Robodoc tag line.
    #
    # Print everything on each line but replace /** and */
    # with begin and end paragraph tags.
    #
    sub( /\/\*\*/, "\n\\noindent\n");
    sub( /\*\//, "\n" );
    sub( /^[ |\t]*\*[ |\t]*/, "" ); # eliminate leading * (and spaces)
    gsub( /_/, "\\_"); # LaTeX doesn't like underscores outside math mode
    sub( /[ \t]*/, "" ); # strip off leading whitespace
    if ( length($0) == 0 ) {
        print "\n\\noindent"; # preserve blank lines in LaTeX output
    } else {
	print $0;
    }
} # end of action for general documentation comment */


/(\/\*( |\t)*@function( |\t)*\*\/|\*\*\*\*f)/ { 
    print "% Start of function prototype.";
    print "\\medskip \\noindent";
}


/(\/\*( |\t)*@function( |\t)*\*\/|\*\*\*\*f)/, /)/ { 
    # Found a function prototype which may be spread over several lines.
    # The beginning is marked with a comment of the form /* @function */
    # or of the Robodoc form.
    # The end is the right parenthesis closing the argument list.
    # At the moment, this does not do brilliant job of handling
    # blank spaces or lines between the marker and the function prototype.

    # Remove the opening brace (beginning the function body)
    # if it appears on the same line as other parts of the
    # function prototype.
    sub( /{/, "" );

    # Replace the begin and end tags with LaTeX markup.
    sub( /\/\*( |\t)*@function( |\t)*\*\/[ |\t]*/, "\\texttt{" );
    sub( /\/\/\*\*\*\*f\*[ |\t]*/, "\\texttt{" );
    sub( /\/\*\*\*\*f\*\/[ |\t]*/, "\\texttt{" );
    sub( /)/, ")}" ); # end of line for function prototype 

    sub( /^[ |\t]*\*[ |\t]*/, "" ); # eliminate leading * (and spaces)
    gsub( /_/, "\\_" ); # LaTeX doesn't like underscores outside math mode
    #gsub( /&/, "\\\\\&" ); # nor does it like the & character
                        # (which it considers a delimiter)
    # do the & -> \& replacement in sed 
    # because \& is just poison to awk's gsub function and
    # no matter what we do, we just get the lone & left in

    sub( /^[ \t]*/, "" ); # strip off leading whitespace
    print $0;
} # end of action for a function prototype


END {
    print "% end of extracted documentation comments for", FILENAME, "\n\n";
    if ( standalone == 1 ) {
        print "\\end{document}";
    }
}
