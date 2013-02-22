#!/usr/bin/env python

import cgi, os, sys, string
import cgitb; cgitb.enable()
import PIL, Image


def write_a_blank_form():
    
    print """<form method="post" action="">
    Number:
    <input type="text" name="number" value="" /><br /><br />
    File:
    <input type="file" name="userfile" value="" size="70" /><br /><br/>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
    <input type="submit" value="Add Record" />&nbsp;
    <input type="reset" value="Clear Record" />
    </form>"""

    return


if __name__ == "__main__":

    print "Content-Type: text/html\n\n"
        
    data = cgi.FieldStorage()

    if data:

        number = data.getvalue("number")

        im_file = data["userfile"]
        fp = open("im_file", "rb")
               
        im = Image.open(fp)

        im.show()
        

##         # upload is the name of the directory.
##         fname = '/upload/' + "up.txt"

##         # Create fp file and allow it to write.
##         fp = open(fname, "w")

##         str = "", number

##         fp.write(str)     # Write all the data in str to the file.
##         fp.close()        # Manual close the file.
        
    else:

        write_a_blank_form()

print "</body></html>"
        
