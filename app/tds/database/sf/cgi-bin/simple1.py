#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 22/08/04

import sys
sys.stderr = sys.stdout

#Import the CGI module.
import cgi

#Turn on CGI debugging info.
import cgitb; cgitb.enable()

#Required header that tells the browser how to render the HTML.
print"Content-Type: text/html\n\n"

#HTML head
print"""<HTML><TITLE>Test CGI</TITLE>
<H2>My First CGI-Script Form!</H2><BODY>"""

#Get the form data, if any
form = cgi.FieldStorage()

#No form data means this is the first access; output the form.
if not form.has_key("data"):
    print"""<FORM METHOD="POST" ACTION="/cgi-bin/sf/simple1.py">

Type your data here: <INPUT NAME=data VALUE="" SIZE=60 TYPE=text><BR>
<INPUT TYPE=submit VALUE="Submit"><BR>
</FORM>"""

else:
    #We do have form data; just show it.
    print"You typed:"
    print form.getvalue("data")

#HTML end
print"</BODY></HTML>"




        
