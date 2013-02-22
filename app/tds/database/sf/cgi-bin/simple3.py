#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 22/08/04
#Modified on: 24/08/04 and 31/08/04
#Finalised on: 30/08/04 (By Dr.PJ)

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
<H2>My First CGI-Script Trial Query Form!</H2>
<H3><I>Please type your query in the text box provided.\n</I></H3><BODY>"""


#Get the form data, if any
form = cgi.FieldStorage()

#No form data means this is the first access; output the form.
if not form.has_key("data"):
    print"""<FORM METHOD="POST" ACTION="/cgi-bin/sf/simple3.py">

<P>
Type your query here:<BR>
<TEXTAREA NAME=data VALUE="" ROWS=3 COLS=60 TYPE=text></TEXTAREA><BR>
<P>
<INPUT TYPE=submit VALUE="Submit"><BR>
</FORM>"""

else:
    #We do have form data; just show it.
    print "You typed this:"
    my_query = form.getvalue("data")
    print my_query
    # It seems that we need to have an appropriate username that matches
    # an entry in the postgresql table of users.
    import os
    username = os.environ.get('USER')
    # print "username: ", username, type(username)
    if username == None:
        # Assume that the web server has started this script and has the
        # username 'apache'.  To allow this to access the database, we had
        # to create a postgresql user of that name and have no password.
        # This new user is also able to create tables and new users,
        # otherwise the access seems to be blocked.  (This might be a
        # security problem but only for our database tables.)
        username = 'apache'

    # Now, we can get to the database...
    import pg
    db = pg.connect("moncdata", user=username, passwd=None)
    qresult = db.query(my_query)
    listOfResults = qresult.dictresult()
    # Have a look at http://www.pygresql.org/README.txt to get documentation
    # on the pgqueryobject methods.
    # Make sure that we have a string.
    resultString = repr(listOfResults)
    print "<P>Raw result obtained from database:</P>"
    print resultString

    print ""
    print "<P>Example of pulling the list of dictionary results apart.</P>"
    for record in listOfResults:
        print "<P><table>"
        for k in record.keys():
            print "<tr> <td>key:</td> <td>", k, "</td> <td>value:</td> <td>", \
                  record[k], "</td> </tr>"
        print "</table></P>"
    db.close()

#HTML end
print"</BODY></HTML>"




        
