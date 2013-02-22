#!/usr/bin/env python
#Created on: 22/09/04

import sys, os

#Import the CGI module.
import cgi

#Turn on CGI debugging info.
import cgitb; cgitb.enable()


#Required header that tells the browser how to render the HTML.
print"Content-Type: text/html\n\n"

#Get the form data, if any
form = cgi.FieldStorage()

#No form data means this is the first access; output the form.

if not form.has_key("type" "data"):
    print"""<form method="post" action="/cgi-bin/sf/try1.py">
    <input type=radio name="type" checked value="shot_number">Shot Number<br>
    <p>
    Type your keywords search here:<br>
    <textarea name="data" value="" rows=2 cols=80 type="text"></textarea><br>
    <p>
    <input type="submit" value="SEARCH"><br>
    </form>"""

else:
    #We do have form data; just show it.
    print "You typed this:"
    
    selection = form.getvalue("type")
    my_query = form.getvalue("data")
    print selection + my_query
    
    import os
    username = os.environ.get('USER')
    # print "username: ", username, type(username)
    if username == None:
        username = 'apache'

    # Now, we can get to the database...
    import pg
    db = pg.connect("moncdata", user=username, passwd=None)
    query = "select * from shot_descriptions where shot_number = my_query"
    qresult = db.query(query)
    listOfResults = qresult.dictresult()
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
