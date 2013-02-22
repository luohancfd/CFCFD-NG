#!/usr/bin/env python
#Author: Chong Soo Fern
#Modify on: 10/09/04

import sys
sys.stderr = sys.stdout

#Import the CGI module.
import cgi

#Turn on CGI debugging info.
import cgitb; cgitb.enable()


#Required header that tells the browser how to render the HTML.
print"Content-Type: text/html\n\n"

#Get the content from the HTML file.
fp=open("/var/www/html/sf/index.html", "r")
listOfLines= fp.readlines()
for eachline in listOfLines:
    print eachline.strip()
fp.close()


#Get the form data, if any
form = cgi.FieldStorage()

#No form data means this is the first access; output the form.

if not form.has_key("data"):
    print"""<FORM METHOD="POST" ACTION="/cgi-bin/sf/trial2.py">
    
<input type=radio name="query" value="shot_number">Shot Number<br>
<input type=radio name="query" value="project">Project<br>
<input type=radio name="query" value="blame">Blame<br>

<p>
Type you search words:<br>
<input typetext name="search_item" size=60>
<br></p>

<input type="submit" value="SEARCH">

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
        username = 'apache'

    # Now, we can get to the database...
    import pg
    db = pg.connect("moncdata", user=username, passwd=None)
    qresult = db.query(my_query)
    listOfResults = qresult.dictresult()
    
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




        
