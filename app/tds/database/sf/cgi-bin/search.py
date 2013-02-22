#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 20/09/04

import sys, os

#Import the CGI module.
import cgi

#Turn on CGI debugging info.
import cgitb; cgitb.enable()


#Required header that tells the browser how to render the HTML.
print"Content-Type: text/html\n\n"

#Open out the html form
try:
    fp=open("/var/www/cgi-bin/sf/form.html", "r")
    listOfLines= fp.readlines()
    for eachline in listOfLines:
        print eachline.strip()
    fp.close()
except Exception:
        print "Sorry, there is an error in opening the html file."

if __name__ == "__main__":
    sys.stderr = sys.stdout
    form = cgi.FieldStorage()
    data = SearchForm(form)
    if data.type == "shot_number":
        print "You typed this:"

        my_query = form.getvalue("text")
        print my_query

        import os
        username = os.environ.get('USER')
        # print "username: ", username, type(username)
        if username == None:
            username = 'apache'
        # Now, we can get to the database...
        import pg
        db = pg.connect("moncdata", user=username, passwd=None)
        query = "select * from shot_descriptions where shot_number = '%(text)d'"
        qresult = db.query(query)
        listOfResults = qresult.dictresult()
        print ""
        print "<P>Example of pulling the list of dictionary results apart.</P>"
        for record in listOfResults:
            print "<P><table>"
            for k in record.keys():
                print "<tr><td>key:</td><td>", k, "</td><td>value:</td><td>", \
                      record[k], "</td></tr>"
        print "</table></P>"
        db.close()
    else:
        print "You typed this:"

        my_query = form.getvalue("text")
        print my_query

        import os
        username = os.environ.get('USER')
        # print "username: ", username, type(username)
        if username == None:
            username = 'apache'
        # Now, we can get to the database...
        import pg
        db = pg.connect("moncdata", user=username, passwd=None)
        query = "select * from shot_descriptions where project = '%(text)s'"
        qresult = db.query(query)
        listOfResults = qresult.dictresult()
        print ""
        print "<P>Example of pulling the list of dictionary results apart.</P>"
        for record in listOfResults:
            print "<P><table>"
            for k in record.keys():
                print "<tr><td>key:</td><td>", k, "</td><td>value:</td><td>", \
                      record[k], "</td></tr>"
            print "</table></P>"
        db.close()

#HTML end
print"</BODY></HTML>"




        
