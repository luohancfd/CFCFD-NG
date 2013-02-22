#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 31/08/04
#Changed from using frame to tables (html) on: 17/09/04
#Last modification on:29/09/04

import sys
sys.stderr = sys.stdout

#Import the CGI module.
import cgi

#Turn on CGI debugging info.
import cgitb; cgitb.enable()


#Required header that tells the browser how to render the HTML.
print"Content-Type: text/html\n\n"

#Get the form data, if any
form = cgi.FieldStorage()

#No form data means this is the first access; output the form.
if not form.has_key("data"):
    #Get the content from the HTML file.
    try:
        fp=open("/var/www/cgi-bin/sf/index.html", "r")
        listOfLines= fp.readlines()
        for eachline in listOfLines:
            print eachline.strip()
        fp.close()
    except Exception:
        print "Sorry, there is an error in opening the html file."
        
    print"""<form method="post" action="/cgi-bin/sf/trial1.py">
         <table width="100%" border="0" cellspacing="0" cellpadding="0">
         <td width="18%" rowspan="7"><br>
         </td>
         <td width="82%" colspan="2">
         Type your query here:<br>
         <textarea name="data" value="" rows="1" cols="90" type="text"
         tabindex="1"></textarea><br>
         <input type="submit" value="Submit Query" tabindex="2">&nbsp;
         <input type="reset" value="Clear" tabindex="3">
         </td>
         <tr><td colspan="2"><br></td></tr>
         <tr><td colspan="2"><br></td></tr>
         <tr>
         <td width="41%" align="center"><a href="http://www.python.org"><img src="/cgi-bin/sf/PythonPoweredAnim.gif" alt="Powered by Python" border="0"></a></td>
         <td width="41%" align="center"><a href="http://www.apache.org"><img src="/cgi-bin/sf/apache_pb.gif" alt="Powered by Apache" border="0"></a></td></tr>
         <tr><td colspan="2"><br></td></tr>
         <tr><td colspan="2"><br></td></tr>
         <tr><td width="82%" colspan="2">
         <hr size="3" align="center" noshade="noshade"></td></tr>
         </form>"""

else:
    #Open the html page to display the results.
    try:
        fp=open("top_bar.html", "r")
        listOfLines= fp.readlines()
        for eachline in listOfLines:
            print eachline.strip()
        fp.close()
    except Exception:
        print "Sorry, there is an error in opening the html file."
        
    print """<body bgcolor=#ffffff text=#000000><font face=Arial, Helvetica,
              sans-serif>"""

    print """<table width="100%">
             <td width="80%"><h2>Results of your search:</h2></td>
             <td width="20%" align="right">
             <a href="http://www.mech.uq.edu.au/cgi-bin/sf/trial1.py">
             <b>Back</b></a></td>
             </table>"""

    print ""
    print "<p><hr size=3 align=left noshade=noshade></p>"
    #We do have form data; just show it.
    print "<b>You typed this:</b>"
    my_query = form.getvalue("data")
    print my_query
    print "<p><hr size=3 align=left noshade=noshade></p>"
    
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
    
    # Make sure that we have a string.
    resultString = repr(listOfResults)
    print "<p><h3>Raw data obtained from database:</h3></p>"
    print resultString
    print ""

    print "<p><hr size=3 align=left noshade=noshade></p>"
    print "<p><h3>List of Results</h3></p>"
    for record in listOfResults:
        print ""
        for k in record.keys():
            print"""<table border="0" width="100%">"""
            print"<tr><td width=7% valign=top><b>Items:</b></td>\
                  <td width=16% valign=top>", k,\
                  "</td><td width=7% valign=top> \
                  <b>Value:</b> </td><td width=70% valign=top>", \
                  record[k],"</td></tr>"
            print"""</table>"""
        print ""

    print "<p><hr size=3 align=left noshade=noshade></p>"
    print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/trial1.py">
             <b>Back</b></a>"""
    print ""    
            
    db.close()

#HTML end
print"</BODY></HTML>"


















































































































































































































































































































































































































































































































































































































    

    
