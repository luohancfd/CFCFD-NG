#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 31/08/04
#Changed from using frame to tables (html) on: 17/09/04
#Last modification: 09/10/04


# Import all necessary modules.
import cgi, os, sys

# Turn on CGI debugging info.
import cgitb; cgitb.enable()

# Import the module that support postgreSQL and python interface.
# Note that this is a old module.
import pg


# Required header that tells the browser how to render the HTML.
print"Content-Type: text/html\n\n"

# Get the form data, if any
form = cgi.FieldStorage()


# No form data means this is the first access; output the form.
if not form.has_key("data"):
    
    # Get the content from the HTML file.
    try:
        fp=open("/var/www/html/sf/manual_search.html", "r")
        listOfLines= fp.readlines()
        for eachline in listOfLines:
            print eachline.strip()
        fp.close()

    # Show an error msg if the HTML file fails to open.
    except Exception:
        print "Sorry, there is an error in opening the html file."

    # Create a blank form using HTML.
    print"""<form method="post" action="">
         <table width="100%" border="0" cellspacing="1" cellpadding="1">
         <td width="100%" colspan="2" align="top">
         <table width="100%" border="0" cellspacing="2" cellpadding="2">
         <td width="20%" rowspan="2" valign="top">
         <b>Type your query here: </b></td>
         <td width="80%">
         <textarea name="data" value="" rows="1" cols="90" type="text"
         tabindex="1"></textarea></td>
         <tr><td width="30%">
         <input type="submit" value="SUBMIT QUERY" tabindex="2">&nbsp;
         <input type="reset" value="CLEAR" tabindex="3">
         </td></tr>
         </table></td>
         <tr><td colspan="2"><br></td></tr>
         <tr><td colspan="2"><br></td></tr>
         <tr><td colspan="2">
         To include the raw data in the output, please tick here:
         <input type=checkbox name="raw_data" value="resultString" />
         </td></tr>
         <tr><td width="100%" colspan="2">
         <hr size="3" align="center" noshade="noshade"></td></tr>
         <tr><td colspan="2"><br></td></tr>
         <tr>
         <td width="40%" align="center"><a href="http://www.python.org">
         <img src="http://www.mech.uq.edu.au/sf/PythonPoweredAnim.gif" alt="Powered by Python" border="0"></a></td>
         <td width="60%" align="center"><a href="http://www.apache.org">
         <img src="http://www.mech.uq.edu.au/sf/apache_pb.gif" alt="Powered by Apache" border="0"></a></td></tr>
         
         </table>
         </form>"""

else:
    # Open the HTML page to display the results.
    try:
        fp=open("/var/www/html/sf/top_bar.html", "r")
        listOfLines= fp.readlines()
        for eachline in listOfLines:
            print eachline.strip()
        fp.close()

    # Show an error msg if the HTML file fails to open.
    except Exception:
        print ""
        print "Sorry, there is an error in opening the html file."


    # Set the background colour, font type, title and etc for the result page.
    # This is the result page that displayed the results requested by the users.       
    print """<body bgcolor=#ffffff text=#000000><font face=Arial, Helvetica,
              sans-serif>"""
    print "<title>Manual Search - Results Page</title>"

    print """<table width="100%">
             <td width="80%"><p><h2>Results of your search:</h2></p></td>
             <td width="20%" align="right">
             <a href="http://www.mech.uq.edu.au/cgi-bin/sf/manual_s.py">
             <b>Back</b></a></td>
             </table>"""

    # Print a empty line.
    print ""

    # Create a horizontal line to separate the text.
    print "<p><hr size=3 align=left noshade=noshade></p>"
    
    # Show the users what data they entered in.
    print "<b>You typed this:</b>"
    
    my_query = form.getvalue("data")
    raw_data = form.getvalue('raw_data')
    
    print my_query
    print "<p><hr size=3 align=left noshade=noshade></p>"

    
    # Now, we can get to the database...
    # Create a username in order to connect to the database.
    # It seems that we need to have an appropriate username that matches
    # an entry in the postgresql table of users.
    username = os.environ.get('USER')

    # Assume that the web server has started this script and has the
    # username 'apache'.  To allow this to access the database, we had
    # to create a postgresql user of that name and have no password.
    # This new user is also able to create tables and new users,
    # otherwise the access seems to be blocked.  (This might be a
    # security problem but only for our database tables.)
    if username == None:
        username = 'apache'

    # Connect to the Moncdata database.
    db = pg.connect("moncdata", user=username, passwd=None)

    # Send the query string to the database and return some results.
    qresult = db.query(my_query)

    # Format the results into dictionary types of results.
    listOfResults = qresult.dictresult()

    # Create a result string.
    resultString = repr(listOfResults)
    
    
    # If the raw data box is not tick, output the following.  
    if not (form.has_key('raw_data')):
        
        print "<p><u><h3>List of Results</h3></u></p>"
        print ""

        for record in listOfResults:
            print ""
            for k in record.keys():
                print """<table border="0" width="100%">"""
                print "<tr><td width=7% valign=top><b>Descriptions:</b></td>\
                <td width=16% valign=top>", k,\
                "</td><td width=7% valign=top> \
                <b>Values:</b> </td><td width=70% valign=top>", \
                record[k],"</td></tr>"
                print """</table>"""
            print ""

            
        print "<p><hr size=3 align=left noshade=noshade></p>"

        print ""
        print "<u><h3>Testing Zone 1</h3></u>"
        print """<table border="0" width="100%">"""
        print """<td width="20%"><b>Descriptions:</b></td>
        <td width="80%"><b>Values:</b></td></table>"""
        for record in listOfResults:
            print ""
            
            for k in record.keys():
                print """<table border="0" width="100%">"""
                print "<tr><td width=20% valign=top>",k,\
                      "</td><td width=80% valign=top>", record[k],"</td></tr>"
                
                print """</table>"""
            print ""

        print "<p><hr size=3 align=left noshade=noshade></p>"

        print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/manual_s.py">
        <b>Back</b></a>"""
        print ""

    
    # If the raw data box is ticked, output the following.
    else:

        print "<p><h3>Raw data obtained from database:</h3></p>"
        print qresult
        print ""
    
        print "<p><hr size=3 align=left noshade=noshade></p>"
        
        print "<p><h3>Raw data after some treatment :</h3></p>"
        print resultString
        print ""
    
        print "<p><hr size=3 align=left noshade=noshade></p>"
        
        print "<p><h3>List of Results</h3></p>"

        print "<p><hr size=3 align=left noshade=noshade></p>"
        
        for record in listOfResults:
            print ""
            print """<table border="0" width="100%">"""
            print """<td width="20%"><b>Descriptions:</b></td>
                     <td width="80%"><b>Values:</b></td></table>"""
            
            for k in record.keys():
                print """<table border="0" width="100%">"""
                print "<tr><td width=20% valign=top>",k,\
                 "</td><td width=80% valign=top>", record[k],"</td></tr>"
                
            print """</table>"""

                
    db.close()     # Close the database connection.


print "</BODY></HTML>"    # End the script.


















































































































































































































































































































































































































































































































































































































    

    
