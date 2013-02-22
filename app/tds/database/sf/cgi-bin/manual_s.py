#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 31/08/04
#Changed from using frame to tables (html) on: 17/09/04
#Last modification: 15/10/04


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


if not form.has_key("data"):
    # No form data means this is the first access; output the form.
    
    try:
        # Get the first part of the content from the HTML file.
        fp=open("/var/www/html/sf/manual_s1.html", "r")
        listOfLines= fp.readlines()
        for eachline in listOfLines:
            print eachline.strip()
        fp.close()

    except Exception:
        # Show an error msg if the HTML file fails to open.
        print "Sorry, there is an error in opening the html file."

    # Create a blank form using HTML.
    print """<form method="post" action="">
         <table width="100%" border="0" cellspacing="1" cellpadding="1">
         <td width="100%" colspan="2" align="top">
         <table width="100%" border="0" cellspacing="2" cellpadding="2">
         <td width="20%" rowspan="2" valign="top">
         <b>Type your query here: </b></td>
         <td width="80%">
         <input type="text" name="data" value="" size="100" tabindex="1" />
         </td>
         <tr><td width="30%">
         <input type="submit" value="SUBMIT QUERY" tabindex="2">&nbsp;
         <input type="reset" value="CLEAR" tabindex="3">
         </td></tr>
         </table></td>
         
         <tr><td colspan="2"><br></td></tr>
         <tr><td width="15%"><b>
         Example of an query:</b></td>
         <td width="85%"><b>
         SELECT project, blame, shot_number FROM shot_descriptions WHERE project ~ 'scramjet'
         (SINGLE QUOTE)         
         </b></td></tr>
         <tr><td colspan="2"><br /></td></tr>
         <tr><td colspan="2">
         To include the raw data in the output, please tick here:
         <input type=checkbox name="raw_data" value="show_raw_data" />
         </td></tr>
         <tr><td colspan="2">
         To display results in table format
         (preferably for small number of results), please tick here:
         <input type=checkbox name="r_in_table" value="result_in_table" />
         </td></tr>
         <tr><td colspan="2">
         For more help on SQL, please click
         <a href="http://www.mech.uq.edu.au/sf/manual_help.html">here</a>.</td>

         <tr><td width="100%" colspan="2">
         <hr size="3" align="center" noshade="noshade"></td></tr>
                           
         </table>
         </form>"""

    try:
        # Get the second part of the content from the HTML file.
        fp=open("/var/www/html/sf/manual_s2.html", "r")
        listOfLines= fp.readlines()
        for eachline in listOfLines:
            print eachline.strip()
        fp.close()

    except Exception:
        # Show an error msg if the HTML file fails to open.
        print "Sorry, there is an error in opening the html file."
    

else:
    try:
        # Open the HTML page to display the results.
        fp=open("/var/www/html/sf/top_bar.html", "r")
        listOfLines= fp.readlines()
        for eachline in listOfLines:
            print eachline.strip()
        fp.close()

    except Exception:
         # Show an error msg if the HTML file fails to open.
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

    # Grab the values entered by the user from the form.
    my_query = form.getvalue("data") 
    raw_data = form.getvalue("raw_data") or ''
    r_in_table = form.getvalue("r_in_table") or ''

    # Show user's query.
    print my_query

    # Create a horizontal line to separate the text.   
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
    
    
    # Define a function to show raw data.
    def show_raw_data():
        # Show all raw data.
        print "<p><h3>Raw data obtained from database:</h3></p>"
        
        # Show the raw data that are obtained from the database.
        print qresult
        print ""
    
        print "<p><hr size=3 align=center noshade=noshade></p>"
        
        # Show a list of raw data in dictionary form.
        print "<p><h3>List of Results (Raw data after treatment):</h3></p>"
        for record in listOfResults:
            print str(record), "<br>"

        return

    # Define a function to display results in a table format.
    def result_in_table():
        print "<p><h3><u>RESULTS</u></h3></p>"
        print """<body bgcolor=#ffffff text=#000000><font face=Arial, Helvetica,
              sans-serif>"""

        print "<table width=100% border=0 cellspacing=2 cellpadding=2>"
        record = listOfResults[0]
        print "<tr>"
        for k in record.keys():
            print "<td><h3>", k, "</h3></td>"
        print "</tr>"
        for record in listOfResults:
            print "<tr>"
            for k in record.keys():
                print "<td valign=top>", record[k], "</td>"
            print "</tr>"
        print "</table>"
    
        return

    # Define a function to display each record individually.
    def individual_result():
        print "<p><h3><u>RESULTS</u></h3></p>"              
        print """<body bgcolor=#ffffff text=#000000><font face=Arial, Helvetica,
              sans-serif>"""
  
        print "<table width=100% border=0 cellspacing=2 cellpadding=2>"
        print """<td width="20%"><b><u>Field Names</u></b></td>
        <td width="80%"><b><u>Values</u></b></td>"""

        for record in listOfResults:
            print ""
            for k in record.keys():
                print "<tr><td width=20% valign=top>",k,\
                      "</td><td width=80% valign=top>", record[k],"</td></tr>"
        print "</table>"

        return

    if (raw_data != ('')) and (r_in_table != ('')):
        # If both raw_data r_in_table boxes are not equal to none.
        # This means that both boxes are ticked, output the following.

        # Show the raw data obtained from the database.
        show_raw_data()
        print ""
        print "<p><hr size=3 align=left noshade=noshade></p>"
        
        # Display results in a table format(columns=fields and rows=records). 
        result_in_table()
        print ""

        # Create a horizontal line to separate the text and to end the page.
        print "<p><hr size=3 align=center noshade=noshade></p>"
        
        # Create a link for user to get back to the original manual search page.  
        print """<p><a href="http://www.mech.uq.edu.au/cgi-bin/sf/manual_s.py">
              <b>Back</b></a></p>"""        
        
    elif (raw_data != ('')):
        # If the raw_data box is not equal to none,
        # which indicates that the box is ticked, output the following.

        # Show the raw data obtained from the database.
        show_raw_data()
        print ""
        print "<p><hr size=3 align=left noshade=noshade></p>"
        
        # Display each record individually.        
        individual_result()
        print ""
        print "<p><hr size=3 align=left noshade=noshade></p>"

        print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/manual_s.py">
        <b>Back</b></a>"""
        print ""

    elif (r_in_table != ('')):
        # If the r_in_table box is not equal to none,
        # which indicates that the box is ticked, output the following.
        
        # Display results in a table format(columns=fields and rows=records). 
        result_in_table()
        print ""
        print "<p><hr size=3 align=center noshade=noshade></p>"

        print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/manual_s.py">
        <b>Back</b></a>"""
        print ""
    
    else:
        # If both raw_data and r_in_table boxes are not tick,
        # output the following.
        
        # Display each record individually.        
        individual_result()
        print ""
        print "<p><hr size=3 align=left noshade=noshade></p>"

        print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/manual_s.py">
        <b>Back</b></a>"""
        print ""        
  
               
    db.close()     # Close the database connection.


print "</BODY></HTML>"    # End the script.


















































































































































































































































































































































































































































































































































































































    

    
