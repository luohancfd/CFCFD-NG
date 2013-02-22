#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 09/10/04


# Import all necessary modules.
import cgi, os, sys, string

# Turn on CGI debugging info.
import cgitb; cgitb.enable()

# Import the module that support postgreSQL and python interface.
# Note that this is a old module.
import pg

# Import strftime from the time module.
from time import strftime

# To create a string representing time under the control of an explicit format string.
strf_time = strftime('%I:%M %p, %A %d %B, %Y')


# Create a blank form using HTML.
def write_the_blank_form():
    print """<form method="post" action="">
    
    <table border="0" width="750" cellspacing="2" cellpadding="2">
    <tr><td width="150" align="left">
    <input type=checkbox name="qtype" value="*" checked />All</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="shot_number" />Shot Number</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="diaphragm1" />Diaphragm 1
    </td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="driver0" />Driver 0</td>
    
    <tr><td width="150" align="left">
    <input type=checkbox name="qtype" value="facility_name" />Facility
    Name</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="shocktube" />Shocktube</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="diaphragm2" />Diaphragm 2
    </td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="driver1" />Driver 1</td>

    <tr><td width="150" align="left">
    <input type=checkbox name="qtype" value="project" />Project</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="nozzle" />Nozzle</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="diaphragm3" />Diaphragm 3
    </td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="driver2" />Driver 2</td>
    
    <tr><td width="150" align="left">
    <input type=checkbox name="qtype" value="blame" />Blame</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="reservoir" />Reservoir</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="dumptank" />Dumptank</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="jackoff" />Jackoff</td>

    <tr><td width="150" align="left">
    <input type=checkbox name="qtype" value="date_string" />Date</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="condition" />Condition</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="acceltube" />Acceltube</td>
    <td width="150" align="left">
    <input type=checkbox name="qtype" value="piston" />Piston</td>
    </table>
    <br />
    
    <br />
    <table width="100%" border="0" cellspacing="2" cellpadding="2">

    <td width="30%"><b>Please select a search field:</b></td>
    <td width="70%" colspan="3">
    <select name="f_type">
    <option value="project">Project</option>
    <option value="blame">Blame</option>
    <option value="shot_number">Shot Number</option>
    <option value="date_string">Date</option>
    <option value="shocktube">Shocktube</option>
    <option value="nozzle">Nozzle</option>
    <option value="reservoir">Reservoir</option>
    <option value="diaphragm1">Diaphragm 1</option>
    <option value="diaphragm2">Diaphragm 2</option>
    <option value="diaphragm3">Diaphragm 3</option>
    <option value="driver0">Driver 0</option>
    <option value="driver1">Driver 1</option>
    <option value="driver2">Driver 2</option>
    <option value="acceltube">Acceltube</option>
    <option value="condition">Condition</option>
    <option value="dumptank">Dumptank</option>
    <option value="jackoff">Jackoff</option>
    <option value="piston">Piston</option>
    </select></td>

    <tr><td></td></tr>
    <tr><td width="30%"><b>Now, type in your search words:</b></td>
    <td width="37%">
    <input type=text name="s_words" value="" size="50"/></td>
    <td width="8%">
    <input type="submit" value="SEARCH" /></td>
    <td width="25%">
    <input type="reset" value="CLEAR" /></td>

    <tr><td><br /></td></tr>
    <tr><td width="100%" colspan="4">
    To include the raw data in the output, please tick here:
    <input type=checkbox name="raw_data" value="resultString" />
    </td></tr>

    <tr><td><br /></td></tr>
    <tr><td><br /></td></tr>
    <tr><td width="100%" colspan="4">
    <hr width="100%" size="3" align="left" noshade="noshade"></td></tr>
    <tr><td></td></tr>
    <tr><td></td></tr>
    </table>

    <table width="100%" border="0" cellspacing="1" cellpadding="1">
    <td width="40%" align="center"><a href="http://www.python.org">
    <img src="http://www.mech.uq.edu.au/sf/PythonPoweredAnim.gif" alt="Powered by Python" border="0"></a></td>
    <td width="60%" align="center"><a href="http://www.apache.org">
    <img src="http://www.mech.uq.edu.au/sf/apache_pb.gif" alt="Powered by Apache" border="0"></a></td></tr>
    </table>
    
    </form>"""
    
    return


if __name__ == "__main__":

    # Required header that tells the browser how to render the HTML.
    print "Content-Type: text/html\n\n"
    
    # Get the data from the form, if any.            
    data = cgi.FieldStorage()
    
    # So if there are inputs in the form, get these values.
    if data:

        f_type = data.getvalue('f_type')
        s_words = data.getvalue('s_words') or ''
        raw_data = data.getvalue('raw_data')

        # Check if there is an input for shot number.
        if (s_words != ''):
            
            #Open the html page to display the results.
            try:
                fp=open("/var/www/html/sf/top_bar.html", "r")
                listOfLines= fp.readlines()
                for eachline in listOfLines:
                    print eachline.strip()
                fp.close()

            # Show an error msg if the html page fails to open.
            except Exception:
                print ""                
                print "Sorry, there is an error in opening the html file."


            # Set the background colour, font type, title and etc for the result page.
            # This is the result page that displayed the results requested by the users.          
            print """<body bgcolor=#ffffff text=#000000><font face=Arial, Helvetica,sans-serif>"""
            print """<title>Quick Search - Results Page</title>"""
                   
            print """<table width="100%" border="0" cellspacing="1" cellpadding="1">
            <td width="100%" colspan="2"><br /></td>
            <tr><td width="60%" rowspan="2" valign="top"><h2>Results of your search:</h2></td>"""

            # Show the time and date.
            print "<td width=40% align=right>Time: </td>", strf_time

            # Create a link for user to get back to the original quick search page.
            print """<tr><td width="40%" align="right"> 
            <a href="http://www.mech.uq.edu.au/cgi-bin/sf/quick_s.py"> 
            <b>Back</b></a></td></tr> 
            </table>"""

            # Create a horizontal line to separate the previous and the next line.
            print "<p><hr size=3 align=left noshade=noshade></p>"
            print "<b>Your selected field is:</b> ", f_type
            print "<b> and you typed in:</b> ", s_words
            print ""
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

            # Create a query string by combining bits and pieces together.
            query = "select "

            # Grab a list of input from the form.
            qtype_list = data.getlist('qtype')

            # This part gets all the input from the checkboxes.
            # And insert into the query string.
            for qtype in qtype_list:
                query += "%s," % (qtype,)
                
            query = query[:-1];  # throw away the last comma
            query += " from shot_descriptions where %s" % (f_type)
            query += " ~ '%s'" % (s_words)

            # Send the query string to the database and return some results.
            qresult = db.query(query)

            # Format the results into dictionary types of results.
            listOfResults = qresult.dictresult()
            
            # Create a result string.
            resultString = repr(listOfResults)
            

            # If the raw data box is not tick, output the following.          
            if not (data.has_key('raw_data')):

                print "<p><h3>List of Results</h3></p>"
                for record in listOfResults:
                    print ""
                    for k in record.keys():
                        print """<table border="0" width="100%">"""
                        print "<tr><td width=7% valign=top><b>Items:</b></td>\
                        <td width=16% valign=top>", k,\
                        "</td><td width=7% valign=top> \
                        <b>Value:</b> </td><td width=70% valign=top>", \
                        record[k],"</td></tr>"
                    print "</table>"
                    print ""

                # Create a horizontal line to separate the text and to end the page.
                print "<p><hr size=3 align=left noshade=noshade></p>"

                # Create a link for user to get back to the original quick search page.               
                print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/quick_s.py">
                <b>Back</b></a>"""

                # Print a blank line.
                print ""
                

            # If the raw data box is ticked, output the following.
            else:

                print "<p><h3>Raw data obtained from database:</h3></p>"
                print resultString
                print ""
                
                print "<p><hr size=3 align=left noshade=noshade></p>"
                print "<p><h3>List of Results</h3></p>"
                for record in listOfResults:
                    print ""
                    for k in record.keys():
                        print """<table border="0" width="100%">"""
                        print "<tr><td width=7% valign=top><b>Items:</b></td>\
                        <td width=16% valign=top>", k,\
                        "</td><td width=7% valign=top> \
                        <b>Value:</b> </td><td width=70% valign=top>", \
                        record[k],"</td></tr>"
                    print "</table>"
                    print ""

                print "<p><hr size=3 align=left noshade=noshade></p>"
                print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/quick_s.py">
                <b>Back</b></a>"""
                print ""    

            db.close()    # Close the database connection.
            

        # No data means this is the first access so output the form.
        else:
            #Open the HTML page to display the results.
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

                
            # Set the background colour, font type, title and etc for the error page.
            # This is the error page that displayed the error msg.   
            print """<body bgcolor=#ffffff text=#000000><font face=Arial, Helvetica,sans-serif>"""
            print "<title>Quick Search - Error Page</title>"
            
            print """<table width="100%">
            <td width="80%">
            <p><h2>There is an Error occurred in your search!</h2></p></td>
            <td width="20%" align="right">
            <a href="http://www.mech.uq.edu.au/sf/index.html">
            <b>Back</b></a></td>
            </table>"""

            print ""
            print "<p><hr size=3 align=left noshade=noshade></p>"
            print "<p>You have no enter a search words!</p>"
            print "<p>Please type at least one search word in order to perform a search!</p>"
            print """<p>Click <a href="http://www.mech.uq.edu.au/cgi-bin/sf/quick_s.py"><b>here</b></a> to go back to the Quick Search Page.</p>"""

    else:
        #Get the content from the HTML file.
        try:
            fp=open("/var/www/html/sf/quick_search.html", "r")
            listOfLines= fp.readlines()
            for eachline in listOfLines:
                print eachline.strip()
            fp.close()
        
        except Exception:
            print ""
            print "Sorry, there is an error in opening the html file."

        # Show a blank form.
        write_the_blank_form()


print "</body></html>"    # End the script.
    

