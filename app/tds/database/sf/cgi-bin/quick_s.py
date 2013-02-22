#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 30/09/04
#Last modified on: 15/10/04


# Import all necessary modules.
import cgi, os, sys, string

# Turn on CGI debugging info.
import cgitb; cgitb.enable()

# Import the module that support postgreSQL and python interface.
# Note that this is a old module.
import pg


# Define a function to display a form.
def write_the_blank_form():
    # Create a blank form using HTML.
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
    <input type=checkbox name="raw_data" value="show_raw_data" />
    </td></tr>
    <tr><td width="100%" colspan="4">
    To display results in a table format
    (preferably for small number of results), please tick here:
    <input type=checkbox name="r_in_table" value="result_in_table" />
    </td></tr>    

    <tr><td width="100%" colspan="4">
    <hr width="100%" size="3" align="left" noshade="noshade"></td></tr>
    </table>

    </form>"""
    
    return


if __name__ == "__main__":

    # Required header that tells the browser how to render the HTML.
    print "Content-Type: text/html\n\n"
    
    # Get the data from the form, if any.            
    data = cgi.FieldStorage()
    
    
    if data:
        # So if there are inputs in the form, get these values.
        f_type = data.getvalue('f_type')
        s_words = data.getvalue('s_words') or ''
        raw_data = data.getvalue('raw_data') or ''
        r_in_table = data.getvalue('r_in_table') or ''

        # Check if there is an input for search words.
        if (s_words != ''):
                        
            try:
                #Open the html page to display the results.
                fp=open("/var/www/html/sf/top_bar.html", "r")
                listOfLines= fp.readlines()
                for eachline in listOfLines:
                    print eachline.strip()
                fp.close()

            except Exception:
                # Show an error msg if the html page fails to open.
                print ""                
                print "Sorry, there is an error in opening the html file."


            # Set the background colour, font type, title and etc for the result page.
            # This is the result page that displayed the results requested by the users.          
            print """<body bgcolor=#ffffff text=#000000><font face=Arial, Helvetica,sans-serif>"""
            print """<title>Quick Search - Results Page</title>"""
                   
            print """<table width="100%" border="0" cellspacing="1" cellpadding="1">
            <td width="100%" colspan="2"><br /></td>
            <tr><td width="60%" rowspan="2" valign="top"><h2>Results of your search:</h2></td>"""

            # Create a link for user to get back to the original quick search page.
            print """<tr><td width="40%" align="right"> 
            <a href="http://www.mech.uq.edu.au/cgi-bin/sf/quick_s.py"> 
            <b>Back</b></a></td></tr> 
            </table>"""

            # Create a horizontal line to separate the previous and the next line.
            print "<p><hr size=3 align=left noshade=noshade></p>"

            # Show the users what they have requested in a particular search.
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

            # Grab the type of fields selected by the user
            # and include into the query string.
            query += " from shot_descriptions where %s" % (f_type)

            # Grab the search words typed in by the user
            # and include into the query string.
            query += " ~ '%s'" % (s_words)

            # Send the query string to the database and return some results.
            qresult = db.query(query)

            # Format the results into dictionary types of results.
            listOfResults = qresult.dictresult()
            
            # Create a result string.
            resultString = repr(listOfResults)

                       
            # Define a function to output the raw data.
            def show_raw_data():
                # Show all raw data.
                print "<p><h3>Raw data obtained from database:</h3></p>"
        
                # Show the raw data that are obtained from the database.
                print qresult
                print ""
    
                print "<p><hr size=3 align=center noshade=noshade></p>"
        
                # Show a list of results in dictionary form.
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
                print "<p><hr size=3 align=left noshade=noshade></p>"

                # Create a link for user to get back to the original quick search page.               
                print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/quick_s.py">
                <b>Back</b></a>"""

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
                
                print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/quick_s.py">
                <b>Back</b></a>"""
                print ""

            elif (r_in_table != ('')):
                # If the r_in_table box is not equal to none,
                # which indicates that the box is ticked, output the following.
                
                # Display results in a table format(columns=fields and rows=records). 
                result_in_table()
                print ""
                print "<p><hr size=3 align=center noshade=noshade></p>"
                
                print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/quick_s.py">
                <b>Back</b></a>"""
                print ""
    
            else:
                # If both raw_data and r_in_table boxes are not tick,
                # output the following.
                
                # Display each record individually.        
                individual_result()
                print ""
                print "<p><hr size=3 align=left noshade=noshade></p>"
                
                print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/quick_s.py">
                <b>Back</b></a>"""
                print ""        

            db.close()    # Close the database connection.
            

        else:
            # There is no input for the field - search words.
            try:
                #Open the HTML page to display the results.
                fp=open("/var/www/html/sf/top_bar.html", "r")
                listOfLines= fp.readlines()
                for eachline in listOfLines:
                    print eachline.strip()
                fp.close()

            except Exception:
                # Show an error msg if the HTML file fails to open.
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
            print """<p>Click <a href="http://www.mech.uq.edu.au/cgi-bin/sf/quick_s.py">
            <b>here</b></a> to go back to the Quick Search Page.</p>"""



    else:
        # No data means this is the first access so output the form.
        try:
            #Get the first part of the content from the HTML file.
            fp=open("/var/www/html/sf/quick_s1.html", "r")
            listOfLines= fp.readlines()
            for eachline in listOfLines:
                print eachline.strip()
            fp.close()
        
        except Exception:
            # Show an error msg if the HTML file fails to open.
            print ""
            print "Sorry, there is an error in opening the html file."

        # Show a blank form.
        write_the_blank_form()

        try:
             #Get the second part of the content from the HTML file.
            fp=open("/var/www/html/sf/quick_s2.html", "r")
            listOfLines= fp.readlines()
            for eachline in listOfLines:
                print eachline.strip()
            fp.close()
        
        except Exception:
            # Show an error msg if the HTML file fails to open.
            print ""
            print "Sorry, there is an error in opening the html file."

            db.close()           # Close database connection.
            
print "</body></html>"    # End the script.
    

