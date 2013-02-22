#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 30/09/04
#Modified on: 02/10/04
#Help from Lee Harr, Danny Yoo and Steve - Python Tutor

import cgi, os, sys, string
import cgitb; cgitb.enable()
import pg

def form():
    print """<form method="post" action="">
           <table border="0" width="750" cellspacing="2" cellpadding="2">
           <tr><td width="150" align="left">
           <input type=checkbox name="qtype" value="*" checked />All</td>
           <td width="150" align="left">
           <input type=checkbox name="qtype" value="shocktube" />Shocktube</td>
           <td width="150" align="left">
           <input type=checkbox name="qtype" value="diaphragm1" />Diaphragm 1
           </td>
           <td width="150" align="left">
           <input type=checkbox name="qtype" value="driver0" />Driver 0</td>

           <tr><td width="150" align="left">
           <input type=checkbox name="qtype" value="facility_name" />Facility
           Name</td>
           <td width="150" align="left">
           <input type=checkbox name="qtype" value="nozzle" />Nozzle</td>
           <td width="150" align="left">
           <input type=checkbox name="qtype" value="diaphragm2" />Diaphragm 2
           </td>
           <td width="150" align="left">
           <input type=checkbox name="qtype" value="driver1" />Driver 1</td>

           <tr><td width="150" align="left">
           <input type=checkbox name="qtype" value="project" />Project</td>
           <td width="150" align="left">
           <input type=checkbox name="qtype" value="reservoir" />Reservoir</td>
           <td width="150" align="left">
           <input type=checkbox name="qtype" value="diaphragm3" />Diaphragm 3
           </td>
           <td width="150" align="left">
           <input type=checkbox name="qtype" value="driver2" />Driver 2</td>

           <tr><td width="150" align="left">
           <input type=checkbox name="qtype" value="blame" />Blame</td>
           <td width="150" align="left">
           <input type=checkbox name="qtype" value="piston" />Piston</td>
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
           <input type=checkbox name="qtype" value="basic_file_name" />Basic
           File Name</td>
           </table>
         </td>
         
   
         <p>Search by shot number:<br>
         <input type=text name="shot_no" value="" tabindex="1" />
         <input type="submit" value="SEARCH" tabindex="2" />&nbsp;
         <input type="reset" value="RESET" tabindex="3" /><br />
         </form>"""
             
          
if __name__ == "__main__":

    print "Content-Type: text/html\n\n"
    print "<head><title>Quick Search using Run No</title></head><body>"
    
    sys.stderr = sys.stdout
    data = cgi.FieldStorage()

    if data:
        shot_no = data.getfirst('shot_no', '')
        qtype = data.getvalue('qtype')

        if qtype == "*":

            print "You entered Shot Number:", shot_no
            
            #Will take care of the security problems later...
            #Defined a username and connect to the database...
            username = os.environ.get('USER')
            if username == "None":
                username = "apache"

            db = pg.connect("moncdata", user=username, passwd=None)
            query = "SELECT * FROM shot_descriptions WHERE shot_number=%s" % (shot_no)
            qresult = db.query(query)
            
            listOfResults = qresult.dictresult()
            print """<p>Example of pulling the list of dictionary results apart.</p>"""
            for record in listOfResults:
                print "<p><table>"
                for k in record.keys():
                    print '<tr>'
                    print '<td>key:</td> <td>', k, '</td>'
                    print '<br />'
                    print '<td>value:</td><td>', record[k], '</td>'
                    print '</tr>'
                print "</table></p>"

            db.close()

        elif qtype != "*":
            
            print "You entered Shot Number:", shot_no
            
            #Will take care of the security problems later...
            #Defined a username and connect to the database...
            username = os.environ.get('USER')
            if username == "None":
                username = "apache"

            db = pg.connect("moncdata", user=username, passwd=None)
            query = "SELECT %s FROM shot_descriptions WHERE shot_number=%s" % (qtype, shot_no)
            qresult = db.query(query)
            
            listOfResults = qresult.dictresult()
            print listOfResults

            db.close()

        elif not qtype:
            print "You have not select an option!\n"
            print "Please tick at least one checkbox provided!\n"
            

        elif not shot_no:
            print "You have no enter a shot number\n!"
            print "Please type a in shot number in order to perform a search\n!"
                             
    else:
        form()

print "</body></html>"
    

    


        


