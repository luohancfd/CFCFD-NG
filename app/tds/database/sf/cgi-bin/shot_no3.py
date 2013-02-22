#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 30/09/04
#Modified on: 02/10/04
#Help from Lee Harr, Danny Yoo and Steve - Python Tutor

import cgi, os, sys, string, time
import cgitb; cgitb.enable()
import pg

def form():
    print """<form method="post" action="">
             <p>
             <input type=checkbox name="qtype" value="*" checked />
             All<br />
             <input type=checkbox name="qtype" value="project" />
             Project<br />
             <input type=checkbox name="qtype" value="date_string" />
             Date<br />
             <input type=checkbox name="qtype" value="blame" />
             Blame<br />
             <input type=checkbox name="qtype" value="nozzle" />
             Nozzle<br />
             <input type=checkbox name="qtype" value="notes" />
             Notes<br /></p>

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

        for qtype in data.getvalue('qtype'):

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
            print """<p>Example of pulling the list of dictionary results apart.</p>"""
            for record in listOfResults:
                print "<p><table>"
                for k in record.keys():
                    print '<tr>'
                    print '<td>key:</td> <td>', k, '</td>'
                    print '<br />'
                    print '<td>value:</td><td>', record[k], '</td>'
                    print '</tr>'
                print '</table></p>'

            db.close()

        else:
            print "You have no enter a shot number!"
            print "Please type a in shot number in order to perform a search!"
                              

    else:
        form()

print "</body></html>"
    

    


        


