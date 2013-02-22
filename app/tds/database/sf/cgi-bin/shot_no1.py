#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 30/09/04
#Help from Lee Harr and Danny Yoo - Python Tutor

import sys, os
import cgi
import cgitb; cgitb.enable()
import pg

def form():
    print """<form method="post" action="">
             <p>
             <input type=checkbox name="qtype" value="all" checked />
             All<br />
             <input type=checkbox name="qtype" value="project" />
             Project<br />
             <input type=checkbox name="qtype" value="date_string" />
             Date<br />
             <input type=checkbox name="qtype" value="blame" />
             Blame<br />
             <input type=checkbox name="qtype" value="notes" />
             Notes<br /></p>

             <p>Search by shot number:<br>
             <input type=text name="qtext" value="" />
             <input type="submit" value="SEARCH"><br>
             </form>"""


print "Content-Type: text/html\n\n"
print "<head><title>Searching by using Shot Number</title></head><body>"

if __name__ == "__main__":
    
    data = cgi.FieldStorage()

    if data:
        qtype = data.getfirst('qtype')
        qtext = data.getfirst('qtext','')

        if qtype == "all":
            print "You typed in Shot Number:", qtext

            #Now, we can get to the database...
            username = os.environ.get('USER')
            if username == None:
                username = 'apache'
                
            db = pg.connect("moncdata", user=username, passwd=None)
            query = "select * from shot_descriptions where shot_number=%(qtext)s"  % {'qtext': qtext}
            qresult = db.query(query)

            listOfResults = qresult.dictresult()

            print """<p>Example of pulling the list of dictionary results apart.</p>"""
            for record in listOfResults:
                print "<p><table>"
                for k in record.keys():
                    print '<tr>'
                    print '<td>key:</td> <td>', k, '</td>'
                    print '<td>value:</td><td>', record[k], '</td>'
                    print '</tr>'
                print '</table></p>'

            db.close()

        elif qtype == "project":
            print "You typed in Shot Number:", qtext

            #Now, we can get to the database...
            username = os.environ.get('USER')
            if username == None:
                username = 'apache'
                
            db = pg.connect("moncdata", user=username, passwd=None)
            query = "select project from shot_descriptions where shot_number=%(qtext)s"  % {'qtext': qtext}
            qresult = db.query(query)

            listOfResults = qresult.dictresult()

            print """<p>Example of pulling the list of dictionary results apart.</p>"""
            for record in listOfResults:
                print "<p><table>"
                for k in record.keys():
                    print '<tr>'
                    print '<td>key:</td> <td>', k, '</td>'
                    print '<td>value:</td><td>', record[k], '</td>'
                    print '</tr>'
                print '</table></p>'

            db.close()

        elif qtype == "date_string":
            print "You typed in Shot Number:", qtext

            #Now, we can get to the database...
            username = os.environ.get('USER')
            if username == None:
                username = 'apache'
                
            db = pg.connect("moncdata", user=username, passwd=None)
            query = "select date_string from shot_descriptions where shot_number=%(qtext)s"  % {'qtext': qtext}
            qresult = db.query(query)

            listOfResults = qresult.dictresult()

            print """<p>Example of pulling the list of dictionary results apart.</p>"""
            for record in listOfResults:
                print "<p><table>"
                for k in record.keys():
                    print '<tr>'
                    print '<td>key:</td> <td>', k, '</td>'
                    print '<td>value:</td><td>', record[k], '</td>'
                    print '</tr>'
                print '</table></p>'

            db.close()

        elif not qtext:
            print "You have no enter a shot number!"
            print "Please type a in shot number in order to perform a search!"

    else:
        form()

print "</body></html>"
            
            
            
            
        


