#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 29/09/04
#Help for Lee Harr - Python Tutor

import sys, os
import cgi
import cgitb; cgitb.enable()
import pg

def form():
    print """<form method="post" action="">
             <p><input type=radio name="qtype"value="shot_number" checked />
             Shot Number<br />
             <input type=radio name="qtype" value="date_string" />
             Date</p>

             <p>Type your keywords search here:<br>
             <textarea name="qtext" value="" rows=2 cols=80 type="text">
             </textarea></p>
             <input type="submit" value="SEARCH"><br>
             </form>"""


print "Content-Type: text/html\n\n"
print '<head><title>The Title</title></head><body>'


if __name__ == "__main__":

    data = cgi.FieldStorage()

    if data:
        qtype = data['qtype'].value
        try:
            qtext = data['qtext'].value
        except KeyError:
            qtext = ''

        if qtext and qtype=="shot_number":
            print "You typed this:", qtext

            # Now, we can get to the database...
            username = os.environ.get('USER')
            if username == None:
                username = 'apache'
                
            db = pg.connect("moncdata", user=username, passwd=None)
            query = "select * from shot_descriptions where shot_number=%(qtext)s" % {'qtext': qtext}
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


        elif qtext and qtype == "date":
            print "You typed this:", qtext
            # Now, we can get to the database...
            username = os.environ.get('USER')
            if username == None:
                username = 'apache'
                
            db = pg.connect("moncdata", user=username, passwd=None)
            query = "select * from shot_descriptions where date_string=%(qtext)s" % {'qtext': qtext}
            qresult = db.query(query)           

            listOfResults = qresult.dictresult()
            resultstring = repr(listOfResults)
            print "<p>Raw results obtained from the database:</p>"
            print resultstring

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
            print 'No text entered!'

    else:
        form()

print '<body></html>'
