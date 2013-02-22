#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 30/09/04
#Modified on: 02/10/04
#Help from Lee Harr and Danny Yoo - Python Tutor

import cgi, os, sys, string
import cgitb; cgitb.enable()
import pg

def write_the_blank_form():
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
             <input type=checkbox name="qtype" value="notes" />
             Notes<br /></p>

             <p>Search by shot number:<br>
             <input type=text name="shot_no" value="" tabindex="1" />
             <input type="submit" value="SEARCH" tabindex="2" />&nbsp;
             <input type="reset" value="RESET" tabindex="3" /><br />
             </form>"""
    return

if __name__ == "__main__":
    print "Content-Type: text/html\n\n"
    print "<head><title>Searching by using Shot Number</title></head><body>"
    
    data = cgi.FieldStorage()

    if data:
        shot_no = data.getfirst('shot_no','')
        if len(shot_no) > 0:
            print "You typed in Shot Number:", shot_no
                
            #Now, we can get to the database...
            username = os.environ.get('USER')
            if username == None:
                username = 'apache'
            db = pg.connect("moncdata", user=username, passwd=None)

            query = "select "
            qtype_list = data.getlist('qtype')
            for qtype in qtype_list:
                query += "%s," % (qtype,)
            query = query[:-1];  # throw away the last comma
            query += " from shot_descriptions where shot_number=%s" % (shot_no,)
            qresult = db.query(query)
                    
            listOfResults = qresult.dictresult()

            # print """<p>Example of pulling the list of dictionary results apart.</p>"""
            for record in listOfResults:
                print "<p><table>"
                for k in record.keys():
                    print '<tr>'
                    print '<td>key:</td> <td>', k, '</td>'
                    print '<td>value:</td><td>', record[k], '</td>'
                    print '</tr>'
                print '</table></p>'

            db.close()

        else:
            print "You have no enter a shot number!"
            print "Please type a in shot number in order to perform a search!"

    else:
        write_the_blank_form()


print "</body></html>"
    

