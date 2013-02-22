#!/usr/bin/env python
# Soo Fern,
#    Success at last!
# This script connects to the database OK and sends its output to
# the web browser.  I have been accessing it via the URL:
# http://127.0.0.1/cgi-bin/sf/test2.py
# and it results in some plain text rendered in my browser window.
#                                          PJ 11-aug-04

# CGI header lines to tell the browser what to expect.
print "Content-type: text/plain"
# print "Content-length: ", len(resultString)
print ""

# We redirect the standard-error stream to the standard-out stream
# so that we can see the error text in the browser.
import sys
sys.stderr = sys.stdout

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
qresult = db.query("select * from shot_descriptions where shot_number = 7399")
listOfResults = qresult.dictresult()
# Have a look at http://www.pygresql.org/README.txt to get documentation
# on the pgqueryobject methods.
# Make sure that we have a string.
resultString = repr(listOfResults)
print "Raw result obtained from database:"
print resultString

print ""
print "Example of pulling the list of dictionary results apart."
for record in listOfResults:
    print "------- start of record -----------"
    for k in record.keys():
        print "key: ", k, "  value:", record[k]
    print "--------- end of record -----------"
db.close()
