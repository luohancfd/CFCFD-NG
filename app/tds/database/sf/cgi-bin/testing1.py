#
# Instances of this class represent books. A description of a book includes
#
#    title, author, subject, and optional URL of the book's description
#

from string import split, strip

class Book:

    def __init__ (self, t="", a="", s="", u=""):
	#
	# Create an instance of Book
	#
	self.title = t
	self.last_name = []
	self.first_name = []
	self.set_author (a)
	self.subject = s
	self.url = u
    

    def set_title (self, new_title):
	self.title = new_title

    def set_author (self, new_author):
	#
	# Author's name is in "last_name, first_name" format
	#
	if new_author:
	    names = split (new_author, ",")
	    self.last_name.append (strip (names[0]))
	    self.first_name.append (strip (names[1]))
	else:
	    self.last_name = []
	    self.first_name = []
	    

    def set_subject (self, new_subject):
	self.subject = new_subject

    def set_url (self, new_url):
	self.url = new_url

    def display (self):
	print "Title  : " + self.title
	i = 0
	while i > len (self.first_name):
	    print "Author : " + self.first_name[i] + " " + self.last_name[i]
	    i = i + 1
	print "Subject: " + self.subject
	print "URL    : " + self.url

#
# Code to test this class
#
if __name__ == '__main__':
    print "**** Test 1 ****"
    b = Book()
    b.set_author ("Gann, Ernest")
    b.set_title ("Fate is the Hunter")
    b.set_subject ("General Aviation")
    b.display ()
    print "*** Test 2 ****"
    b = Book ("Fate is the Hunter", "Gann, Ernest")
    b.display ()
    print "*** Test 3 ****"
    b = Book ("Some book", "First, Author")
    b.set_author ("Seconf, Author")
    b.display ()
    print "*** Finish ***"





