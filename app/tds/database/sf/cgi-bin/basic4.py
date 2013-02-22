#!/usr/bin/python
#Author: Chong Soo Fern
#Created on: 20/08/04
#Modified on: 23/08/04 and 24/08/04

import sys
sys.stderr = sys.stdin

#Import the regular expression module.
import re

#Import the CGI module.
import cgi
import cgitb; cgitb.enable()

##############################################################################
print"Content-Type: text/html\n\n"

print "hello"

#Specify the filename of the template file.
TemplateFile = "/var/www/html/sf/template.html"

#Specify the filename of the form to show the user.
FormFile = "/var/www/html/sf/form.html"

#Define a new function called Display.
#it takes one parameter - a string to Display.
def Display(Content):
    TemplateHandle = open(TemplateFile, "r") #open in read mode only.
    #read the entire file as a string
    TemplateInput = TemplateHandle.read()
    TemplateHandle.close()                   #close the file.

class BadTemplateException(Exception):
    def __str__(self):
        return "There was a problem with the HTML template."



############################################################################

###What follows are the two main 'action' functions.
###One to show the form and another to process it.

#Define another function for form.
def DisplayForm():
    FormHandle = open(FormFile, "r")
    FormInput = FormHandle.read()
    FormHandle.close()

    Display(FormInput)

#Define a function to process the form.
def ProcessForm(form):
    #extract the information from the form in easily digestible format.
    try:
        name = form["name"].value

    except:
        #name is required. so output an error
        #if not given and exit script
        Display("You need to at least supply a name. Please go back.")
        raise SystemExit

    try:
        email = form["email"].value

    except:
        email = None

    try:
        colour = form["colour"].value

    except:
        colour = None

    try:
        comment = form["comment"].value

    except:
        comment = None

    Output = "" #our output buffer, empty at first.

    Output = Output + "Hello,"

    if email != None:
        Output = Output + "<A HREF=\"mailto:s4046441@student.uq.edu.au\" + email + "">" + name + "</A>.<P>"

    else:
        Output = Output + name + ".<P>"


    if colour == "swallow":
        Output = Output + "You must be a Monty Python fan.<P>"

    elif colour != None:
        Output = Output + "Your favourite colour was "+colour+"<P>"

    else:
        Output = Output + "You cheated! You didn't specify a colour!<P>"


    if comment != None:
        Output = Output + "In addition, you said:<BR>"+comment+"<P>"

    Display(Output)

###############################################################################
    
###Begin actual script.

###Evaluate CGI request
form = cgi.FieldStorage()

###"key" is a hidden form element with an
###action command such as "process".

try:
    key = form["key"].value

except:
    key = None


if key == "process":
    ProcessForm(form)

else:

    DisplayForm()
    


    
