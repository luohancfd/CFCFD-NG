#!/usr/bin/env python

#Redirect the standard-error stream to the standard-out stream
#so that we are able to see the error text in the browser.
import sys
sys.stderr = sys.stdout

#Import the CGI module
import cgi

#Required header that tells the browser how to render the HTML.
print"Content-Type: text/html\n\n"

#Define function to generate HTML form.
def generate_form():
    print "hello world"
    print"<HTML>\n"
    print"<HEAD>\n"
    print"\t<TITLE>Info Form</TITLE>\n"
    print"</HEAD>\n"
    print"<BODY BGCOLOR = white>\n"
    print"\t<H3>Please, enter your name and age.</H3>\n"
    print"\t<TABLE BORDER = 0>\n"
    print"\t\t<FORM METHOD = post ACTION = \
    \"basic.py\">\n"
    print"\t\t<TR><TH>Age:</TH><TD><INPUT type = text name = \
    \"age\"</TD></TR>\n>"
    print"\t</TABLE>\n"
    print"\t<INPUT TYPE = hidden NAME = \"action\" VALUE = \
         \"display\">\n"
    print"\t<INPUT TYPE = submit VALUE = \"Enter\">\n"
    print"\t</FORM>\n"
    print"</BODY>\n"
    print"</HTML>\n"

#Define function display data.
#def display_data(name, age):
    print"<HTML>\n"
    print"<HEAD>\n"
    print"\t<TITLE>Info Form</TITLE>\n"
    print"</HEAD>\n"
    print"<BODY BGCOLOR = white>\n"
    print name, "you are", age, "years old."
    print"</BODY>\n"
    print"</HTML>\n"

#Define main function.
def main():
    form = cgi.FieldStorage()
    if (form.has_key("action") and form.has_key("name")\
        and form.has_key("age")):

        if (form["action"].value == "display"):
            display_data(form["name"].value, form['age'].value)

        else:
            generate_form()

#call main function.
main()
          
          
