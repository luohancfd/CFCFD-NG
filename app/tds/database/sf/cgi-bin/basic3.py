#!/usr/bin/python
#Author: Chong Soo Fern
#Created on: 20/08/04

#Import the CGI module.
import cgi

#Import the regular expression module.
import re

#Specify the filename of the template file.
TemplateFile = "template.html"

#Define a new function called Display.
#it takes one parameter - a string to Display.
def Display(Content):
    TemplateHandle = open(TemplateFile, "r") #open in read mode only.
    #read the entire file as a string
    TemplateInput = TemplateHandle.read()
    TemplateHandle.close()                   #close the file.

    #This defines an exception string in case
    #our template file is messed up.
    BadTemplateException = "There was a problem with the HTML template."

    SubResult = re.subn("<!--***INSERT CONTENT HERE***--",
                        Content, TemplateInput)
    if SubResult[1] == 0:
        raise BadTemplateException

    #Required header that tells the browser how to render the HTML.
    print"Content-Type: text/html\n\n"

    print SubResult[0]

    
