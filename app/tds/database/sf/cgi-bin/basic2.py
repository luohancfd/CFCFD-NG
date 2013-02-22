#!/usr/bin/python
#Author: Chong Soo Fern
#Created on: 20/08/04

#Import the CGI module
import cgi

#Required header that tells the browser how to render the HTML.
print"Content-Type: text/plain\n\n"

  
form = cgi.FieldStorage()

for name in form.keys():
    print"Input:" + name + "value:" + form[name].value + "<BR>"

print"Done!!!"
