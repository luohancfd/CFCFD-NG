#A signal to the shell that This Is A Python Program.
#!/usr/bin/python
#Author: Chong Soo Fern
#Created on: 16/08/04

#Import the CGI module
import cgi

#Required header that tells the browser how to render the HTML.
print"Content-Type: text/html\n\n"
print""

#Title of the output page.
print"<TITLE>CGI script output</TITLE>"
#Heading 1 in HTML format.
print"<H1>This is my first CGI script</H1>"      

#Print a test string.
print"Hello, world!"


