#!/usr/bin/python
#Author: Chong Soo Fern
#Created on: 20/08/04


#Import the CGI module.
import cgi

form = cgi.FieldStorage()

def values(fields):
    for value in fields:
        if value in form: print form[value].value + "<BR>"

if __name__ == "__main__":

    print"Content-Type: text/html\n\n"

    if "submit" and "done" in form:
        values(("name", "age", "email"))

    else:
        print"If you were directed here in error please vist here.com"

else:
    main()





        
