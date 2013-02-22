#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 04/10/04
#Modified on: 15/10/04


# Import all necessary modules.
import cgi, os, sys, string

# Turn on CGI debugging info.
import cgitb; cgitb.enable()

# Import strftime from the time module.
# To create a string representing time under the control of an explicit format string.
from time import strftime


# Define a function to output a form.
def write_a_blank_form():
    # Create a blank form using HTML.
    print """<form method="post" action="">
    <table width="100%" border="0" cellspacing="2" cellpadding="2">
    """
    print '<tr><td align="right" colspan="2">%s</td></tr>' % strftime("%I:%M %p, %A %d %B, %Y")    
    print """
    <td width="20%"><b>Facility:</b></td>
    <td>
    <select name="facility_name">
    <option value="T4 Shock Tunnel">T4 Shock Tunnel</option>
    <option value="X1 Expansion Tube">X1 Expansion Tube</option>
    <option value="X2 Expansion Tube">X2 Expansion Tube</option>
    <option value="X3 Expansion Tube">X3 Expansion Tube</option>
    <option value="Drummond Shock Tunnel">Drummond Shock Tunnel</option>
    <option value="Supersonic Blowdown Tunnel">Supersonic Blowdown Tunnel</option>
    </select><br /></td>
    
    <table width="100%" valign="top" border="0" cellspacing="2" cellpadding="2">
    <td width="20%"><b>Updated By:</b></td>
    <td width="80%"><input type="text" name="name" value="" /></td>
    <tr><td width="20%"><b>Shot Number:</b></td>
    <td width="80%"><input type="text" name="shot_number" value="" /></td></tr>
    <tr><td width="20%"><b>Wavelength:</b></td>
    <td width="80%"><input type="text" name="wave_len" value="" />
    nano-metres</td></tr>
    <tr><td width="20%"><b>Image Time:</b></td>
    <td width="80%"><input type="text" name="image_time" value="" />
    micro-second (measured)</td></tr>
    <tr><td width="20%" valign="top"><b>Programmed Delays:</b></td>
    <td width="80%">
    <textarea name="delays" value="" rows="2" cols="90">
    </textarea></td></tr>    

    <tr><td width="20%" valign="top"><b>Remarks:</b></td>
    <td width="80%">
    <textarea name="remarks" value="" rows="5" cols="90">
    </textarea></td></tr>

    <tr><td width="20%" valign="top"><b>Flow Image:</b></td>
    <td width="80%">
    <input type="file" name="image_file" value="" size="80" /></td></tr>
    
    <tr><td width="20%"></td>
    <td width="80%"><input type="submit" value="Add Record" />&nbsp;
    <input type="reset" value="Clear Record" />
    </td></tr>

    <tr><td colspan="2"><br /></td></tr>
    <tr><td colspan="2"><br /></td></tr>

    <tr><td width="100%" colspan="2">
    <hr size="3" noshade="noshade" />
    </td><tr>
    </table>

    <table width="100%" border="0" cellspacing="1" cellpadding="1">
    <td width="40%" align="center"><a href="http://www.python.org">
    <img src="http://www.mech.uq.edu.au/sf/PythonPoweredAnim.gif" alt="Powered by Python" border="0"></a></td>
    <td width="60%" align="center"><a href="http://www.apache.org">
    <img src="http://www.mech.uq.edu.au/sf/apache_pb.gif" alt="Powered by Apache" border="0"></a></td></tr>
    </table>
    </form>"""

    return

if __name__ == "__main__":
    
    # Required header that tells the browser how to render the HTML.
    print "Content-Type: text/html\n\n"

    # Get the data from the form, if any.    
    data = cgi.FieldStorage()

    if data:
        # So if there are inputs in the form, get these values.
        facility_name = data.getfirst('facility_name', '')
        name = data.getvalue('name') or ''
        shot_number = data.getvalue('shot_number') or ''
        wave_len = data.getvalue('wave_len') or ''
        image_time = data.getvalue('image_time')  or ''
        delays = data.getvalue('delays') or ''
        remarks = data.getvalue('remarks') or ''

        
        # Set the background colour, font type, title and etc for the display page.
        # This is the page where users will see after they submit their data.
        # This page show the users what actually they have entered.
        print """<body bgcolor=#ffffff text=#000000><font face=Arial, Helvetica,sans-serif>"""
        print "<title>Flow Image Data Description Text Area</title>"
        print "<h3><u>Flow Image Data Description Text</u></h3><br />"

        # Show the selected option in this field.        
        print "<b>Facility Name:</b> ", facility_name

        
        if name != (''):
            # If there is an input in this field, show it.
            print "<br><b>Updated By:</b> ", name
        
        if shot_number != (''):
            print "<br><b>Shot Number:</b> ", shot_number

        if wave_len != (''):
            print "<br><b>Wavelength:</b> ", wave_len

        if image_time != (''):
            print "<br><b>Total Measured Image Time:</b> ", image_time

        if delays != (''):
            print "<br><b>Programmed Delays:</b> ", delays

        if remarks != (''):
            print "<br><b>Remarks:</b> ", remarks

        
        # Create two break between the previous and the next line.
        print "<br /><br />"

        # Create a horizontal line to separate the text and to end the page.
        print """<hr size="3" noshade="noshade" />"""
        print "<br />"    

        # Create a link for user to get back to the original page.
        print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/flow_image.py">
                 <b>Back</b></a>"""
               
        # Create a filename by combining the required string together.
        filename = "%s" % (facility_name)
        filename = filename[:2];          # Dump away the last two words.
        filename += "_"
        filename += "%s" % (shot_number)
        filename += "_"
        filename += "f_image"
        filename += ".xml"

        # upload is the name of the directory where the file will be saved to.
        fname = 'upload/' + filename

        # Create fp file and allow it to write.
        fp = open(fname, "w")
        
        
        # Create a string that includes all the input data with xml tags.
        str = "<facility_name>" + facility_name + "</facility_name>" + "\n"
        str += "<name>" + name + "</name>" + "\n"
        str += "<shot_number>" + shot_number + "</shot_number>" + "\n"
        str += "<wave_len>" + wave_len + "nano-metres" + "</wave_len>" + "\n"
        str += "<image_time>" + image_time + "micro-second" + "</image_time>" + "\n"
        str += "<delays>" + delays + "</delays>" + "\n"
        str += "<remarks>" + remarks + "</remarks>" + "\n"

       
        fp.write(str)     # Write the string into the file.
        fp.close()        # Manual close the file.

       
    else:
        # No data means this is the first access so output the form.
        try:
            #Get the content from the HTML file.
            fp=open("/var/www/html/sf/f_image.html", "r")
            listOfLines= fp.readlines()
            for eachline in listOfLines:
                print eachline.strip()
            fp.close()
        
        except Exception:
            # Show an error msg if the HTML file fails to open.     
            print ""
            print "Sorry, there is an error in opening the html file."

        # Show a blank form.
        write_a_blank_form()


print "</body></html>"     # End the script.
    
        

        

        

        

    


