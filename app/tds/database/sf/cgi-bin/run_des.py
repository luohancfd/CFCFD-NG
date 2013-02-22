#!/usr/bin/env python
#Author: Chong Soo Fern
#Created on: 05/10/04
#Modified on: 15/10/04


# Import all necessary modules.
import cgi, os, sys, string

# Turn on CGI debugging info.
import cgitb; cgitb.enable()    

# Import strftime from the time module.
# To create a string representing time under the control of an explicit format string.
from time import strftime


# Define a function to display a form.
def write_a_blank_form():
    # Create a blank form using HTML.
    print """<form method="post" action="">
    <table width="100%" border="0" cellspacing="2" cellpadding="2">
    """
    print '<tr><td align="right" colspan="2">%s</td></tr>' % strftime("%I:%M %p, %A %d %B, %Y")
    print """
    <tr><td></td><tr>
    <tr><td width="15%"><b>Facility:</b></td>
    <td>
    <select name="facility_name">
    <option value="T4 Shock Tunnel">T4 Shock Tunnel</option>
    <option value="X1 Expansion Tube">X1 Expansion Tube</option>
    <option value="X2 Expansion Tube">X2 Expansion Tube</option>
    <option value="X3 Expansion Tube">X3 Expansion Tube</option>
    <option value="Drummond Shock Tunnel">Drummond Shock Tunnel</option>
    <option value="Supersonic Blowdown Tunnel">Supersonic Blowdown Tunnel</option>
    </select></td></tr></table>
    
    <br />
    <table width="100%" valign="top" border="0" cellspacing="2" cellpadding="2">
    <td width="15%"><b>Shot Number:</b></td>
    <td width="85%"><input type="text" name="shot_number" value="" /></td>
    <tr><td width="15%"><b>Project Name:</b></td>
    <td width="85%"><input type="text" name="project" value="" size="90" /></td></tr>
    <tr><td width="15%"><b>Blame:</b></td>
    <td width="85%"><input type="text" name="blame" value="" size="50" /></td></tr>
    <tr><td width="15%"><b>Date:</b></td>
    <td width="85%">
    <table width="25%" valign="top" border="0" cellspacing="1" cellpadding="1">
    <td width="5%">
    <select name="day">
    <option value="1">1</option><option value="2">2</option><option value="3">3</option>
    <option value="4">4</option><option value="5">5</option><option value="6">6</option>
    <option value="7">7</option><option value="8">8</option><option value="9">9</option>
    <option value="10">10</option><option value="11">11</option><option value="12">12</option>
    <option value="13">13</option><option value="14">14</option><option value="15">15</option>
    <option value="16">16</option><option value="17">17</option><option value="18">18</option>
    <option value="19">19</option><option value="20">20</option><option value="21">21</option>
    <option value="22">22</option><option value="23">23</option><option value="24">24</option>
    <option value="25">25</option><option value="26">26</option><option value="27">27</option>
    <option value="28">28</option><option value="29">29</option><option value="30">30</option>
    <option value="31">31</option>
    </select></td>

    <td width="12%">
    <select name="month">
    <option value="Jan">January</option><option value="Feb">February</option>
    <option value="Mar">March</option><option value="Apr">April</option>
    <option value="May">May</option><option value="Jun">June</option>
    <option value="Jul">July</option><option value="Aug">August</option>
    <option value="Sept">September</option><option value="Oct">October</option>
    <option value="Nov">November</option><option value="Dec">December</option>
    </select></td>

    <td width="8%">
    <select name="year">
    <option value="1985">1985</option><option value="1986">1986</option>
    <option value="1987">1987</option><option value="1988">1988</option>
    <option value="1989">1989</option><option value="1990">2005</option>
    <option value="1991">1991</option><option value="1992">1992</option>
    <option value="1993">1993</option><option value="1994">1994</option>
    <option value="1995">1995</option><option value="1996">1996</option>
    <option value="1997">1997</option><option value="1998">1998</option>
    <option value="1999">1999</option><option value="2000">2000</option>
    <option value="2001">2001</option><option value="2002">2002</option>
    <option value="2003">2003</option><option value="2004">2004</option>
    <option value="2005">2005</option><option value="2006">2006</option>
    <option value="2007">2007</option><option value="2008">2008</option>
    <option value="2009">2009</option><option value="2010">2010</option>
    <option value="2011">2011</option><option value="2012">2012</option>
    <option value="2013">2013</option><option value="2014">2014</option>
    <option value="2015">2015</option><option value="2016">2004</option>
    <option value="2016">2016</option><option value="2017">2017</option>
    <option value="2018">2018</option><option value="2019">2019</option>
    <option value="2020">2020</option><option value="2021">2021</option>
    <option value="2022">2022</option><option value="2023">2023</option>
    <option value="2024">2024</option><option value="2025">2025</option>
 
    </select></td>

    </table></td>
    <tr><td></td><tr></table>
    
    <table width="100%" valign="top" border="0" cellspacing="3" cellpadding="3">
    <td width="15%"><b>Nozzle:</b></td>
    <td width="35%"><input type="text" name="nozzle" value="" size="40" /></td>
    <td width="15%"><b>Shocktube:</b></td>
    <td width="35%"><input type="text" name="shocktube" value="" size="40" /></td>

    <tr><td width="15%"><b>Reservoir:</b></td>
    <td width="35%"><input type="text" name="reservoir" value="" size="40" /></td>
    <td width="15%"><b>Diaphragm 1:</b></td>
    <td width="35%"><input type="text" name="diaphragm1" value="" size="40" />
    </td></tr>
    <tr><td></td></tr>
    </table></td></tr>

    <table width="100%" valign="top" border="0" cellspacing="3" cellpadding="3">
    <td width="15%"><b>Driver 0:</b></td>
    <td width="85%"><input type="text" name="driver0" value="" size="60" /></td>
    <tr><td></td></tr>
    </table>

    <table width="100%" valign="top" border="0" cellspacing="3" cellpadding="3">
    <td width="15%"><b>Driver 1:</b></td>
    <td width="35%"><input type="text" name="driver1" value="" size="40" /></td>
    <td width="15%"><b>Diaphgram 2:</b></td>
    <td width="35%"><input type="text" name="diaphgram2" value="" size="40" /></td>

    <tr><td width="15%"><b>Driver 2:</b></td>
    <td width="35%"><input type="text" name="driver2" value="" size="40" /></td>
    <td width="15%"><b>Diaphragm 3:</b></td>
    <td width="35%"><input type="text" name="diaphragm3" value="" size="40" />
    </td></tr>
    <tr><td></td></tr>
    </table></td></tr>

    <table width="100%" valign="top" border="0" cellspacing="3" cellpadding="3">
    <td width="15%"><b>Acceltube:</b></td>
    <td width="35%"><input type="text" name="acceltube" value="" size="40" /></td>
    <td width="15%"><b>Dumptank:</b></td>
    <td width="35%"><input type="text" name="dumptank" value="" size="40" /></td>

    <tr><td width="15%"><b>Jackoff:</b></td>
    <td width="35%"><input type="text" name="jackoff" value="" size="40" /></td>
    <td width="15%"><b>Piston:</b></td>
    <td width="35%"><input type="text" name="piston" value="" size="40" />
    </td></tr>
    <tr><td></td></tr>
    </table></td></tr>

    <table width="100%" valign="top" border="0" cellspacing="3" cellpadding="3">
    <td width="15%"><b>Condition:</b></td>
    <td width="85%"><input type="text" name="condition" value="" size="60" /></td>
    <tr><td></td></tr>
    </table>    
    
    <table width="100%" valign="top" border="0" cellspacing="3" cellpadding="3">    
    <tr><td width="15%" valign="top"><b>Notes:</b></td>
    <td width="85%">
    <textarea name="notes" value="" rows="3" cols="90">
    </textarea></td></tr>
    
    <tr><td width="15%"></td>
    <td width="85%"><input type="submit" value="Add Record" accesskey="a" />&nbsp;
    <input type="reset" value="Clear Record" accesskey="c" />
    </td></tr>

    <tr><td colspan="2"><br /></td></tr>
    <tr><td colspan="2"><br /></td></tr>

    <tr><td width="100%" colspan="2">
    <hr size="3" noshade="noshade" />
    </td><tr></table>

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
        shot_number = data.getvalue('shot_number') or ''
        project = data.getvalue('project') or ''
        blame = data.getvalue('blame') or ''
        day = data.getvalue('day') or ''
        month = data.getvalue('month') or ''
        year = data.getvalue('year') or ''
        nozzle = data.getvalue('nozzle') or ''
        shocktube = data.getvalue('shocktube') or '' 
        reservoir = data.getvalue('reservoir') or ''
        diaphragm1 = data.getvalue('diaphragm1') or ''
        driver0 = data.getvalue('driver0') or ''
        driver1 = data.getvalue('driver1') or ''
        driver2 = data.getvalue('driver2') or ''
        diaphragm2 = data.getvalue('diaphragm2') or '' 
        diaphragm3 = data.getvalue('diaphragm3') or ''
        acceltube = data.getvalue('acceltube') or ''
        dumptank = data.getvalue('dumptank') or ''
        jackoff = data.getvalue('jackoff') or ''
        piston = data.getvalue('piston') or ''
        condition = data.getvalue('condition') or ''
        notes = data.getvalue('notes')  or ''

        # To combine the day, month and year into a string.
        date_string = "%s " % (day)
        date_string += "%s " % (month)
        date_string += "%s" % (year)

        
        # Set the background colour, font type, title and etc for the display page.
        # This is the page where users will see after they submit their data.
        # This page show the users what actually they have entered.
        print """<body bgcolor=#ffffff text=#000000><font face=Arial, Helvetica,sans-serif>"""
        print "<title>Run Description Text Area</title>"
        
        print "<p><h3><u>Run Description Data</u></h3></p>"
    

        # Show the selected option in this field.
        print "<b>Facility Name:</b> ", facility_name

        
        if shot_number != (''):
            # If there is an input in this field, show it.
            print "<br><b>Shot Number:</b> ", shot_number

        if project != (''):
            print "<br><b>Project Name: </b> ", project
            
        if blame !=  (''):
            print "<br><b>Blame: </b> ", blame

        if date_string != (''):
            print "<br><b>Date: </b> ", date_string 

        if nozzle != (''):
            print "<br><b>Nozzle: </b> ", nozzle
            
        if shocktube != (''):
            print "<br><b>Shocktube: </b> ", shocktube
            
        if reservoir !=(''):
            print "<br><b>Reservoir: </b> ", reservoir
            
        if diaphragm1 != (''):
            print "<br><b>Diaphragm1: </b> ", diaphragm1
            
        if driver0 != (''):
            print "<br><b>Driver0: </b> ", driver0
            
        if driver1 != (''):
            print "<br><b>Driver1: </b> ", driver1
           
        if driver2 != (''):
            print "<br><b>Driver2: </b> ", driver2
            
        if diaphragm2 != (''):
            print "<br><b>Diaphragm2: </b> ", diaphragm2
            
        if diaphragm3 != (''):
            print "<br><b>Diaphragm3: </b> ", diaphragm3
            
        if acceltube != (''):
            print "<br><b>Acceltube: </b> ", acceltube
            
        if dumptank != (''):
            print "<br><b>Dumptank: </b> ", dumptank
            
        if jackoff != (''):
            print "<br><b>Jackoff: </b> ", jackoff
          
        if piston != (''):
            print "<br><b>Piston: </b> ", piston
          
        if condition != (''):
            print "<br><b>Condition: </b> ", condition
            
        if notes != (''):
            print "<br><b>Notes: </b> ", notes

        # Create two break between the previous and the next line.
        print "<br /><br />"
        
        # Create a horizontal line to separate the text and to end the page.
        print """<hr size="3" noshade="noshade" />"""
        print "<br />"

        # Create a link for user to get back to the original page.
        print """<a href="http://www.mech.uq.edu.au/cgi-bin/sf/run_des.py">
                 <b>Back</b></a>"""
               
        # Create a filename by combining the required string together.
        filename = "%s" % (facility_name)
        filename = filename[:2];          # Dump away the last two words.
        filename += "_"
        filename += "%s" % (shot_number)
        filename += "_"
        filename += "run_des"
        filename += ".xml"

        # upload is the name of the directory where the file will be saved to.
        fname = 'upload/' + filename

        # Create fp file and allow it to write.
        fp = open(fname, "w")


        # Create a string that includes all the input data with xml tags.
        str = "<facility_name>" + facility_name + "</facility_name>" + "\n"
        str += "<shot_number>" + shot_number + "</shot_number>" + "\n"
        str += "<project>" + project + "</project>" + "\n"
        str += "<blame>" + blame + "</blame>" + "\n"
        str += "<date_string>" + date_string + "</date_string>" + "\n"
        str += "<nozzle>" + nozzle + "</nozzle>" + "\n"
        str += "<shocktube>" + shocktube + "</shocktube>" + "\n"
        str += "<reservoir>" + reservoir + "</reservoir>" + "\n"
        str += "<diaphragm1>" + diaphragm1 + "</diaphragm1>" + "\n"
        str += "<driver0>" + driver0 + "</driver0>" + "\n"
        str += "<driver1>" + driver1 + "</driver1>" + "\n"
        str += "<driver2>" + driver2 + "</driver2>" + "\n"
        str += "<diaphragm2>" + diaphragm2 + "</diaphragm2>" + "\n"
        str += "<diaphragm3>" + diaphragm3 + "</diaphragm3>" + "\n"
        str += "<acceltube>" + acceltube + "</acceltube>" + "\n"
        str += "<dumptank>" + dumptank + "</dumptank>" + "\n"
        str += "<jackoff>" + jackoff + "</jackoff>" + "\n"
        str += "<piston>" + piston + "</piston>" + "\n"
        str += "<condition>" + condition + "</condition>" + "\n"
        str += "<notes>" + notes + "</notes>" + "\n"

        
        fp.write(str)     # Write the string into the file.
        fp.close()        # Manual close the file.
        

    else:
        # No data means this is the first access so output the form.
        
        try:
            # Open the HTML file that contains the header of this page.
            fp=open("/var/www/html/sf/run_des.html", "r")
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


print "</body></html>"    # End the script.
    
        

        

        

        

    


