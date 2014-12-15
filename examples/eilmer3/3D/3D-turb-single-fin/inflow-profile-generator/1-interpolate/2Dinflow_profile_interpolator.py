#!/usr/bin/env python
# 2Dinflow_profile_interpolator.py
# Accepts 2D high-res slice .dat, low res 3D slice .dat, interpolates values and saves new low res .dat
# Based on compute_viscous_data_simple.py by Peter Jacobs and Wilson Chan
# Created 3 November 2014 by Samuel Stennett

#EXECUTE USING python 2Dinflow_profile_interpolator.py OTHERWISE WILL NOT WORK

print "Begin Profile Interpolator"

def interpolate(b1,b3,a1,a2,a3):	#a values are distances, b values are flow properties
    b2 = ((a2-a3)/(a1-a3))*(b1-b3)+b3;
    return b2;

lowresIN = open("./tc2update-inflow-slice.dat", "r")	#CHANGE AS REQUIRED
highresIN = open("./flat_runup-204mm-slice.dat","r")	#CHANGE AS REQUIRED
line = lowresIN.readline()  # read away the first (comment) line containing notes on format
line2 = highresIN.readline()
# Assuming .dat format is ...
# Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0]

outfile = open("interpolated_profile_ghost2.dat", "w")
outfile.write("# pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu k[0] mu_t k_t S tke omega massf[0] e[0] T[0] \n")


#Read first line as first previous line, second line as first next line.
linePREV = highresIN.readline().strip()
lineNEXT = highresIN.readline().strip()

while True:
    #Extract line from lowres file, check if at end of file.
    found=0
    currentVAL = lowresIN.readline().strip()
    if len(currentVAL) == 0: break
    data_arrayZ = [float(word) for word in currentVAL.split()]
    hcurrent = data_arrayZ[2]
    
    while found==0:
        #Iterate through high res data, searching for surrounding values.
        data_arrayY1 = [float(word) for word in linePREV.split()]
        data_arrayY2 = [float(word) for word in lineNEXT.split()]
        
        hprev = data_arrayY1[1]
        hnext = data_arrayY2[1]

        if hcurrent > hprev and hcurrent <= hnext:
            #Calculate interpolated flow properties.
            print "Match found!"
            xval = 0.0
            yval = 0.0
            zval = hcurrent
            vol = interpolate(data_arrayY1[3],data_arrayY2[3],hprev,hcurrent,hnext)
            rho = interpolate(data_arrayY1[4],data_arrayY2[4],hprev,hcurrent,hnext)
            velx = interpolate(data_arrayY1[5],data_arrayY2[5],hprev,hcurrent,hnext)
            vely = interpolate(data_arrayY1[6],data_arrayY2[6],hprev,hcurrent,hnext)
            velz = interpolate(data_arrayY1[7],data_arrayY2[7],hprev,hcurrent,hnext)
            press = interpolate(data_arrayY1[8],data_arrayY2[8],hprev,hcurrent,hnext)
            soundspeed = interpolate(data_arrayY1[9],data_arrayY2[9],hprev,hcurrent,hnext)
            mu = interpolate(data_arrayY1[10],data_arrayY2[10],hprev,hcurrent,hnext)
            kval = interpolate(data_arrayY1[11],data_arrayY2[11],hprev,hcurrent,hnext)
            mut = interpolate(data_arrayY1[12],data_arrayY2[12],hprev,hcurrent,hnext)
            kvalt = interpolate(data_arrayY1[13],data_arrayY2[13],hprev,hcurrent,hnext)
            Sval = interpolate(data_arrayY1[14],data_arrayY2[14],hprev,hcurrent,hnext)
            tke = interpolate(data_arrayY1[15],data_arrayY2[15],hprev,hcurrent,hnext)
            omega = interpolate(data_arrayY1[16],data_arrayY2[16],hprev,hcurrent,hnext)
            mdot = interpolate(data_arrayY1[17],data_arrayY2[17],hprev,hcurrent,hnext)
            intenergy = interpolate(data_arrayY1[18],data_arrayY2[18],hprev,hcurrent,hnext)
            temp = interpolate(data_arrayY1[19],data_arrayY2[19],hprev,hcurrent,hnext)
#            mach = interpolate(data_arrayY1[20],data_arrayY2[20],hprev,hcurrent,hnext)
#            pitot = interpolate(data_arrayY1[21],data_arrayY2[21],hprev,hcurrent,hnext)
            outfile.write("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n" % (xval, yval, zval, vol, rho, velx, vely, velz, press, soundspeed, mu, kval, mut, kvalt, Sval, tke, omega, mdot, intenergy, temp))
            found=1
        else:
            #Continue searching for values either side of selected value.
            print "Still searching for match..."
            linePREV = lineNEXT
            lineNEXT = highresIN.readline().strip()
  
#Close infiles and outfiles.  
lowresIN.close()
highresIN.close()
outfile.close()

print "Done."
