def doarm():
    f=open("/dev/tranrec","r+")
    f.write('2')
    f.close()
    print "Trigger units armed"
def dotrigger():
    f=open("/dev/tranrec","r+")
    f.write('3')
    f.close()
    print "Trigger units triggered"
def doreadtr():
    f=open("/dev/tranrec","r+")
    f.write('4')
    f.close()
    f=open("/dev/tranrec","r+")
    temp=f.readline()
    f.close()
    print "Inputted data: 0x%s" % (temp[0:8])
    print "Sampling:        ", (int(temp[0:2],16)&0x80)/0x80
    print "\nTrig 1 slope:    ", (int(temp[0:2],16)&0x02)/0x02
    print "Trig 1 coupling: ", (int(temp[0:2],16)&0x01)/0x01
    print "Trigger unit 1 Trigger Level: ", int(temp[2:4],16)
    print "\nTrig 2 slope:    ", (int(temp[0:2],16)&0x08)/0x08
    print "Trig 2 coupling: ", (int(temp[0:2],16)&0x04)/0x04
    print "Trigger unit 2 Trigger Level: ", int(temp[4:6],16)
    print "\nTrig 3 slope:    ", (int(temp[0:2],16)&0x20)/0x20
    print "Trig 3 coupling: ", (int(temp[0:2],16)&0x10)/0x10
    print "Trigger unit 3 Trigger Level: ", int(temp[6:8],16)
def doreadtb():
    f=open("/dev/tranrec","r+")
    f.write('5')
    f.close()
    f=open("/dev/tranrec","r+")
    temp=f.readline()
    f.close()
    print "Inputted data: 0x%s" % (temp[0:12])
    for i in range(3):
        print "\nTimebase %d:" % (i+1)
        print "Buffer size", 2**((int(temp[2*i:2*i+2],16)&0xc0)/0x40)
        print "Triggering No.", ((int(temp[2*i:2*i+2],16)&0x30)/0x10)+1
        print "Timebase multiplier", (((int(temp[2*i:2*i+2],16)&0x08)/0x08)*99)+1
        n=(int(temp[2*i:2*i+2],16)&0x07)
        print "Sample period %d us" % int(round(0.3125*n**4 - 5.2917*n**3 + 33.562*n**2 - 98.119*n + 119.5)) #quartic regression. should probably be a switch statement to reduce confusion
        print "Pre Trigger Delay", temp[6+2*i:8+2*i]
def doreadads():
    f=open("/dev/tranrec","r+")
    f.write('6')
    f.close()
    f=open("/dev/tranrec","r+")
    temp=f.readline()
    f.close()
    print "Inputted data: 0x%s" % (temp[0:4])
    print "Sampling: " , (int(temp[0:2],16)&0x80)/0x80
    print "Timebase: " , (int(temp[2:4],16)&0x60)/0x20
    print "Buffer Pointer: ", 256*(int(temp[2:4],16)&0x1F)+int(temp[0:2],16)
def doreadadbuf():
    f=open("/dev/tranrec","r+")
def dodirectread(port):
    f=open("/dev/tranrec","r+")
    f.write('8')
    f.close()
    f=open("/dev/tranrec","r+")
    f.write("%d" % port)
    f.close()
    f=open("/dev/tranrec","r+")
    temp=f.readline()
    f.close()
    print "0x%x" % temp
def dodirectwrite(port, data):
    f=open("/dev/tranrec","r+")
    f.write('8')
    f.close()
    f=open("/dev/tranrec","r+")
    f.write("%d" % port)
    f.close()
    f=open("/dev/tranrec","r+")
    f.write("%d" % data)
    f.close()
    print "0x%x written to port 0x%x" % (data, port)
def docardselect(channel, card):
    f=open("/dev/tranrec","r+")
    f.write('11')
    f.close()
    f=open("/dev/tranrec","r+")
    temp = "%d" % (16*channel+card)
    f.write(temp)
    f.close()
def dohelp():
    print """Valid commands:
    help:                   Displays this list
    quit:                   Quits the program
    arm:                    Arms the trigger units
    trigger:                Triggers the trigger units
    read trigger:           Reads the trigger unit status registers
    read timebase:          Reads the timebase unit status registers
    read status:            Reads the S register of the current A/D card
    read direct PORT:       Reads the byte at PORT
    write direct PORT DATA: Writes the byte DATA to port PORT
    select CHANNEL CARD:    Selects channel CHANNEL on data acquisition card CARD"""
    
quit=0
while(not quit):
    n=raw_input(":: ")
    n=n.lower()
    if (n=="help" or n=='?'):
        dohelp()
    elif (n=="arm"):
        doarm()
    elif (n=="trigger"):
        dotrigger()
    elif (n=="read trigger"):
        doreadtr()
    elif (n=="read timebase"):
        doreadtb()
    elif (n=="read status"):
        doreadads()
    elif (n[0:11]=="read direct"):
        a=n.split()
        dodirectread(int(a[2]))
    elif (n[0:12]=="write direct"):
        a=n.split()
        dodirectwrite(int(a[2]),int(a[3]))
    elif (n[0:6]=="select"):
        a=n.split()
        if (int(a[1]<1) or int(a[1]>3)):
            print "Error: Invalid channel"
        else:
            if (int(a[2]<1) or int(a[2]>7)):
                print "Error: Invalid card"
            else:
                docardselect(int(a[1]),int(a[2]))
    elif (n=='q' or n=="quit"):
        print "Goodbye"
        quit=1
    else:
        print "Syntax error: Unknown command"
