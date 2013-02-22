# Transient Recorder debug program
# GUI version
# Zane Smith
# July 2004

# Destroys the input window
def killInput(event = None):
    global passData
    global iframe
    passData = iframe.value.get()
    inproot.destroy()
    inproot.quit()

# Creates a window to get input from the user
def getInput(prompt):
    global iframe
    global inproot

    inproot = Tk()

    inproot.title(prompt)

    iframe = Frame(inproot)
    iframe.pack()

    iframe.label = Label(inproot, text=prompt).pack(side=TOP,padx=5,pady=10)

    iframe.field = Entry(inproot)
    iframe.field.pack(padx=5,pady=10)

    iframe.value = StringVar(inproot)
    iframe.field["textvariable"] = iframe.value

    iframe.button = Button(inproot, text=" OK ", command = killInput)
    iframe.button.pack(padx=5,pady=10)

    iframe.field.bind('<Key-Return>', killInput)
    iframe.field.focus_set()
    iframe.mainloop()

    return passData

# Destroys the message window
def killMessage(event = None):
    global msgroot
    msgroot.destroy()
    msgroot.quit()

# Creates a window to pass a message to the user
def doMessage(message,title):
    global msgroot
    msgroot = Tk()
    msgroot.title(title)

    msgFrame = Frame(msgroot, width=300, height=30)
    msgFrame.pack(side=TOP,padx=10,pady=10)

    msg = Label(msgFrame, text=message,wraplength=300)
    msg.pack()

    button = Button(msgroot,text=" OK ",defaul=ACTIVE,command=killMessage)
    button.pack(padx=5,pady=10)

    msgroot.bind('<Key-Return>', killMessage)

    msgroot.mainloop()

# Arms the trigger units
def doArm(event=None):
    global race
    race = 1
    f=open("/dev/tranrec","r+")
    f.write('2')
    f.close()
    print "Trigger units armed"
    race = 0

# Triggers the trigger units
def doTrigger(event=None):
    global race
    race = 1
    f=open("/dev/tranrec","r+")
    f.write('3')
    f.close()
    print "Trigger units triggered"
    race = 0

# Displays debug data for the trigger units
def doReadTr(event=None):
    global race
    race = 1
    passMessage = ""
    f=open("/dev/tranrec","r+")
    f.write('4')
    f.close()
    f=open("/dev/tranrec","r+")
    temp=f.readline()
    f.close()
    race = 0
    print "Inputted data: 0x%s" % (temp[0:8])
    if ((int(temp[0:2],16)&0x80)/0x80):
        tempstring="Sampling"
    else:
        tempstring="Valid Data"
    passMessage = passMessage + tempstring
    if ((int(temp[0:2],16)&0x02)/0x02):
        tempstring="Rising"
    else:
        tempstring="Falling"
    passMessage = passMessage + "\n\nTrig 1 slope:    " + tempstring
    if ((int(temp[0:2],16)&0x01)/0x01):
        tempstring="D.C."
    else:
        tempstring="A.C."
    passMessage = passMessage + "\nTrig 1 coupling: " + tempstring
    passMessage = passMessage + "\nTrig 1 Level: %d" % int(temp[2:4],16)
    if ((int(temp[0:2],16)&0x08)/0x08):
        tempstring="Rising"
    else:
        tempstring="Falling"
    passMessage = passMessage + "\n\nTrig 2 slope:    "+ tempstring
    if ((int(temp[0:2],16)&0x04)/0x04):
        tempstring="D.C."
    else:
        tempstring="A.C."
    passMessage = passMessage + "\nTrig 2 coupling: " + tempstring
    passMessage = passMessage + "\nTrig 2 Level: %d" % int(temp[4:6],16)
    if ((int(temp[0:2],16)&0x20)/0x20):
        tempstring="Rising"
    else:
        tempstring="Falling"
    passMessage = passMessage + "\n\nTrig 3 slope:    " + tempstring
    if ((int(temp[0:2],16)&0x10)/0x10):
        tempstring="D.C."
    else:
        tempstring="A.C."
    passMessage = passMessage + "\nTrig 3 coupling: " + tempstring
    passMessage = passMessage + "\nTrig 3 Level:  %d" % int(temp[6:8],16)

    doMessage(passMessage, "Trigger Unit")

# Displays debug data for the timebase units
def doReadTb(event=None):
    global race
    race = 1
    passMessage = ""
    f=open("/dev/tranrec","r+")
    f.write('5')
    f.close()
    f=open("/dev/tranrec","r+")
    temp=f.readline()
    f.close()
    race = 0
    print "Inputted data: 0x%s" % (temp[0:12])
    for i in range(3):
        bufs = "ERROR"
        trign = "ERROR"
        ptd=0
        x=((int(temp[2*i:2*i+2],16)&0xc0)/0x40)
        if x==0:
            bufs = "Armed"
        elif x==1:
            bufs = "2k"
        elif x==3:
            bufs = "4k"
        elif x==2:
            bufs = "8k"
            ptd=1
        x=((int(temp[2*i:2*i+2],16)&0x30)/0x10)
        if x==0:
            trign = "#1"
        elif x==1:
            trign = "#2"
        elif x==3:
            trign = "#3"
        passMessage = passMessage + "\nTimebase %d:" % (i+1)
        passMessage = passMessage + "\nBuffer size " + bufs
        passMessage = passMessage + "\nTriggering No. " + trign
        passMessage = passMessage + "\nTimebase multiplier %d" % ((((int(temp[2*i:2*i+2],16)&0x08)/0x08)*99)+1) #dodgy function. make it if/else
        n=(int(temp[2*i:2*i+2],16)&0x07)
        passMessage = passMessage + "\nSample period %d us" % int(round(0.3125*n**4 - 5.2917*n**3 + 33.562*n**2 - 98.119*n + 119.5)) #quartic regression. should probably be a switch statement to reduce confusion
        if ptd:
            passMessage = passMessage + "\nPre Trigger Delay " + temp[6+2*i:8+2*i]
        passMessage = passMessage + "\n"
    doMessage(passMessage, "Timebase Unit")

# Displays debug information for the A/D buffer
def doReadAdS(event=None):
    global race
    race = 1
    f=open("/dev/tranrec","r+")
    f.write('6')
    f.close()
    f=open("/dev/tranrec","r+")
    temp=f.readline()
    f.close()
    race = 0
    passMessage = "Inputted data: 0x%s" % (temp[0:4])
    passMessage = passMessage + "\nSampling: %d" % ((int(temp[0:2],16)&0x80)/0x80)
    passMessage = passMessage + "\nTimebase: %d" % ((int(temp[2:4],16)&0x60)/0x20)
    passMessage = passMessage + "\nBuffer Pointer: %d" % ( 256*(int(temp[2:4],16)&0x1F)+int(temp[0:2],16))
    doMessage(passMessage,"A/D S Register")

# Reads the contents of the A/D buffer into memory
def doReadAdBuf(event=None):
    global data
    global race
    race = 1
    channel=var.get()
    card=varb.get()
    f=open("/dev/tranrec","r+")
    f.write('7')
    f.close()
    f=open("/dev/tranrec","r+")
    data=f.read()
    f.close()
    race = 0
    low=(ord(data[0])&0xF0)/0x10
    sc=(ord(data[0])&0x0F)
    high=ord(data[1])
    print "high",high,"\nlow",low,"\nsc",sc
    print ((((16.0*high)+low)/2048.0)-1.0)
    dataStatus.set("Data: Card %d, Channel %d" % (card,channel))
    saveStatus.set("Not Saved")

# Saves the current data to a file
def doSaveFile(event=None):
    global data
    if (data==""):
        doMessage("No data to be saved","Error")
    else:
        fname = getInput("Save to file: ")
        f=open(fname,"w")
        f.write(data)
        f.close()
        saveStatus.set("Saved to: %s" % fname)

# Graphs the current data
def doGraph(event=None):
    global data
    if (data==""):
        doMessage("No data to be graphed","Error")
    else:
        graphroot=Tk()
        canvas = Canvas(graphroot, width = 1050, height = 420, bg = "white")
        canvas.pack()
        sampledata = []
        for i in range(1000):
            L = (ord(data[2*i])&0xF0)/0x10
            N = ord(data[2*i])&0x0F
            H = ord(data[(2*i)+1])
            if N==2:
                FS=5
                #coupling="AC"
            elif N==3:
                FS=5
                #coupling="DC"
            elif N==4:
                FS=2
                #coupling="AC"
            elif N==5:
                FS=2
                #coupling="DC"
            elif N==6:
                FS=1
                #coupling="AC"
            elif N==7:
                FS=1
                #coupling="DC"
            elif N==8:
                FS=0.5
                #coupling="AC"
            elif N==9:
                FS=0.5
                #coupling="DC"
            elif N==10:
                FS=0.2
                #coupling="AC"
            elif N==11:
                FS=0.2
                #coupling="DC"
            elif N==12:
                FS=0.1
                #coupling="AC"
            elif N==13:
                FS=0.1
                #coupling="DC"
            elif N==14:
                FS==0
                #coupling="GND"
            sampledata.append((i+50, 215.0-40.0*(((16.0*H+L)/2048.0)-1.0)*FS))
        canvas.create_line([(50,215),(1050,215)],fill = "#000000")
        canvas.create_line([(50,5),(50,415)],fill="#000000")
        for i in range(-5,6,1):
            canvas.create_text(35,215-40.0*i,text='%d'%i,anchor=E)
        canvas.create_line(sampledata, fill = "royalblue")

# Reads from a port directly
def doDirectRead(event=None):
    global race
    race = 1
    ports= getInput("Port: ")
    if (ports[:2]=="0x"):
        port = int(ports[2:],16)
    else:
        port = int(ports)
    f=open("/dev/tranrec","r+")
    f.write('8')
    f.close()
    f=open("/dev/tranrec","r+")
    f.write("%d" % port)
    f.close()
    f=open("/dev/tranrec","r+")
    temp=f.readline()
    f.close()
    race = 0
    passMessage = "0x%x" % int(temp,16)
    doMessage(passMessage,"Direct Read Port 0x%x (%d)" % (port,port))

# Writes to a port directly
def doDirectWrite(event=None):
    global race
    race = 1
    ports= getInput("Port: ")
    if (ports[:2]=="0x"):
        port = int(ports[2:],16)
    else:
        port = int(ports)
    datas= getInput("Data: ")
    if (datas[:2]=="0x"):
        data = int(datas[2:],16)
    else:
        data = int(datas)
    f=open("/dev/tranrec","r+")
    f.write('9')
    f.close()
    f=open("/dev/tranrec","r+")
    f.write("%d" % port)
    f.close()
    f=open("/dev/tranrec","r+")
    f.write("%d" % data)
    f.close()
    race = 0
    print "0x%x written to port 0x%x" % (data, port)

# Selects a card and channel
def doCardSelect(event=None):
    global race
    race = 1
    channel=var.get()
    card=varb.get()
    f=open("/dev/tranrec","r+")
    f.write('11')
    f.close()
    f=open("/dev/tranrec","r+")
    temp = "%d" % (16*channel+card)
    f.write(temp)
    f.close()
    race = 0

# Quits the application
def doQuit(event=None):
    root.destroy()
    root.quit()

# If polling, update the trigger status
def doUpdate(event=None):
    global trigStatus
    global trigLabel
    global race
    global poll
    if (not poll.get()):
        trigStatus.set("Disabled")
        trigLabel.configure(bg="#BBBBBB")
    else:
        if (not race):
            f=open("/dev/tranrec","r+")
            f.write('8')
            f.close()
            f=open("/dev/tranrec","r+")
            f.write("786")
            f.close()
            f=open("/dev/tranrec","r+")
            temp=f.readline()
            f.close()
            bufs = (int(temp,16)&0xC0)/0x40
            if (bufs==0):
                trigStatus.set("Armed")
                trigLabel.configure(bg="#11EE11")
            else:
                trigStatus.set("Triggered")
                trigLabel.configure(bg="#EE1111")
    root.after(1000,doUpdate)



from Tkinter import *

root=Tk()

root.title("Databox Debug Menu")

data=""
passData = 0
frame = Frame(root, width=30, height=300)
frame.pack(side=LEFT)
btns = [("(A)rm","#11EE11",doArm,'a'),
        ("(T)rigger","#EE1111",doTrigger,'t'),
        ("Read tr(i)gger unit","#888888",doReadTr,'i'),
        ("Read ti(m)ebase unit","#888888",doReadTb,'m'),
        ("Read A/(D) register","#888888",doReadAdS,'d'),
        ("Read A/D (b)uffer","#888888",doReadAdBuf,'b'),
        ("(R)ead direct port","#888888",doDirectRead,'r'),
        ("(W)rite direct port","#888888",doDirectWrite,'w'),
        ("(G)raph data","#888888",doGraph,'g'),
        ("(S)ave data","#888888",doSaveFile,'s'),
        ("(Q)uit","#888888",doQuit,'q')]
for name, colour, evnt, key in btns:
    btn = Button(frame, text=name, bg=colour, command=evnt, width = 30)
    root.bind(key,evnt)
    btn.pack(side=TOP,expand=NO, pady=5)



framea = Frame(root)
framea.pack(side=TOP,padx=10,pady=10)

frameb = Frame(framea)
frameb.pack(side=LEFT,padx=10,pady=10)

varb = IntVar()
varb.set(4)
for text, value in [('Card 1',1),
                    ('Card 2',2),
                    ('Card 3',3),
                    ('Card 4',4),
                    ('Card 5',5),
                    ('Card 6',6),
                    ('Card 7',7)]:
    rdbtn = Radiobutton(frameb, text=text,value=value,variable=varb,command=doCardSelect)
    rdbtn.pack(expand=NO)

framec = Frame(framea)
framec.pack(side=LEFT,padx=10,pady=10)

var = IntVar()
var.set(1)
for text, value in [('Channel 1',1),
                    ('Channel 2',2),
                    ('Channel 3',3)]:
    rdbtn = Radiobutton(framec, text=text,value=value,variable=var,command=doCardSelect)
    rdbtn.pack(expand=NO)

dataStatus = StringVar()
dataLabel = Label(root,textvariable=dataStatus)
dataLabel.pack(side=TOP,padx=10,pady=10)

saveStatus = StringVar()
saveLabel = Label(root,textvariable=saveStatus)
saveLabel.pack(side=TOP,padx=10,pady=10)

trigStatus = StringVar()
trigLabel = Label(root,textvariable = trigStatus)
trigLabel.pack(side=TOP,padx=10,pady=10)

root.after(10,doUpdate)

poll=IntVar()
poll.set(1)

pollCheck = Checkbutton(root,state=NORMAL,text="Poll trigger state",variable=poll)
pollCheck.pack(side=TOP,padx=10,pady=10)

doCardSelect()

root.mainloop()
