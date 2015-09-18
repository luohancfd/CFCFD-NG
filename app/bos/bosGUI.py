#! /usr/bin/env python
""" bosgui.py

Starter GUI for BOS functions

.. Authors:: PJ, 17-Dec-05 starting sample
             Dwishen Ramanah, 12-Jan-2007 Version 1.0
"""

import sys
import Tkinter, tkFileDialog, tkSimpleDialog, tkMessageBox
import Image, ImageTk, ImageChops, ImageOps
import numpy
try:
    from scipy.signal.signaltools import correlate2d
    from scipy.signal.signaltools import medfilt
except:
    print "scipy did not load correctly..."
    print "To get correlation2d and medfilt functions, we need Scipy."

#---------------------------------------------------------------
# Global data
refImage = None; refCanvas = None
x0 = 400; y0 = 400; x1 = 420; y1 = 420

# upper-left and lower-right coordinates of the sample box
# The numbers above suit one of Sreekanth's images.

testImage = None; testCanvas = None
diffImage = None; diffCanvas = None

#---------------------------------------------------------------

def main():
    global root
    root = Tkinter.Tk()
    setupGUImenus(root)
    layoutMainWindow(root)
    root.mainloop()
    return

def layoutMainWindow(root):
    """
    The main window consists of a menubar and some labeled entry boxes.
    """
    global refFileName, testFileName
    global xshiftVar, yshiftVar, sampleBoxVar, sxyVar, wxyVar
    global lowerLevelVar, upperLevelVar, numcontourVar, modelgreylevelVar
    global goodvectorsVar, colourlevelVar
    global x0, y0, x1, y1
    
    refFileName = Tkinter.StringVar()
    testFileName = Tkinter.StringVar()
    sampleBoxVar = Tkinter.StringVar()
    sampleBoxVar.set("%d %d %d %d" % (x0, y0, x1, y1))
    
    xshiftVar = Tkinter.IntVar(); xshiftVar.set(0)
    yshiftVar = Tkinter.IntVar(); yshiftVar.set(0)
    
    wxyVar = Tkinter.IntVar(); wxyVar.set(32)
    sxyVar = Tkinter.IntVar(); sxyVar.set(16)

    lowerLevelVar = Tkinter.DoubleVar(); lowerLevelVar.set(-2)
    upperLevelVar = Tkinter.DoubleVar(); upperLevelVar.set(2)
    numcontourVar =  Tkinter.IntVar(); numcontourVar.set(80)
    modelgreylevelVar =  Tkinter.DoubleVar(); modelgreylevelVar.set(-1.9)
    goodvectorsVar = Tkinter.DoubleVar(); goodvectorsVar.set(0.7)
    colourlevelVar = Tkinter.DoubleVar(); colourlevelVar.set(10)

    #
    root.title("BOS Image Processing")
    f = Tkinter.Frame(root, borderwidth=1)
    f.pack(expand=1, fill="both")
    #
    lab = Tkinter.Label(f, text="Ref image")
    entry = Tkinter.Entry(f, width=60, textvariable=refFileName)
    lab.grid(row=0, column=0, sticky='e')
    entry.grid(row=0, column=1, sticky='w')
    #
    lab = Tkinter.Label(f, text="Test image")
    entry = Tkinter.Entry(f, width=60, textvariable=testFileName)
    lab.grid(row=1, column=0, sticky='e')
    entry.grid(row=1, column=1, sticky='w')
    #
    lab = Tkinter.Label(f, text="xshift")
    entry = Tkinter.Entry(f, width=5, textvariable=xshiftVar)
    lab.grid(row=2, column=0, sticky='e')
    entry.grid(row=2, column=1, sticky='w')
    #
    lab = Tkinter.Label(f, text="yshift")
    entry = Tkinter.Entry(f, width=5, textvariable=yshiftVar)
    lab.grid(row=3, column=0, sticky='e')
    entry.grid(row=3, column=1, sticky='w')
    #
    lab = Tkinter.Label(f, text="Sample box")
    entry = Tkinter.Entry(f, width=20, textvariable=sampleBoxVar)
    lab.grid(row=4, column=0, sticky='e')
    entry.grid(row=4, column=1, sticky='w')
    #
    lab = Tkinter.Label(f, text="Window size")
    entry = Tkinter.Entry(f, width=5, textvariable=wxyVar)
    lab.grid(row=5, column=0, sticky='e')
    entry.grid(row=5, column=1, sticky='w')
    #
    lab = Tkinter.Label(f, text="Search size")
    entry = Tkinter.Entry(f, width=5, textvariable=sxyVar)
    lab.grid(row=6, column=0, sticky='e')
    entry.grid(row=6, column=1, sticky='w')
    #
    lab = Tkinter.Label(f, text="Good Vectors (0-1)")
    entry = Tkinter.Entry(f, width=5, textvariable=goodvectorsVar)
    lab.grid(row=7, column=0, sticky='e')
    entry.grid(row=7, column=1, sticky='w')
    #
    lab = Tkinter.Label(f, text="Lower level value for filled contours")
    entry = Tkinter.Entry(f, width=5, textvariable=lowerLevelVar)
    lab.grid(row=8, column=0, sticky='e')
    entry.grid(row=8, column=1, sticky='w')
    #
    lab = Tkinter.Label(f, text="Upper level value for filled contours")
    entry = Tkinter.Entry(f, width=5, textvariable=upperLevelVar)
    lab.grid(row=9, column=0, sticky='e')
    entry.grid(row=9, column=1, sticky='w')    
    #
    lab = Tkinter.Label(f, text="Number of contour levels")
    entry = Tkinter.Entry(f, width=5, textvariable=numcontourVar)
    lab.grid(row=10, column=0, sticky='e')
    entry.grid(row=10, column=1, sticky='w')
    #
    lab = Tkinter.Label(f, text="Model Grey Level")
    entry = Tkinter.Entry(f, width=5, textvariable=modelgreylevelVar)
    lab.grid(row=11, column=0, sticky='e')
    entry.grid(row=11, column=1, sticky='w')
    #
    lab = Tkinter.Label(f, text="Colour Level")
    entry = Tkinter.Entry(f, width=5, textvariable=colourlevelVar)
    lab.grid(row=12, column=0, sticky='e')
    entry.grid(row=12, column=1, sticky='w')
    #    
    return

def setupGUImenus(root):
    "Establish the menu options and their function calls"
    menuBar = Tkinter.Menu(root)
    #
    fileMenu = Tkinter.Menu(menuBar)
    fileMenu.add_command(label="Open Reference Image", command=openReferenceImage)
    fileMenu.add_command(label="Open Test Image", command=openTestImage)
    fileMenu.add_command(label="Save Difference Image", command=saveDifferenceImage)
    fileMenu.add_separator()
    fileMenu.add_command(label="Quit", command=quitDialog)
    menuBar.add_cascade(label="File", menu=fileMenu)
    #
    actionMenu = Tkinter.Menu(menuBar)    
    actionMenu.add_command(label="1. Estimate Coarse Shift",
                           command=estimateCoarseShift)
    actionMenu.add_command(label="2. Apply Shift to Test Image",
                           command=shiftTestImage)
    actionMenu.add_command(label="3. Compute Difference Image",
                           command=showImageDifference)
    actionMenu.add_command(label="4. Determine horizontal and vertical displacement vectors",
			   command=showVectorDisplacements)
    actionMenu.add_command(label="5. Plot Histogram to represent frequency distribution of horizontal and vertical displacement values",
			   command=histogram_u_v)    
    actionMenu.add_command(label="6. Show Vertical Knife-Edge Schlieren (Contour Filled Plot of horizontal vectors)",
			   command=showContours_u)
    actionMenu.add_command(label="7. Show Horizontal Knife-Edge Schlieren (Contour Filled Plot of vertical vectors)",
			   command=showContours_v)
    actionMenu.add_command(label="8. Show Colour Schlieren",
			   command=showColourSchlieren)
    menuBar.add_cascade(label="Action", menu=actionMenu)
    #
    root.config(menu=menuBar)
    return

#--------------------------------------------------------------------

def quitDialog():
    if tkMessageBox.askokcancel("Quit","Really quit?"):
        sys.exit()

def set_x0_y0_via_mouse(event):
    global sampleBoxVar, x0, y0, x1, y1
    canvas = event.widget
    x1 = x0 = canvas.canvasx(event.x)
    y1 = y0 = canvas.canvasy(event.y)
    sampleBoxVar.set("%d %d %d %d" % (x0, y0, x1, y1))
    canvas.delete("sample-box")
    try:
        testCanvas.delete("sample-box")
    except:
        pass
    return

def set_x1_y1_via_mouse(event):
    global sampleBoxVar, x0, y0, x1, y1
    canvas = event.widget
    x1 = canvas.canvasx(event.x)
    y1 = canvas.canvasy(event.y)
    if (x1 - x0) % 2 != 0: x1 += 1 # make even size
    if (y1 - y0) % 2 != 0: y1 += 1

    dx = x1 - x0  #make sqaure samplebox
    dy = y1 - y0
    if dx > dy: dy=dx
    if dy > dx: dx=dy
    y1 = y0 + dy
    x1 = x0 + dx
    
    sampleBoxVar.set("%d %d %d %d" % (x0, y0, x1, y1))
    canvas.delete("sample-box")
    canvas.create_line(x0, y0, x1, y0, x1, y1, x0, y1, x0, y0,
                       fill="yellow", tags="sample-box")
    try:
        testCanvas.delete("sample-box")
        testCanvas.create_line(x0, y0, x1, y0, x1, y1, x0, y1, x0, y0,
                               fill="yellow", tags="sample-box")
    except:
        pass # the testImage canvas may not exist
    return

def openReferenceImage():    
    """
    Interactive dialog to select the reference image.
    """
    global refImage, refCanvas, refFileName
    import os
    fileName = tkFileDialog.askopenfilename()
    if fileName:
        refImage = Image.open(fileName).convert('L') # to grayscale
        
        print "Finished reading image from file", fileName
        refCanvas = openImageWindow(refImage, 'BOS Reference Image-Filename: '+os.path.basename(fileName), refCanvas)
        refFileName.set(fileName)
    # We will bind a couple of mouse events to the reference canvas so that
    # we can identify the sample area used in the cross-correlation.
    refCanvas.bind('<ButtonPress-1>', set_x0_y0_via_mouse) # anchor box
    refCanvas.bind('<B1-Motion>', set_x1_y1_via_mouse) # draw box
    refCanvas.bind('<ButtonRelease-1>', set_x1_y1_via_mouse) # draw box, final    
    return

def openTestImage():
    """
    Interactive dialog to select the test image.
    """
    global testImage, testCanvas, testFileName
    import os
    fileName = tkFileDialog.askopenfilename()
    if fileName:
        testImage = Image.open(fileName).convert('L') # to grayscale
        print "Finished reading image from file", fileName
        testCanvas = openImageWindow(testImage, 'BOS Test Image-Filename: '+os.path.basename(fileName), testCanvas)
        testFileName.set(fileName)
    return

def saveDifferenceImage():
    """
    Interactive dialog to select the test image.
    """
    global diffImage
    if diffImage:
        fileName = tkFileDialog.asksaveasfilename()
        if fileName:
            testImage = diffImage.save(fileName, 'JPEG')
    return


def openImageWindow(image, title, canvas):
    """
    Displays the supplied image on a scrolled canvas in a new top-level window.

    Returns a reference to the canvas.
    This allows the canvas to reused next time and it allows access to the
    cursor coordinates on the canvas.
    """
    if canvas == None:
        # The canvas does not yet exist.
        t = Tkinter.Toplevel()
        t.title(title)
        f = Tkinter.Frame(t, borderwidth=1)
        f.pack(expand=1, fill="both")
        # We want a large canvas that can be scrolled.
        # Since the Python Megawidgets version seems to be broken,
        # build our own with standard Tkinter widgets.
        scf = Tkinter.Frame(f)
        scf.pack(expand=1, fill="both")
        sh = Tkinter.Scrollbar(scf, orient=Tkinter.HORIZONTAL)
        sv = Tkinter.Scrollbar(scf)
        canvasWidth, canvasHeight = image.size
        # The screen view of the canvas with be fixed, smaller than the image.
        canvas = Tkinter.Canvas(scf, width=800, height=800,
                                scrollregion=(0,0,canvasWidth,canvasHeight) )
        canvas.grid(row=0, column=0, sticky="nsew")
        sv.grid(row=0, column=1, sticky="ns")
        sh.grid(row=1, column=0, sticky="ew")
        sv.config(command=canvas.yview)
        sh.config(command=canvas.xview)
        canvas.config(yscrollcommand=sv.set, xscrollcommand=sh.set)
    # Put the image onto the canvas.
    canvas.im = ImageTk.PhotoImage(image)
    # this c.im reference lasts as long as the canvas
    canvas.create_image(0, 0, image=canvas.im, anchor="nw")
    return canvas

def locateMaximum(mat):
    
    rowsize, colsize = mat.shape
    biggest = 0; ibig = 0; jbig = 0
    for i in range(rowsize):
        for j in range(colsize):
            if mat[i,j] > biggest:
                biggest = mat[i,j]
                ibig = i; jbig = j
                
    return jbig, ibig, biggest

def ind2sub(row,col,ind):
    """ Converts row,col indices into one index for .flat """
    i = ind/col 
    j = ind - i* col
    return i,j  


def gausspeak(cmap,start):
    from scipy.stats import tvar as var
    from scipy.stats.distributions import norm
    from scipy import cumsum,log,exp,zeros,r_,dot,transpose,array
    from scipy import linalg,nan

    #(pseudo-)Gaussian peak fit: p(x,y) = exp(a + b*x + c*y + d*x*x + e*x*y + f*y*y)
    global QMAT
    rhs = zeros((9), dtype=float)
    ind= 0
    for j in r_[-1:2]:
        for i in r_[-1:2]:            
            rhs[ind] = log(max(cmap[start[1]+j,start[0]+i],1e-12))
            ind = ind + 1       
    if var(rhs) == 0:
        
    #right hand side was clipped everywhere: no good
        peak= 0
        shift= array([0,0])
        ok= -1

    #solve normal equations for least squares problem
    qr= dot(transpose(QMAT),transpose(rhs))
    coeffs= linalg.solve(QQ,qr)  
   
    #unpack solution vector; find peak position from zero derivative
    mmat= array([[2*coeffs[3],coeffs[4]],[coeffs[4],2*coeffs[5]]])
    mrhs= array([[-coeffs[1]],[-coeffs[2]]])
    qvec= transpose(linalg.solve(mmat,mrhs))  #qvec = (mmat\mrhs)'
    if norm(qvec) > 1 :
    #interpolated displacement is too large: no good
        peak= 0
        shift= array([0,0])
        ok= -1   
    ok= 1
    shift= qvec + start
    qvec = transpose(qvec)
    peak= exp(coeffs[0] + coeffs[1]*qvec[0] + coeffs[2]*qvec[1] + coeffs[3]*qvec[0]*qvec[0] + coeffs[4]*qvec[0]*qvec[1] + coeffs[5]*qvec[1]*qvec[1])
    return peak,shift,ok

def slidesum(a,n,m):
    from scipy import hstack,vstack,cumsum,zeros,r_,array
    #sliding window summation; averaging window size is [n,m]
    na,ma= a.shape
    aa= vstack(((zeros((1,ma+1),dtype=float)),(hstack(((zeros((na,1),dtype=float)),a)))))
    a1= cumsum(aa,0)
    a2 = a1[n:na+1,:]-a1[0:na+1-n,:]
    a3 = cumsum(a2,1)              
    ss = a3[:,m:ma+1]-a3[:,0:ma+1-m]
    return ss

def simplepiv(I1,I2,wxy,mxy,sxy):
    from scipy.stats import tvar
    from scipy.signal.signaltools import correlate2d
    from scipy import shape,reshape,zeros,r_,dot,transpose,real,conj,sqrt
    from scipy.fftpack import fft2,ifft2
    global QMAT, QQ
    import numpy 

    #convert Image to numpy arrays    
    I1data = numpy.array(I1.getdata())
    I2data = numpy.array(I2.getdata())
    xsize ,ysize = I1.size
    I1data = reshape(I1data,(ysize,xsize))
    I2data = reshape(I2data,(ysize,xsize))

    QMAT = numpy.array([[1,-1,-1,1,1,1],[1,0,-1,0,0,1],[1,1,-1,1,-1,1],
                 [1,-1,0,1,0,0],[1,0,0,0,0,0],[1,1,0,1,0,0],
                 [1,-1,1,1,-1,1],[1,0,1,0,0,1],[1,1,1,1,1,1]])
    QQ = dot(transpose(QMAT),QMAT)


    #-----------------------------------------------------------------------------

    if I1data.shape != I2data.shape:
        print "Images should be same size ..."
    if wxy < 1:
        print "Wrong window size wxy ..."
    if mxy < 1:
        print "Search range mxy too small ..."
    if 2*mxy >= wxy:
        print "Search range mxy too large ..."
    if sxy < 1:
        print "Wrong window spacing sxy ..."        
    height,width = I1data.shape   
    #start and end coordinates for interrogation
    left= 0
    right= width - wxy
    top= 0
    bot= height - wxy
    print top, bot, left, right, sxy
    # size of correlation map, always odd
    CSZ= 2*mxy +1
    #size of sliding template
    wsz= wxy-CSZ+1
    
    x = zeros((int((bot-top)/sxy)+1,int((right-left)/sxy)+1),dtype=float)
    y = zeros((int((bot-top)/sxy)+1,int((right-left)/sxy)+1),dtype=float)
    u = zeros((int((bot-top)/sxy)+1,int((right-left)/sxy)+1),dtype=float)
    v = zeros((int((bot-top)/sxy)+1,int((right-left)/sxy)+1),dtype=float)
    cor = zeros((int((bot-top)/sxy)+1,int((right-left)/sxy)+1),dtype=float)
    vld = zeros((int((bot-top)/sxy)+1,int((right-left)/sxy)+1),dtype=float)

    nx= -1
    for xoff in r_[left:right+1:sxy]:
        nx = nx +1
        ny = -1
        for yoff in r_[top:bot+1:sxy]:
            ny = ny + 1

            A= I1data[yoff+mxy:yoff-mxy+wxy,xoff+mxy:xoff-mxy+wxy]
            B= I2data[yoff:yoff+wxy,xoff:xoff+wxy]
            sumA= sum(sum(A))
            varA= tvar(A)*(wsz*wsz-1)
            sumB= slidesum(B,wsz,wsz)
            varB= slidesum(B**2,wsz,wsz)-(sumB**2)/(wsz*wsz)                           
          
            #compute covariance for windows           
            fft_A= fft2(A,shape=(wxy,wxy)) 
            fft_B= fft2(B,shape=(wxy,wxy))
            corrAB= real(ifft2(conj(fft_A)*fft_B))
            covAB = corrAB[0:CSZ,0:CSZ] - (sumB*sumA)/(wsz*wsz)
            #compute normalized cross correlation coefficient
            cmap= covAB/sqrt(varA*varB)
            xpeak, ypeak, corrbig = locateMaximum(cmap)            
            #compute pattern shift based on correlation peak position
            utemp= xpeak - mxy
            vtemp= ypeak - mxy
            #check for "edge peak": no sub-pixel interpolation possible (or fast processing)

            if xpeak == 0 or xpeak == CSZ-1  or  ypeak==0 or ypeak == CSZ-1:
                valid= 0
            else  :
                
                #subpixel interpolation
                peak,shift,ok = gausspeak(cmap,[xpeak,ypeak])
      
                if ok != 1 or peak/corrbig < 0.98 or peak/corrbig > 2.0:
                    #no good: no subpixel fit available, saddle point topology or inconsistent peak height
                    valid= 0
                else:
                    #good data: use interpolated shifts and new correlation peak value
                    valid= 1
                    corrbig = peak[0]
                    utemp= shift[0,0] - mxy
                    vtemp= shift[0,1] - mxy

            #compute pattern shift based on correlation peak position
            x[ny,nx]= xoff + wxy/2 - 0.5
            y[ny,nx]= yoff + wxy/2 - 0.5
            u[ny,nx]= utemp
            v[ny,nx]= vtemp
            cor[ny,nx]= corrbig
            vld[ny,nx]= valid
    return x,y,u,v,cor,vld

def estimateCoarseShift():
    from scipy import floor
    global testImage, refImage, xshiftVar, yshiftVar
    global sampleBoxVar
    text = sampleBoxVar.get()
    x0, y0, x1, y1 = text.split()
    sampleBox = tuple([int(x0), int(y0), int(x1), int(y1)])
    print "sampleBox=", sampleBox
    refSmallImage = refImage.crop(sampleBox).convert('L')
    testSmallImage = testImage.crop(sampleBox).convert('L')
    # refSmallImage.show(); testSmallImage.show()
    xsize, ysize = refSmallImage.size
    # Convert to numpy arrays and do cross correlation.
    refData = numpy.array(refSmallImage.getdata())
    refData.shape = refSmallImage.size
    testData = numpy.array(testSmallImage.getdata())
    testData.shape = testSmallImage.size
    try:
        # Requires Scipy.
        corrData = correlate2d(refData,testData)
        # print "corrData=", corrData
        ibig, jbig, corrbig = locateMaximum(corrData)
        # Might it be better to shift by a fraction of a pixel?
        # This is possible by doing bilinear interpolation, say.
        print "Max correlation=", corrbig, "  at i=", ibig, "j=", jbig
        xshift, yshift = xsize - 1 - ibig, ysize - 1 - jbig
       # xshift, yshift = ibig - xsize/2 + 1, jbig - ysize/2 + 1
    except:
        xshift, yshift = 0, 0
    xshiftVar.set(xshift); yshiftVar.set(yshift)
    return

    
def shiftTestImage():
    global testImage, testCanvas
    xsize, ysize = testImage.size
    xshift = int(xshiftVar.get())
    yshift = int(yshiftVar.get())
    shiftedImage = testImage.crop((xshift, yshift, xsize, ysize))
    testCanvas.im = ImageTk.PhotoImage(shiftedImage)
    testCanvas.create_image(0, 0, image=testCanvas.im, anchor="nw")
    testImage = shiftedImage
    print "After shifting, new testImage size:", testImage.size
    return

def showImageDifference():
    global testImage, refImage, diffImage, diffCanvas,croppedRefImage 
    xsize, ysize = testImage.size
    croppedRefImage = refImage.crop((0, 0, xsize, ysize))
    diffImage = ImageChops.difference(testImage, croppedRefImage)
    diffCanvas = openImageWindow(diffImage, 'BOS Difference', diffCanvas)
    return

def showVectorDisplacements():
    
    global testImage, croppedRefImage, u, v,valid,q1,umean,vmean,x,y,sxyVar,wxyVar, goodvectorsVar  
    from scipy import where,compress,logical_and,median,logical_or,nan
    from pylab import resize,transpose,quiver,title,show,find,imshow,hist,figure,clf,draw,save,load,xlabel,ylabel,flipud
        
    mxy = 3     
    wxy = int(wxyVar.get())
    sxy = int(sxyVar.get())
    goodvectors = float(goodvectorsVar.get())   
    #process to find PIV-style displacements
    x,y,u,v,q1,valid = simplepiv(croppedRefImage,testImage,wxy,mxy,sxy)    
    good = where(logical_and(q1 > goodvectors,valid > 0),True,False)
    umean= median(compress(good.flat,u.flat))
    vmean= median(compress(good.flat,v.flat))    
    u = where(logical_or(q1 < goodvectors,valid < 0),0,u)
    v = where(logical_or(q1 < goodvectors,valid < 0),0,v)
    u = u-umean
    v = v-vmean
    save('vecx.out',x)
    save('vecy.out',y)
    save('vecu.out',u)
    save('vecv.out',v)
    save('vecq1.out',q1)
    save('vecvalid.out',valid)
    u = flipud(u)
    v = -flipud(v)    
    quiver (x,y,u,v)        
    title ('Vector displacements')
    xlabel('Pixels')
    ylabel('Pixels')
    show()
    return

def histogram_u_v():
   
    from pylab import figure,hist,title,load,xlabel,ylabel
    from scipy import where,logical_or,flipud
    
    u = load('vecu.out.npy')
    v = load('vecv.out.npy')
    u = flipud(u)
    v = -flipud(v)       
    figure()
    hist(u,100)
    title ('Frequency Distribution of horizontal displacements data')
    xlabel('Horizontal Displacement in Pixels')
    ylabel('Frequency')
    figure()
    hist(v,100)
    title ('Frequency Distribution of vertical displacements data')
    xlabel('Vertical Displacement in Pixels')
    ylabel('Frequency')
    
def showContours_u():

    global q1,valid,umean,vmean,upperLevelVar,lowerLevelVar,numcontourVar,goodvectorsVar,modelgreylevelVar,lowerLevel,upperLevel,numcontour,tick
    from scipy import zeros,where,logical_or,argmin,shape,ravel,nan,compress,r_,flipud
    from pylab import imshow,clf,title,save,load,figure,contourf,cm,hold,contour,xlabel,ylabel
    

    x = load('vecx.out.npy')
    y = load('vecy.out.npy')
    u = load('vecu.out.npy')
    q1= load('vecq1.out.npy')
    valid=load('vecvalid.out.npy')    
    modelgreylevel = float(modelgreylevelVar.get())
    goodvectors = float(goodvectorsVar.get())    
    u = where(logical_or(q1 < goodvectors,valid < 0),modelgreylevel,u)    
    u = flipud(u)       
    lowerLevel = float(lowerLevelVar.get())    
    upperLevel = float(upperLevelVar.get())
    numcontour = int(numcontourVar.get())    
    tick = float(upperLevel-lowerLevel)/numcontour    
    uu=r_[lowerLevel:upperLevel:tick]    
    figure()    
    contourf(x,y,u,uu,cmap=cm.gray)
    contourf(x,y,u,uu,cmap=cm.gray)    
    xlabel('Pixels')
    ylabel('Pixels')
    return

def showContours_v():

    global umean,vmean,modelgreylevelvar,goodvectorsVar,lowerLevelVar,upperLevelVar,tickvar,numcontourVar
    from scipy import zeros,where,logical_or,r_,argmin,shape,ravel,nan,compress,flipud
    from pylab import imshow,clf,title,save,load,figure,contourf,cm,hold,contour,xlabel,ylabel
    
    x = load('vecx.out.npy')
    y = load('vecy.out.npy')
    u = load('vecu.out.npy')
    v = load('vecv.out.npy')
    q1= load('vecq1.out.npy')
    valid=load('vecvalid.out.npy')    
    modelgreylevel = float(modelgreylevelVar.get())
    goodvectors = float(goodvectorsVar.get())    
    v = where(logical_or(q1 < goodvectors,valid < 0),modelgreylevel,v)
    v = -flipud(v)    
    lowerLevel = float(lowerLevelVar.get())
    upperLevel = float(upperLevelVar.get())
    numcontour = int(numcontourVar.get())
    tick = float(upperLevel-lowerLevel)/numcontour    
    vv=r_[lowerLevel:upperLevel:tick]    
    figure()    
    contourf(x,y,v,vv,cmap=cm.gray)
    contourf(x,y,v,vv,cmap=cm.gray)    
    xlabel('Pixels')
    ylabel('Pixels')
    return

def showColourSchlieren():

    global colourlevelVar
    from scipy import arctan2, shape,r_,pi,sqrt,zeros,shape
    from pylab import contourf,quiver,load,figure,show,imshow
    from colorsys import hsv_to_rgb

    x = load('vecx.out.npy')
    y = load('vecy.out.npy')
    u = load('vecu.out.npy')
    v = load('vecv.out.npy')

    angle=arctan2(v,u)
    colourlevel = float(colourlevelVar.get())
    #convert angle to go from 0 to 2*pi.
    height, width = x.shape
    for i in r_[0:height-1]:
        for j in r_[0:width-1]:
            if angle[i,j]<0:
                angle[i,j]=2*pi+angle[i,j]        

    magnitude=sqrt(u**2+v**2)
    max_magnitude= magnitude.max()

    C = zeros((int(height),int(width),int(3)),dtype=float)
    for i in r_[0:height-1]:
        for j in r_[1:width-1]:
            h=angle[i,j]/2/pi
            s=colourlevel*magnitude[i,j]/max_magnitude
            if s>1:
                s=1
            vv=1
            h,s ,vv = hsv_to_rgb(h,s,vv)
            C[i,j,0]=h
            C[i,j,1]=s
            C[i,j,2]=vv
    figure()     
    imshow(C)
    show()
    return 
#-------------------------------------------------------------------
if __name__ == '__main__':
    main()
