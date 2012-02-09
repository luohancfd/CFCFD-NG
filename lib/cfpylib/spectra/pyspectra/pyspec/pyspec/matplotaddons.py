import math
import pylab
import matplotlib
import os
import tempfile

__version__   = "$Revision: 176 $"
__author__    = "Stuart B. Wilkins <stuwilkins@mac.com>"
__date__      = "$LastChangedDate: 2011-02-14 05:08:16 +1000 (Mon, 14 Feb 2011) $"
__id__        = "$Id: matplotaddons.py 176 2011-02-13 19:08:16Z stuwilkins $"

def printfig(*args, **kwargs):
  """Print pyplot figure.

  Identical to savefig, but calles the system printer (unix) to print the current plot.
  specifying the keyword argument printer = 'name' will send to a printer which is not
  the system default. """

  if kwargs.has_key['printer']:
    printer = '-d' + kwargs['printer']
    del kwargs['printer']
  else:
    printer = ''

  fo = tempfile.NamedTemporaryFile(delete = False)
  matplotlib.pyplot.savefig(fo, *args, **kwargs)
  fo.close()
  os.system('lp %s%s' % (printer ,fo.name))

class AnnoteFinder:
  """
  callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.
    
  Register this function like this:
    
  scatter(xdata, ydata)
  af = AnnoteFinder(xdata, ydata, annotes)
  connect('button_press_event', af)
  """

  def __init__(self, xdata, ydata, annotes, axis=None, xtol=None, ytol=None):
    self.data = zip(xdata, ydata, annotes)
    if xtol is None:
      xtol = ((max(xdata) - min(xdata))/float(len(xdata)))/2
    if ytol is None:
      ytol = ((max(ydata) - min(ydata))/float(len(ydata)))/2
    self.xtol = xtol
    self.ytol = ytol
    if axis is None:
      self.axis = pylab.gca()
    else:
      self.axis= axis
    self.drawnAnnotations = {}
    self.links = []
    
  def distance(self, x1, x2, y1, y2):
    """
    return the distance between two points
    """
    return(math.sqrt( (x1 - x2)**2 + (y1 - y2)**2 ))

  def __call__(self, event):
    if event.inaxes:
      clickX = event.xdata
      clickY = event.ydata
      if (self.axis is None) or (self.axis==event.inaxes):
        annotes = []
        for x,y,a in self.data:
          if  (clickX-self.xtol < x < clickX+self.xtol) and  (clickY-self.ytol < y < clickY+self.ytol) :
            annotes.append((self.distance(x,clickX,y,clickY),x,y, a) )
        if annotes:
          annotes.sort()
          distance, x, y, annote = annotes[0]
          self.drawAnnote(event.inaxes, x, y, annote) 
          for l in self.links:
            l.drawSpecificAnnote(annote)
    
  def drawAnnote(self, axis, x, y, annote):
    """
    Draw the annotation on the plot
    """
    if self.drawnAnnotations.has_key((x,y)):
      markers = self.drawnAnnotations[(x,y)]
      for m in markers:
        m.set_visible(not m.get_visible())
      self.axis.figure.canvas.draw()
    else:
      t = axis.text(x,y, "(%3.2f, %3.2f) - %s"%(x,y,annote), )
      m = axis.scatter([x],[y], marker='d', c='r', zorder=100)
      self.drawnAnnotations[(x,y)] =(t,m)
      self.axis.figure.canvas.draw()

  def drawSpecificAnnote(self, annote):
    annotesToDraw = [(x,y,a) for x,y,a in self.data if a==annote]
    for x,y,a in annotesToDraw:
      self.drawAnnote(self.axis, x, y, a)
