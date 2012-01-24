## \file svg_render.py
## \ingroup geom
## \brief Render a drawing in Scalable Vector Graphics format.
##
## \author P.Jacobs
## \version 23-Oct-2005 first cut to suit scriptit.py in mb_cns
"""
This module provides a few convenient functions for rendering a drawing
in Scalable Vector Graphics format.

The main transformation is from a user-space coordinate system
with (0,0) at the lower-left corner
to the SVG coordinate system with (0,0) in the upper-left corner
of the page.

Along the way, user-space units are converted to points because
Inkscape seems to behave better if everything is specified in points.
"""

from math import sqrt

class SvgEnvironment(object):
    def __init__(self, width=120.0, height=120.0, unitLength="mm",
                 title="Untitled", desc="No description"):
        """
        Creates a SVG environment with a particular canvas size.
        """
        self.width = width
        self.height = height
        self.unitLength = unitLength
        self.unitsToPoints = 1.0
        if unitLength == "in":
            self.unitToPoints = 72.0
        elif unitLength == "mm":
            self.unitToPoints = 72.0 / 25.4
        elif unitLength == "cm":
            self.unitToPoints = 72.0 / 2.54
        elif unitLength == "m":
            self.unitToPoints = 72.0 / 0.0254
        else:
            raise ValueError, "SvgEnvironment: Unknown units %s" % unitLength
        # print "unitToPoints=", self.unitToPoints
        self.fileHandle = None
        self.lineColour = "black"
        self.fillColour = "none"
        self.setLineWidth(0.25)
        self.lineCap = "round"
        self.lineJoin = "round"
        self.setDashArray()
        self.title = title
        self.description = desc
        return

    def toPointsX(self, x):
        """
        Transforms x-coordinate from user-space to SVG space.
        
        @returns: points for SVG
        """
        return x * self.unitToPoints

    def toPointsY(self, y):
        """
        Transforms y-coordinate from user-space to SVG space.

        @param y: y-coordinate in units (with the origin in the lower-left corner)
        @type y: float or int
        @returns: points in SVG coordinate system (with the origin in the upper left)
        """
        return (self.height - y) * self.unitToPoints

    def setLineWidth(self, w):
        """
        Sets line width.
        
        @param w: line-width in mm.
        @type w: float or int
        """
        self.lineWidth = w * 72.0 / 25.4
        return

    def setDashArray(self, dashLength=2.0, gapLength=2.0):
        """
        Sets length of dashes and gaps.
        
        @param dashLength: in mm
        @type dashLength: float or int
        @param gapLength: in mm
        @type gapLength: float or int
        """
        self.dashLength = dashLength * 72.0/25.4
        self.gapLength = gapLength * 72.0/25.4
        return
    
    def getLineStyle(self, dashed=0):
        """
        Assembles a suitable style specification string.
        
        @param dashed: flag to indicate that the line is dashed
        @type dashed: boolean or int
        """
        style = "stroke:%s;stroke-width:%.2f;stroke-linecap:%s;fill:%s" % \
                (self.lineColour, self.lineWidth, self.lineCap, self.fillColour)
        if dashed:
            style += ";stroke-dasharray: %.2f %.2f" % \
                     (self.dashLength, self.gapLength)
        return style
    
    def open(self, fileName="drawing.svg"):
        """
        Opens the SVG file and writes the preamble.
        """
        f = open(fileName, "w")
        self.fileHandle = f
        f.write("<?xml version=\"1.0\"?>\n")
        f.write("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\"\n")
        f.write("\"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n")
        f.write("<svg width=\"%.2f\" height=\"%.2f\">\n" % \
                (self.width*self.unitToPoints, self.height*self.unitToPoints) )
        f.write("<title>%s</title>\n" % self.title)
        f.write("<desc>%s</desc>\n" % self.description)
        f.write("<!-- drawing begins here -->\n")
        return

    def add_comment(self, text):
        """
        Inserts a comment into the SVG file.
        """
        f = self.fileHandle
        f.write("<!-- ")
        f.write(text)
        f.write(" -->\n")
        return
    
    def close(self):
        """
        Finishes off the SVG file and closes it.
        """
        f = self.fileHandle
        if f:
            f.write("<!-- end of drawing -->\n")
            f.write("</svg>\n")
            f.close()
        return

    def line(self, x1, y1, x2, y2, dashed=0):
        """
        Render a line from point 1 to point 2.
        """
        x1p = self.toPointsX(x1); y1p = self.toPointsY(y1)
        x2p = self.toPointsX(x2); y2p = self.toPointsY(y2)
        f = self.fileHandle
        f.write("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" " \
                 % (x1p, y1p, x2p, y2p) )
        f.write("style=\"%s\" />\n" % self.getLineStyle(dashed=dashed) )
        return

    def polyline(self, xlist, ylist, dashed=0):
        """
        Render a polyline from lists of x and y coordinates.
        """
        x0 = self.toPointsX(xlist[0]); y0 = self.toPointsY(ylist[0])
        f = self.fileHandle
        f.write("<path d=\"M %.2f %.2f" % (x0, y0) )
        for i in range(1, min(len(xlist), len(ylist))):
            x = self.toPointsX(xlist[i]); y = self.toPointsY(ylist[i])
            f.write(", L %.2f %.2f" % (x, y) )
        f.write("\"") # finish the coordinate string
        f.write(" style=\"%s\"" % self.getLineStyle(dashed=dashed) )
        f.write("/>\n")
        return

    def arc(self, x0, y0, x1, y1, xc, yc, dashed=0):
        """
        Render a circular arc from point 0 to pint 1 with centre at c.
        """
        # in user-space coordinates, comput the properties of the arc.
        dx0 = x0 - xc; dy0 = y0 - yc
        dx1 = x1 - xc; dy1 = y1 - yc
        r0 = sqrt(dx0*dx0 + dy0*dy0)
        r1 = sqrt(dx1*dx1 + dy1*dy1)
        assert abs(r0 - r1) < 1.0e-6*(r0+1.0), "Radii don't match"
        crossz = (dx0 * dy1 - dx1 * dy0) # z-component of vector product
        clockwise = (crossz < 0.0)
        # now, do the rendering in SVG space (in points)
        x0p = self.toPointsX(x0); y0p = self.toPointsY(y0)
        x1p = self.toPointsX(x1); y1p = self.toPointsY(y1)
        rp = r0 * self.unitToPoints
        x_axis_rotation = 0.0
        large_arc_flag = 0
        if clockwise:
            sweep_flag = 1
        else:
            sweep_flag = 0
        f = self.fileHandle
        f.write("<path d=\"M %.2f %.2f A %.2f %.2f, %.2f, %d, %d, %.2f %.2f\" " \
                % (x0p, y0p, rp, rp, x_axis_rotation, large_arc_flag,
                   sweep_flag, x1p, y1p) )
        f.write("style=\"%s\"" % self.getLineStyle(dashed=dashed) )
        f.write("/>\n")
        return
    
    def circle(self, x, y, r, dashed=0):
        """
        Render a circle of radius r at centre (x,y).
        """
        xp = self.toPointsX(x); yp = self.toPointsY(y)
        rp = r * self.unitToPoints
        f = self.fileHandle
        f.write("<circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" " \
                % (xp, yp, rp) )
        f.write("style=\"%s\" />\n" % self.getLineStyle(dashed=dashed) )
        return

    def bezier3(self, x0, y0, x1, y1, x2, y2, x3, y3, dashed=0):
        """
        Render a thrid-order Bezier curve.
        """
        x0p = self.toPointsX(x0); y0p = self.toPointsY(y0)
        x1p = self.toPointsX(x1); y1p = self.toPointsY(y1)
        x2p = self.toPointsX(x2); y2p = self.toPointsY(y2)
        x3p = self.toPointsX(x3); y3p = self.toPointsY(y3)
        f = self.fileHandle
        f.write("<path d=\"M %.2f %.2f " % (x0p, y0p) )
        f.write("C %.2f %.2f %.2f %.2f %.2f %.2f\" " \
                % (x1p, y1p, x2p, y2p, x3p, y3p) )
        f.write("style=\"%s\" />\n" % self.getLineStyle(dashed=dashed) )
        return
    
    def text(self, x, y, textString, angle=0.0, fontSize=10,
             anchor="start", colour="black",
             fontFamily="sanserif"):
        """
        Render the textString at point (x,y).
        @param x: x-coordinate of anchor in user-space
        @type x: float or int
        @param y: y-coordinate of anchor in user-space
        @type y: float or int
        @param textString: string of characters to render
        @param angle: angle (in degrees) of text line wrt x-axis
                      (counterclockwise is positive)
        @type angle: float or int
        @param fontSize: size of font in points
        @type fontSize: int or float
        @param anchor: one of 'start', 'middle' or 'end'
        @type anchor: string
        @param colour: of the text
        @type colour: string
        """
        xp = self.toPointsX(x); yp = self.toPointsY(y)
        f = self.fileHandle
        f.write("<text x=\"%.2f\" y=\"%.2f\"" % (xp, yp) )
        if angle != 0.0:
            # my angle is positive counterclockwise.
            f.write(" transform=\"rotate(%.2f,%.2f,%.2f)\" " % (-angle, xp, yp) )
        f.write(" style=\"font-size:%g;text-anchor:%s;fill:%s;stroke:none" \
                % (fontSize,anchor,colour) )
        f.write(";font-weight:normal;font-family:%s\" " % fontFamily)
        f.write(">") # end of opening tag
        f.write(textString)
        f.write("</text>\n")
        return

    def dotlabel(self, x, y, label=None, anchor="middle", dotSize=2.0,
                 fontSize=10, colour="black",
                 fontFamily="sanserif"):
        """
        Render a dot with a text label.

        @param x: x-coordinate in user-space
        @param y: y-coordinate in user-space
        @param label: label text
        @param anchor: anchor location on label
        @param dotSize: dot diameter in mm
        @param textSize: font size in points
        @param colour: of both the label and the dot
        """
        xp = self.toPointsX(x); yp = self.toPointsY(y)
        rp = dotSize/2.0 * 72.0/25.4
        f = self.fileHandle
        f.write("<circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" " \
                % (xp, yp, rp) )
        f.write("style=\"fill:%s\" />\n" % colour )
        if label:
            # put the label slightly above the dot.
            f.write("<text x=\"%.2f\" y=\"%.2f\"" % (xp, yp-1.5*rp) )
            f.write(" style=\"font-size:%g;text-anchor:%s;fill:%s;stroke:none" \
                    % (fontSize,anchor,colour) )
            f.write(";font-weight:normal;font-family:%s\" " % fontFamily)
            f.write(">") # end of opening tag
            f.write(label)
            f.write("</text>\n")
        return
    
#-----------------------------------------------------------------------

if __name__ == "__main__":
    print "Begin test of svg_render..."
    s = SvgEnvironment(width=100, title="test run")
    s.open("test.svg")
    s.line(0.0, 0.0, 100.0, 120.0)
    s.close()

    s.open("test2.svg")
    s.line(0.0, 0.0, 90.0, 120.0)
    s.setLineWidth(0.5)
    s.circle(25.0, 85.0, 12.3)
    s.setLineWidth(0.25)
    s.polyline([0.0, 10.0, 20.0, 30.0], [50.0, 60.0, 50.0, 60.0], dashed=1)
    s.text(25.0, 85.0, "Circle", -30.0, anchor="middle")
    s.arc(90.0, 0.0, 60.0, 30.0, 60.0, 0.0)
    s.setLineWidth(0.75)
    s.arc(30.0, 0.0, 60.0, 30.0, 60.0, 0.0, dashed=1)
    s.bezier3(80.0, 80.0, 40.0, 80.0, 80.0, 100.0, 40.0, 100.0)
    s.setLineWidth(0.01)
    s.dotlabel(70.0,20.0,"a")
    s.close()
    print "Done."
    
