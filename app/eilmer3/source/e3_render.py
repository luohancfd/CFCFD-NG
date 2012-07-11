"""
e3_render.py -- Classes for rendering 2D and 3D blocks.

It is expected that this module be imported by the application program e3prep.py.

.. Versions:
   16-Mar-2008 : extracted from e3prep.py (formerly mbcns_prep.py)
"""

from libprep3 import *
from e3_defs import *
from cfpylib.geom.svg_render import *
from e3_block import *


def rad_to_degrees(rad):
    """
    Convert radians to degrees.
    
    This is a convenience function for the writing of the MetaPost file
    but may also be used in the user script.

    :param rad: angle in radians.
    :returns: angle in degrees
    """
    return rad * 180.0 / math.pi


class SketchEnvironment(object):
    """
    A place to put the sketch settings for both Metapost rendering
    and SVG rendering.

    A sketch file will contain a rendering of the geometry objects
    defining the simulation domain.
    Since the coordinates for mb_cns are in metres, you will probably
    have to apply suitable scale factors to get a drawing that can
    be printed to an A4 page conveniently.
    Also, because the origin of a postsript figure is in the bottom-left
    corner of the page, you may have to reset the origin to see all
    of the geometry elements if some of them have negative coordinates.
    """

    def __init__(self):
        self.root_file_name = ""
        # We might like to shift the plot origin so that the postscript
        # picture as produced my metapost is not off the bottom-left
        # of the plotted page.
        self.xorigin_plot = 0.0
        self.yorigin_plot = 0.0

        self.xmin = 0.0
        self.xmax = 1.0
        self.xtic = 0.2
        self.xaxis_offset = -0.05

        self.ymin = 0.0
        self.ymax = 1.0
        self.ytic = 0.2
        self.yaxis_offset = -0.04

        # Usually we will be plotting onto an A4 page, so we do not
        # want the plotted picture to be too big.
        self.xscale = 0.12
        self.yscale = 0.12

        # Text is usually desired so will be turned on by default
        self.withText = True
        return

    def do_labels(self, label_state=True):
        """
        Set the state of label printing: True or False
        
        :param label_state: State of label printing: True or False
        """
        assert type(label_state) == bool
        self.withText = label_state
        return

    def xaxis(self, xmin=None, xmax=None, xtic=None, xaxis_offset=None):
        """
        Set the x-axis scale parameters.

        :param xmin: Minimum value for x-axis scale. (float)
        :param xmax: Maximum value for x-axis scale. (float)
        :param xtic: Interval between tic marks and labels. (float)
        :param xaxis_offset: The vertical offset (from ymin) for
            drawing the length of the scale. (float)
            Negative values will lower the x-axis scale.
        """
        if xmin != None:
            self.xmin = xmin
        if xmax != None:
            self.xmax = xmax
        if xtic != None:
            self.xtic = xtic
        if xaxis_offset != None:
            self.xaxis_offset = xaxis_offset
        return

    def yaxis(self, ymin=None, ymax=None, ytic=None, yaxis_offset=None):
        """
        Set the y-axis scale parameters.

        :param ymin: Minimum value for y-axis scale.
        :param ymax: Maximum value for y-axis scale.
        :param ytic: Interval between tic marks and labels.
        :param yaxis_offset: The horizontal offset (from xmin) for
            drawing the length of the scale.
            Negative values will move the y-axis scale to the left.
        """
        if ymin != None:
            self.ymin = ymin
        if ymax != None:
            self.ymax = ymax
        if ytic != None:
            self.ytic = ytic
        if yaxis_offset != None:
            self.yaxis_offset = yaxis_offset
        return

    def scales(self, xscale=None, yscale=None):
        """
        Set the scale factors for the drawing.

        Model coordinates are multiplied by these scales to get page coordinates.
        """
        if xscale != None:
            self.xscale = xscale
        if yscale != None:
            self.yscale = yscale
        return

    def set_drawing_size(self, width=None, height=None):
        """
        Set the scale factors for the drawing based on a length specified by the user.

        Depending on which scale is specified, the function searches the list
        of defined nodes and finds the max distance from the origin and
        scales the drawing by setting that distance as that specified in this
        function by the user.

        :param width: The size of the drawing in the x-direction. (Units = metres)
        :param height: The size of the drawing in the y-direction. (Units = metres)

        Note that only one scale needs to be specified. 
        In the case that one scale is specifed, 
        the drawing will hold true shape and scale in both directions 
        by the same factor (as computed by the function, not as specified by 
        the user). If both are specified, the drawing will be distorted to 
        meet the users demands.
        """
        # search through nodeList for the maximum x and y coordinates on which
        # to scale the drawing
        maxx = 0.0
        maxy = 0.0
        minx = 0.0
        miny = 0.0
        for node in Node.nodeList:
            if node.x > maxx:
                maxx = node.x
            if node.y > maxy:
                maxy = node.y
            if node.x < minx:
                minx = node.x
            if node.y < miny:
                miny = node.y
        if width != None and height != None:
            # scale both x and y with their respective length scales
            assert abs(maxx) > 0.0
            assert abs(maxy) > 0.0
            self.xscale = width / abs(maxx - minx)
            self.yscale = height / abs(maxy - miny)
        elif height != None:
            # scale both x and y with the y-length-scale
            self.xscale = height / abs(maxy - miny)
            self.yscale = height / abs(maxy - miny)
        elif width != None:
            # scale both x and y with the x-length scale
            self.xscale = width / abs(maxx - minx)
            self.yscale = width / abs(maxx - minx)
        return

    def set_length_scale(self, length=None, ref_length=0.1):
        """
        Sets the scale of the drawing given a user specified length.

        Reference length is 10cm on the page by default, but can also be 
        changed by user.
        If no parameters are specified the function searches for the largest
        orthogonal distance and scales that to the reference length on the 
        page while maintaining true shape.

        :param length: User specified length. (Units = metres)
        :param ref_length: Reference length. (Units = metres)
        """
        if length != None:
            self.xscale = ref_length / length
            self.yscale = ref_length / length
        else:
            # search through nodeList for the maximum x and y coordinates on which
            # to scale the drawing
            maxx = 0.0
            maxy = 0.0
            minx = 0.0
            miny = 0.0
            for node in Node.nodeList:
                if node.x > maxx:
                    maxx = node.x
                if node.y > maxy:
                    maxy = node.y
                if node.x < minx:
                    minx = node.x
                if node.y < miny:
                    miny = node.y
            if abs(maxx - minx) > abs(maxy - miny):
                length = abs(maxx - minx)
            else:
                length = abs(maxy - miny)
            self.xscale = length / ref_length
            self.yscale = length / ref_length
        return

    def origin(self, x=0.0, y=0.0):
        """
        Set the origin on the page for the rendered picture.

        For example, it is sometimes good to select an origin
        of (0.05, 0.05) to get the origin 5 centimetres up and
        right from the bottom-left corner of the page.
        """
        self.xorigin_plot = x
        self.yorigin_plot = y
        return

    def window(self, xmin=0.0, ymin=0.0, xmax=1.0, ymax=1.0,
               page_xmin=0.05, page_ymin=0.05, page_xmax=0.17, page_ymax=0.17):
        """
        Define a world-view window and a corresponding window on the rendered page.

        These windows are defined by the coordinates of their lower-left and
        upper-right corners.
        Distances, in both views, are in metres.
        """
        self.xscale = (page_xmax - page_xmin) / (xmax - xmin)
        self.yscale = (page_ymax - page_ymin) / (ymax - ymin)
        self.xorigin_plot = page_xmin - self.xscale * xmin
        self.yorigin_plot = page_ymin - self.yscale * ymin
        return
        
    def transform(self, x, y):
        x_plot = x * self.xscale + self.xorigin_plot
        y_plot = y * self.yscale + self.yorigin_plot
        return (x_plot, y_plot)

    def render_xaxis(self, fp):
        max_tic = 100
        x0, y0 = self.transform(self.xmin, self.ymin+self.xaxis_offset)
        x1, y1 = self.transform(self.xmax, self.ymin+self.xaxis_offset)
        xmid = 0.5 * (x0 + x1)
        if isinstance(fp, SvgEnvironment):
            fp.setLineWidth(0.5)
            fp.line(x0, y0, x1, y1)
            fp.text(xmid, y0-0.009, "x")  # user-space is measured in metres
        else:
            # Assume file handle for the Metapost file
            fp.write("pickup pencircle scaled 0.5mm;\n")
            fp.write("draw (%fu,%fu)--(%fu,%fu);\n" % (x0, y0, x1, y1) )
            fp.write("label.bot(\"%s\", (%fu,%fu) shifted (0,-0.6cm));\n" %
                     ("x", xmid, y0) )
        count = 0
        x = self.xmin
        while x < self.xmax + 1.0e-6:
            xu, yu = self.transform(x, self.ymin+self.xaxis_offset)
            if isinstance(fp, SvgEnvironment):
                fp.line(xu, yu, xu, yu+0.002)
            else:
                # Assume file handle for the Metapost file
                fp.write("draw (%fu,%fu)--(%fu,%fu) shifted (0,0.2cm);\n" %
                         (xu, yu, xu, yu) )
            if math.fabs(x) < 1.0e-6:
                xstring = "0"
            else:
                xstring = "%g" % x
            if isinstance(fp, SvgEnvironment):
                fp.text(xu, yu-0.005, xstring, anchor="middle")
            else:
                # Assume file handle for the Metapost file
                fp.write("label.bot(\"%s\", (%fu,%fu) shifted (-0.1cm,-0.1cm));\n" %
                         (xstring, xu, yu) )
            count += 1
            if count > max_tic: break
            x += self.xtic
        return

    def render_yaxis(self, fp):
        max_tic = 100
        x0, y0 = self.transform(self.xmin+self.yaxis_offset, self.ymin)
        x1, y1 = self.transform(self.xmin+self.yaxis_offset, self.ymax)
        ymid = 0.5 * (y0 + y1)
        if isinstance(fp, SvgEnvironment):
            fp.setLineWidth(0.5)
            fp.line(x0, y0, x1, y1)
            fp.text(x0-0.008, ymid, "y", anchor="end")  # user-space is measured in metres
        else:
            # Assume file handle for the Metapost file
            fp.write("pickup pencircle scaled 0.5mm;\n")
            fp.write("draw (%fu,%fu)--(%fu,%fu);\n" % (x0, y0, x1, y1) )
            fp.write("label.lft(\"%s\", (%fu,%fu) shifted (-0.8cm,0));\n" %
                     ("y", x0, ymid) )
        count = 0
        y = self.ymin
        while y < self.ymax + 1.0e-6:
            xu, yu = self.transform(self.xmin+self.yaxis_offset, y)
            if isinstance(fp, SvgEnvironment):
                fp.line(xu, yu, xu+0.002, yu)
            else:
                # Assume file handle for the Metapost file
                fp.write("draw (%fu,%fu)--(%fu,%fu) shifted (0.2cm,0);\n" %
                         (xu, yu, xu, yu) )
            if math.fabs(y) < 1.0e-6:
                ystring = "0"
            else:
                ystring = "%g" % y
            if isinstance(fp, SvgEnvironment):
                fp.text(xu-0.003, yu, ystring, anchor="end")
            else:
                # Assume file handle for the Metapost file
                fp.write("label.lft(\"%s\", (%fu,%fu) shifted (-0.3cm,0));\n" %
                         (ystring, xu, yu) )
            count += 1
            if count > max_tic: break
            y += self.ytic
        return
    
    def render_face(self, block_object, fp, which_face):
        # Label for boundary condition should be along the path,
        # at the middle of the path.
        nni = block_object.nni
        nnj = block_object.nnj
        psurf = block_object.psurf
        grid = block_object.grid
        if which_face == NORTH:
            offset_suffix = ".top"
            if( psurf != None ):
                p1 = psurf.eval(0.45, 1.0)
                p2 = psurf.eval(0.55, 1.0)
            else:
                p1 = Vector3(grid.x[int(nni/2)-1][nnj][0], grid.y[int(nni/2)-1][nnj][0])
                p2 = Vector3(grid.x[int(nni/2)][nnj][0], grid.y[int(nni/2)][nnj][0])
        elif which_face == SOUTH:
            offset_suffix = ".top"
            if( psurf != None ):
                p1 = psurf.eval(0.45, 0.0)
                p2 = psurf.eval(0.55, 0.0)
            else:
                p1 = Vector3(grid.x[int(nni/2)-1][0][0], grid.y[int(nni/2)-1][0][0])
                p2 = Vector3(grid.x[int(nni/2)][0][0], grid.y[int(nni/2)][0][0])
        elif which_face == WEST:
            offset_suffix = ".lft"
            if( psurf != None ):
                p1 = psurf.eval(0.0, 0.45)
                p2 = psurf.eval(0.0, 0.55)
            else:
                p1 = Vector3(grid.x[0][int(nnj/2)-1][0], grid.y[0][int(nnj/2)-1][0])
                p2 = Vector3(grid.x[0][int(nnj/2)][0], grid.y[0][int(nnj/2)][0])
        elif which_face == EAST:
            offset_suffix = ".lft"
            if( psurf != None ):
                p1 = psurf.eval(1.0, 0.45)
                p2 = psurf.eval(1.0, 0.55)
            else:
                p1 = Vector3(grid.x[nni][int(nnj/2)-1][0], grid.y[nni][int(nnj/2)-1][0])
                p2 = Vector3(grid.x[nni][int(nnj/2)][0], grid.y[nni][int(nnj/2)][0])
        else:
            print "render_face(): Unknown face type."
        x_mid, y_mid = self.transform(0.5*(p1.x+p2.x), 0.5*(p1.y+p2.y))
        theta = rad_to_degrees(math.atan2((p2.y-p1.y), (p2.x-p1.x)))
        bc_label = bcName[block_object.bc_list[which_face].type_of_BC]
        if block_object.bc_list[which_face].type_of_BC == ADJACENT:
            bc_label = None
        if bc_label and self.withText:
            fp.write( ("label%s(\"%s\" infont defaultfont scaled " +
                       "defaultscale rotated %f, (%fu, %fu));\n") %
                      (offset_suffix, bc_label, theta, x_mid, y_mid) )
        # Render the curve as seen by the grid-generator..
        if which_face == NORTH:
            N = nni  
            dt = 1.0/N
            if ( psurf != None ):
                p = psurf.eval(0.0, 1.0)
            else:
                p = Vector3(grid.x[0][nnj][0], grid.y[0][nnj][0])
            x, y = self.transform(p.x, p.y)
            fp.write("draw (%fu, %fu)" % (x, y))
            for i in range(1,N+1):
                if( psurf != None ):
                    p = psurf.eval(i*dt, 1.0)
                else:
                    p = Vector3(grid.x[i][nnj][0], grid.y[i][nnj][0])
                x, y = self.transform(p.x, p.y)
                fp.write("--(%fu, %fu)" % (x, y))
        elif which_face == SOUTH:
            N = nni  
            dt = 1.0/N
            if ( psurf != None ):
                p = psurf.eval(0.0, 0.0)
            else:
                p = Vector3(grid.x[0][0][0], grid.y[0][0][0])
            x, y = self.transform(p.x, p.y)
            fp.write("draw (%fu, %fu)" % (x, y))
            for i in range(1,N+1):
                if( psurf != None ):
                    p = psurf.eval(i*dt, 0.0)
                else:
                    p = Vector3(grid.x[i][0][0], grid.y[i][0][0])
                x, y = self.transform(p.x, p.y)
                fp.write("--(%fu, %fu)" % (x, y))
        elif which_face == WEST:
            N = nnj  
            dt = 1.0/N
            if ( psurf != None ):
                p = psurf.eval(0.0, 0.0)
            else:
                p = Vector3(grid.x[0][0][0], grid.y[0][0][0])
            x, y = self.transform(p.x, p.y)
            fp.write("draw (%fu, %fu)" % (x, y))
            for j in range(1,N+1):
                if( psurf != None ):
                    p = psurf.eval(0.0, j*dt)
                else:
                    p = Vector3(grid.x[0][j][0], grid.y[0][j][0])
                x, y = self.transform(p.x, p.y)
                fp.write("--(%fu, %fu)" % (x, y))
        elif which_face == EAST:
            N = nnj  
            dt = 1.0/N
            if ( psurf != None ):
                p = psurf.eval(1.0, 0.0)
            else:
                p = Vector3(grid.x[nni][0][0], grid.y[nni][0][0])
            x, y = self.transform(p.x, p.y)
            fp.write("draw (%fu, %fu)" % (x, y))
            for j in range(1,N+1):
                if( psurf != None ):
                    p = psurf.eval(1.0, j*dt)
                else:
                    p = Vector3(grid.x[nni][j][0], grid.y[nni][j][0] )
                x, y = self.transform(p.x, p.y)
                fp.write("--(%fu, %fu)" % (x, y))
        if block_object.bc_list[which_face].type_of_BC == ADJACENT:
            line_style = "dashed evenly"
        else:
            line_style = "";
        fp.write("%s;\n" % line_style)
        return
    
    def svg_render_face(self, block_object, svg, which_face):
        # Label should be along the path, at the middle of the path.
        nni = block_object.nni
        nnj = block_object.nnj
        psurf = block_object.psurf
        grid = block_object.grid
        if which_face == NORTH:
            if( psurf != None ):
                p1 = psurf.eval(0.45, 1.0)
                p2 = psurf.eval(0.55, 1.0)
            else:
                p1 = Vector3(grid.x[int(nni/2)-1][nnj][0], grid.y[int(nni/2)-1][nnj][0])
                p2 = Vector3(grid.x[int(nni/2)][nnj][0], grid.y[int(nni/2)][nnj][0])
        elif which_face == SOUTH:
            if( psurf != None ):
                p1 = psurf.eval(0.45, 0.0)
                p2 = psurf.eval(0.55, 0.0)
            else:
                p1 = Vector3(grid.x[int(nni/2)-1][0][0], grid.y[int(nni/2)-1][0][0])
                p2 = Vector3(grid.x[int(nni/2)][0][0], grid.y[int(nni/2)][0][0])
        elif which_face == WEST:
            if( psurf != None ):
                p1 = psurf.eval(0.0, 0.45)
                p2 = psurf.eval(0.0, 0.55)
            else:
                p1 = Vector3(grid.x[0][int(nnj/2)-1][0], grid.y[0][int(nnj/2)-1][0])
                p2 = Vector3(grid.x[0][int(nnj/2)][0], grid.y[0][int(nnj/2)][0])
        elif which_face == EAST:
            if( psurf != None ):
                p1 = psurf.eval(1.0, 0.45)
                p2 = psurf.eval(1.0, 0.55)
            else:
                p1 = Vector3(grid.x[nni][int(nnj/2)-1][0], grid.y[nni][int(nnj/2)-1][0])
                p2 = Vector3(grid.x[nni][int(nnj/2)][0], grid.y[nni][int(nnj/2)][0])
        else:
            print "svg_render_face(): Unknown face type."
        x_mid, y_mid = self.transform(0.5*(p1.x+p2.x), 0.5*(p1.y+p2.y))
        theta = rad_to_degrees(math.atan2((p2.y-p1.y), (p2.x-p1.x)))
        bc_label = bcName[block_object.bc_list[which_face].type_of_BC]
        if block_object.bc_list[which_face].type_of_BC == ADJACENT:
            bc_label = None
        if bc_label and self.withText:
            svg.text(x_mid, y_mid, bc_label, theta, anchor="middle")
        # Render the curve as seen by the grid-generator..
        if which_face == NORTH:
            N = nni  
            dt = 1.0/N
            if( psurf != None ):
                p = psurf.eval(0.0, 1.0)
            else:
                p = Vector3(grid.x[0][nnj][0], grid.y[0][nnj][0])
            x, y = self.transform(p.x, p.y)
            xlist = [x,]; ylist = [y,]
            for i in range(1,N+1):
                if( psurf != None ):
                    p = psurf.eval(i*dt, 1.0)
                else:
                    p = Vector3(grid.x[i][nnj][0], grid.y[i][nnj][0])
                x, y = self.transform(p.x, p.y)
                xlist.append(x); ylist.append(y)
        elif which_face == SOUTH:
            N = nni  
            dt = 1.0/N
            if( psurf != None ):
                p = psurf.eval(0.0, 0.0)
            else:
                p = Vector3(grid.x[0][0][0], grid.y[0][0][0])
            x, y = self.transform(p.x, p.y)
            xlist = [x,]; ylist = [y,]
            for i in range(1,N+1):
                if( psurf != None ):
                    p = psurf.eval(i*dt, 0.0)
                else:
                    p = Vector3(grid.x[i][0][0], grid.y[i][0][0])
                x, y = self.transform(p.x, p.y)
                xlist.append(x); ylist.append(y)
        elif which_face == WEST:
            N = nnj  
            dt = 1.0/N
            if( psurf != None ):
                p = psurf.eval(0.0, 0.0)
            else:
                p = Vector3(grid.x[0][0][0], grid.y[0][0][0])
            x, y = self.transform(p.x, p.y)
            xlist = [x,]; ylist = [y,]
            for j in range(1,N+1):
                if( psurf != None ):
                    p = psurf.eval(0.0, j*dt)
                else:
                    p = Vector3(grid.x[0][j][0], grid.y[0][j][0])
                x, y = self.transform(p.x, p.y)
                xlist.append(x); ylist.append(y)
        elif which_face == EAST:
            N = nnj  
            dt = 1.0/N
            if( psurf != None ):
                p = psurf.eval(1.0, 0.0)
            else:
                p = Vector3(grid.x[nni][0][0], grid.y[nni][0][0])
            x, y = self.transform(p.x, p.y)
            xlist = [x,]; ylist = [y,]
            for j in range(1,N+1):
                if( psurf != None ):
                    p = psurf.eval(1.0, j*dt)
                else:
                    p = Vector3(grid.x[nni][j][0], grid.y[nni][j][0] )
                x, y = self.transform(p.x, p.y)
                xlist.append(x); ylist.append(y)
        svg.polyline(xlist, ylist,
                     dashed=(block_object.bc_list[which_face].type_of_BC == ADJACENT))
        return
    
    def write_file(self, g, flowList, blockList):
        print "Begin write MetaPost file."
        fp = open(self.root_file_name+".mpost", "w")
        fp.write("%% %s\n" % gdata.title)
        fp.write("% Metapost file generated by scriptit.py\n")
        fp.write("prologues:=2;\n")
        fp.write("beginfig(1);\n")
        fp.write("u = 100cm;\n")  # mb_cns works in metres
        #
        self.render_xaxis(fp)
        self.render_yaxis(fp)
        # Render nodes
        fp.write("%% There are %d nodes\n" % len(Node.nodeList))
        fp.write("pickup pencircle scaled 0.25mm;\n")
        for node in Node.nodeList:
            x, y = self.transform(node.x, node.y)
            if node.label and self.withText:
                fp.write("dotlabel.lrt(\"%s\", (%fu, %fu));\n" %
                         (node.label, x, y))
            else:
                fp.write("dotlabel.lrt(\"\", (%fu, %fu));\n" % (x, y))
        # Render the blocks
        fp.write("%% There are %d blocks\n" % len(blockList) )
        fp.write("pickup pencircle scaled 0.25mm;\n")
        for blk in blockList:
            fp.write("%% Render block %d\n" % blk.blkId)
            # Put the label a little above the middle of the block
            # to (hopefully) avoid label clashes with the bounding paths.
            if ( blk.psurf != None ):
                pos = blk.psurf.eval(0.5, 0.67)
            else:
                pos = Vector3(blk.grid.x[int(blk.nni/2)][int(2*blk.nnj/3)][0],
                              blk.grid.y[int(blk.nni/2)][int(2*blk.nnj/3)][0])
            x_mid, y_mid = self.transform(pos.x, pos.y)
            if blk.label and self.withText:
                fp.write("label(\"%s\", (%fu, %fu));\n" %
                         (blk.label, x_mid, y_mid))
            # Render the block boundary.
            for face in faceList:
                self.render_face(blk, fp, face)
        # Finish Metapost file
        fp.write("endfig;\n")
        fp.write("end\n");
        fp.close()
        print "End write Metapost file."
        return
    
    def write_svg_file(self, g, flowList, blockList, faceList):
        """
        Render the Block2D objects to a scalable-vector-graphics file.

        This SVG image is useful for documentation as well as debugging.
        It may be later edited with a suitable program, such as inkscape.
        """
        if g.dimensions != 2:
            print "Will not render a 2D SVG image when dimensions=", g.dimensions
            return
        svgName = self.root_file_name+".svg"
        print "Begin write SVG file:", svgName
        svg = SvgEnvironment(0.20, 0.20, unitLength="m",
                             title=g.title,
                             desc="SVG file written by e3prep.py")
        svg.open(svgName)
        #
        self.render_xaxis(svg)
        self.render_yaxis(svg)
        # Render nodes
        svg.add_comment("There are %d nodes\n" % len(Node.nodeList))
        for node in Node.nodeList:
            x, y = self.transform(node.x, node.y)
            if self.withText:
                svg.dotlabel(x, y, node.label, anchor="end")
            else:
                svg.dotlabel(x, y, None, anchor="end")
        # Render the blocks
        svg.add_comment("There are %d blocks\n" % len(blockList) )
        svg.setLineWidth(0.25)
        for blk in blockList:
            svg.add_comment("Render block %d\n" % blk.blkId)
            # Put the label a little above the middle of the block
            # to (hopefully) avoid label clashes with the bounding paths.
            if( blk.psurf != None ):
                pos = blk.psurf.eval(0.5, 0.67)
            else:
                pos = Vector3(blk.grid.x[int(blk.nni/2)][int(2*blk.nnj/3)][0],
                              blk.grid.y[int(blk.nni/2)][int(2*blk.nnj/3)][0])
            x_mid, y_mid = self.transform(pos.x, pos.y)
            if blk.label and self.withText:
                svg.text(x_mid, y_mid, blk.label, anchor="middle")
            # Render the block boundary
            for face in faceList:
                self.svg_render_face(blk, svg, face)
        svg.close()
        print "End write SVG file."
        return
