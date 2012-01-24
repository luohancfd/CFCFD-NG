## \file cluster_points.py
## \ingroup geom
## \brief Top-level clustering function that can use a range of
##        distribution schemes.
##
## \author P.Jacobs and A.Window
## \version 29-Nov-2005
"""
Clustering that can use a range of distribution schemes.

The intention is to cover the differences in PJ's and Adriaan's
clustering functions and thus make the grid generation code
a little easier to read.
"""
#----------------------------------------------------------------------

from roberts import distribute_points_1
import hyptan
import valliammai

#----------------------------------------------------------------------

def distribute_parameter_values(npoints, cluster_tuple=(0,0,0.0),
                                path_length=1.0):
    """
    Produce an array of parameter values distributed in the range [0.0,1.0].

    @param npoints: number of points, including the end points
    @type npoints: int
    @param cluster_tuple: the tuple of clustering parameters
        If there are three elements in the tuple, assume Robert's clustering.
        If there are two elements, use Adriaan's hyperbolic-tangent clustering.
        Default to uniform distribution otherwise.
    @type tuple
    @param path_length: Dimensional length of the path along which the parameter
        values are to be distributed.  This is needed by Adriaan's hyperbolic tangent
        clustering because the parameters specify the cell size at each end.
    @type path_length: float
    """
    try:
        if len(cluster_tuple) == 3:
            end0, end1, beta = cluster_tuple
            values = distribute_points_1(0.0, 1.0, npoints-1, end0, end1, beta)
        elif len(cluster_tuple) == 2:
            s0, s1 = cluster_tuple
            values = valliammai.make_cluster(path_length, npoints, s0, s1, iters=10, midFrac=0.2)
        elif len(cluster_tuple) == 4:
            s0, s1, end0, end1 = cluster_tuple
            values = hyptan.make_cluster(path_length, npoints, s0, s1, end0, end1, \
                iters=10, midFrac=0.2)
        else:
            raise ValueError, "Invalid cluster_tuple"
    except:
        # Default to uniform spacing.
        end0, end1, beta = 0, 0, 0.0
        values = distribute_points_1(0.0, 1.0, npoints-1, end0, end1, beta)
    return values
