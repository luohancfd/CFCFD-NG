***************************
Introduction to CCD Package
***************************

This chapter gives a short overview about the steps which are provided
by CCD package and the brief ideas.

The CCD rountines are broken down into three main classes or
objects. Their function is logically defined by the steps which one
performs to analyze data, from the raw CCD images to the resulting
data. These steps and their corresponding classes are:

* Reading of images - ``FileProcessor``
* Processing images and gridding - ``ImageProcessor``
* Plotting of data - ``PlotGrid``

These classes and their functions are described below.

Reading and compiling images
============================

The class ``FileProcessor`` reads all the images for a given dataset
and does both darkimage subtraction and image normalization. The
``FileProcessor`` obtains a list of CCD image file names, either
supplied by the user, or from a scan object. Upon processing, each sub
image is summed and corrected for the dark image then normalized to
the monitor. The ``FileProcessor`` returns a 3D array, or stack of
images.

In order to intelligently handle dark images, the latest dark image in
the stack is stored and used to correct the current image being
processed. This allows the user to only take periodic dark images at a
lesser frequency of data to speed up data collection.

Processing Images
=================

The class ``ImageProcessor`` takes a ``FileProcessor`` object and
further processes the data. 

In this class, first the data is transformed
into the coordinate system chosen by the user. This coordinate system
can be chosen for convenience to be anything from the diffractometer
theta coordinates to the HKL reciprocal space coordinates. 
  
Secondly, the data can be gridded onto a regular :math:`N x 3` grid to
allow efficient processing. This grid can be obtained as a 3D array of
average intensities for each voxel along with the corresponding
standard error for each voxel.

Finally, 1D cuts and sums can be performed to obtain line cuts of the
data.

Plotting Images
===============

In order to easily visualize the data a number of plotter routines are
provided to aid in quickly graphically representing the data. 

Important note on indexing
==========================

