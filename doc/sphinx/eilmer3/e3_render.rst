e3_render
=========

.. automodule:: e3_render

Rendering in 2D Scalable Vector graphics (SVG)
----------------------------------------------

To get an SVG rendering set up, it is usually convenient to set all of the
length scales with a call to the window method and then set the parameters
for the x and y axes so that they have suitable ranges, tick marks and positions. 

.. autoclass:: e3_render.SketchEnvironment

.. automethod:: e3_render.SketchEnvironment.do_labels

.. automethod:: e3_render.SketchEnvironment.xaxis

.. automethod:: e3_render.SketchEnvironment.yaxis

.. automethod:: e3_render.SketchEnvironment.scales

.. automethod:: e3_render.SketchEnvironment.set_drawing_size

.. automethod:: e3_render.SketchEnvironment.set_length_scale

.. automethod:: e3_render.SketchEnvironment.origin

.. automethod:: e3_render.SketchEnvironment.window


Miscellaneous
-------------

.. autofunction:: e3_render.rad_to_degrees
