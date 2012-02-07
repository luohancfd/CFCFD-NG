Geometry module
===============

.. automodule:: cfpylib.geom

.. currentmodule:: cfpylib/geom

minimal_geometry
----------------

.. automodule:: cfpylib.geom.minimal_geometry

.. autoclass:: Vector
   :members: __init__, __str__, __repr__, __add__, __neg__, __pos__, __sub__, __mul__, __rmul__, __div__, sum, __abs__, unit

.. autofunction:: dot

.. autofunction:: cross

.. autofunction:: quad_properties

.. autofunction:: quad_centroid

.. autofunction:: tetrahedron_properties

.. autofunction:: wedge_properties

.. autofunction:: hexahedron_properties

.. autofunction:: hexahedron_volume


svg_render
----------

.. automodule:: cfpylib.geom.svg_render

.. autoclass:: SvgEnvironment
   :members: __init__, toPointsX, toPointsY, setLineWidth, setDashArray, getLineStyle, open, add_comment, close, line, polyline, arc, circle, bezier3, text, dotlabel


transform_pyfunc
----------------

.. automodule:: cfpylib.geom.transform_pyfunc

.. autofunction:: rotate_pyfunc

.. autofunction:: rotate_pyfunc

