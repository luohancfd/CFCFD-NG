# __init__.py
"""
Flow-solution manipulation funtions.
"""

import blockflow2d
try:
    import vtk_xml_writer
except ImportError:
    pass
try:
    import shock_layer_surface
except:
    pass
