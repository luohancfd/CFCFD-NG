"""
Apply a matrix transformation such as rotation or translation to a Python function.

The functions provided by this module are used to manipulate a python function
prior to using the function to create a path with libprep3's PyFunctionPath.

Peter Blyton
    March 2011: Rotate and translate functions created.
"""

import math

def rotate_pyfunc(original_pyfunc, axis, angle):
    """
    Rotate a python function.
    
    Return a function that is the "original_pyfunc" function rotated about an
    axis "axis" through an angle of "angle".
    
    Arguments:
    original_pyfunc: (function) The name of the python function to be rotated.
        This python function must return a tuple of the 3 co-ordinates.
    axis: (str) The axis of rotation ["x", "y" or "z"].
    angle: (float) Counter-clockwise (right hand rule) angle of rotation [radians].
    
    Return Value:
    (function): Rotated version of the original python function.
    """
    def rotated_pyfunc_x(t):
        """
        original_pyfunc rotated about the x axis.
        
        This function is not called directly, it will be returned when
        calling rotate_pyfunc. The 3D x axis rotation matrix is used:
        
        |x'|   |1    0     0| |x|
        |y'| = |0 cos@ -sin@| |y|
        |z'|   |0 sin@  cos@| |z|
        """
        original_tuple = original_pyfunc(t)
        x = original_tuple[0]
        y = original_tuple[1]
        z = original_tuple[2]
        
        y_rotated = y*math.cos(angle) - z*math.sin(angle)
        z_rotated = y*math.sin(angle) + z*math.cos(angle)
        
        return (x, y_rotated, z_rotated)
    
    def rotated_pyfunc_y(t):
        """
        original_pyfunc rotated about the y axis.
        
        This function is not called directly, it will be returned when
        calling rotate_pyfunc. The 3D y axis rotation matrix is used:
        
        |x'|   | cos@ 0 sin@| |x|
        |y'| = |    0 1    0| |y|
        |z'|   |-sin@ 0 cos@| |z|
        """
        original_tuple = original_pyfunc(t)
        x = original_tuple[0]
        y = original_tuple[1]
        z = original_tuple[2]
        
        x_rotated = x*math.cos(angle) + z*math.sin(angle)
        z_rotated = -x*math.sin(angle) + z*math.cos(angle)
        
        return (x_rotated, y, z_rotated)
    
    def rotated_pyfunc_z(t):
        """
        original_pyfunc rotated about the z axis.
        
        This function is not called directly, it will be returned when
        calling rotate_pyfunc. The 3D z axis rotation matrix is used:
        
        |x'|   |cos@ -sin@ 0| |x|
        |y'| = |sin@  cos@ 0| |y|
        |z'|   |   0     0 1| |z|
        """
        original_tuple = original_pyfunc(t)
        x = original_tuple[0]
        y = original_tuple[1]
        z = original_tuple[2]
        
        x_rotated = x*math.cos(angle) - y*math.sin(angle)
        y_rotated = x*math.sin(angle) + y*math.cos(angle)
        
        return (x_rotated, y_rotated, z)
    
    if axis == "x":
        return rotated_pyfunc_x
    elif axis == "y":
        return rotated_pyfunc_y
    elif axis == "z":
        return rotated_pyfunc_z

def translate_pyfunc(original_pyfunc, new_origin):
    """
    Translate a python function.
    
    Return a function that is the "original_pyfunc" function with it's origin
    tranlated to "new_origin".
    
    Arguments:
    original_pyfunc: (function) The name of the python function to be rotated.
        This python function must return a tuple of the 3 co-ordinates.
    new_origin: (tuple) The new origin of original_pyfunc.
    
    Return Value:
    (function): Translated version of the original python function.
    """
    def translated_pyfunc(t):
        """
        Translated version of original_pyfunc.
        
        This function is not called directly, it will be returned when
        calling translate_pyfunc.
        """
        original_tuple = original_pyfunc(t)
        x = original_tuple[0]
        y = original_tuple[1]
        z = original_tuple[2]
        
        x_translated = x + new_origin[0]
        y_translated = y + new_origin[1]
        z_translated = z + new_origin[2]
        
        return (x_translated, y_translated, z_translated)
        
    return translated_pyfunc

