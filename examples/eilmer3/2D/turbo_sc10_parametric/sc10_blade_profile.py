"""
Standard Condition 10 blade profile and camber functions.

Peter Blyton
    March 2011: Blade profile defined using functional form.

"""

import math

def thickness(s):
    """
    Modified NACA0006 aerofoil thickness distribution.
    
    Standard NACA0006 aerofoil equation modified to give a zero thickness
    at the trailing edge. Return the full aerofoil thickness
    from top to bottom surface, not just centerline to top surface.
    Equation source: http://rpmturbo.com/testcases/sc10/index.html
    
    Arguments:
    s: (float) The distance along the chord of aerofoil [0 <= s <= 1].
    
    Return Value:
    (float) Full aerofoil thickness.
    
    """
    return 0.06*(2.969*s**0.5 - 1.26*s - 3.516*s**2 + 2.843*s**3 - 1.036*s**4)

def camber(s):
    """
    Standard Configuration 10 camber line.
    
    Equation for the circular arc for the camber line of the Standard
    Configuration 10 blade profile. This is the upper arc of a circle
    where camber(0) = camber(1) = 0. Return the y coordinate of the arc,
    and the angle that a tangent makes above the horizontal.
    Equation source: http://rpmturbo.com/testcases/sc10/index.html
    
    Arguments:
    s: (float) The distance along the chord of aerofoil [0 <= s <= 1].
    
    Return Value:
    (tuple(float, float)) The y coordinate and angle of a tangent line above
    the horizontal at the location "s".
    
    """
    # Camber line is equation of a circle.
    # (x - a)**2 + (y - b)**2 = r**2
    a = 0.5
    b = -2.475
    r = 2.525
    
    y = b + math.sqrt(r**2 - (s - a)**2)
    
    # First derivative of camber line.
    dy_ds = -(s - a)/(math.sqrt(r**2 - (s - a)**2))
    
    # Angle that tangent makes above horizontal
    phi = math.atan(dy_ds)
    
    return (y, phi)

def sc10_top(s):
    """
    Standard Configuration 10 upper surface.
    
    Returns a tuple coordinate along the surface to be used with PyFunctionPath.    
    Equation source: http://rpmturbo.com/testcases/sc10/index.html
    
    Arguments:
    s: (float) The distance along the chord of aerofoil [0 <= s <= 1].
    
    """
    camber_data = camber(s)

    x = s - 0.5*thickness(s)*math.sin(camber_data[1])
    y = camber_data[0] + 0.5*thickness(s)*math.cos(camber_data[1])
    
    return (x, y, 0.0)

def sc10_bottom(s):
    """
    Standard Configuration 10 lower surface.
    
    Returns a tuple coordinate along the surface to be used with PyFunctionPath.
    Equation source: http://rpmturbo.com/testcases/sc10/index.html
    
    Arguments:
    s: (float) The distance along the chord of aerofoil [0 <= s <= 1].
    
    """
    camber_data = camber(s)
    
    x = s + 0.5*thickness(s)*math.sin(camber_data[1])
    y = camber_data[0] - 0.5*thickness(s)*math.cos(camber_data[1])
    
    return (x, y, 0.0)

