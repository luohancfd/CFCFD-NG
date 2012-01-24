## \file imp.py
## \author PJ, 19-Mar-2009

job_title = "MNM Implosion Problem."
print job_title
gdata.dimensions = 2

# Use a fudged air model
gas_gamma = 5.0/3.0
select_gas_model(model='ideal gas', species=['air'])
change_ideal_gas_attribute('air', 'gamma', gas_gamma)

L = 1.0
radius = L/2
pL = 100.0e3  # low pressure is 1 atm
pH = 10.0*pL  # high pressure

def my_domain(r, s, t=0.0):
    """
    The overall domain is a square of side L.

    User-defined function for the parametric volume maps from
    parametric space to physical space.
    Note that a (Python) tuple of coordinates is returned.
    """
    global L
    return (L*r, L*s, 0.0)

N = 100   # MNM's guess
dL = L/N  # cell width

def my_gas(x, y, z):
    """
    There is a circular region of low-pressure gas embedded in
    a larger, square region of high-pressure gas.
    Only one quarter of the full problem is simulated.

    User-defined function for the initial gas state
    works in physical space.
    Note that this function returns a dictionary
    of flow properties.
    """
    global dL, radius, pL, pH
    r2 = radius*radius
    x0 = x - 0.5*dL; x1 = x0 + dL
    y0 = y - 0.5*dL; y1 = y0 + dL
    r00 = x0*x0 + y0*y0
    r10 = x1*x1 + y0*y0
    r11 = x1*x1 + y1*y1
    r01 = x0*x0 + y1*y1
    if r00 < r2 and r10 < r2 and r11 < r2 and r01 < r2:
        # Fill the lower-left corner with low-pressure gas.
        p = pL
    elif r00 >= r2 and r10 >= r2 and r11 >= r2 and r01 >= r2:
        # and the outer-part of the field with high-pressure gas.
        p = pH
    else:
        # The cell is cut by the circular boundary.
        # Subdivide the cell to work out how much is inside radius.
        fcount = 0
        ddL = dL/10
        for i in range(10):
            xx = x0 + (i+0.5)*ddL
            for j in range(10):
                yy = y0 + (j+0.5)*ddL
                if xx*xx + yy*yy < r2:
                    fcount += 1
        f = float(fcount)/100.0
        p = f*pL + (1.0-f)*pH
    # We use the FlowCondition object to conveniently set all of
    # the relevant properties. 
    return FlowCondition(p=p, u=0.0, v=0.0, T=300.0, add_to_list=0).to_dict()

# Define a single block for the domain.
Block2D(PyFunctionSurface(my_domain), nni=N, nnj=N, fill_condition=my_gas)

# We can set individual attributes of the global data object.
# These are often used to control the simulation process.
gdata.title = job_title
gdata.flux_calc = AUSMDV
sound_speed = sqrt(gas_gamma*287.1*300.0)
print "sound_speed=", sound_speed
gdata.max_time = 0.296*L/sound_speed  # to match Fig.2 of Macrossan et al.
gdata.max_step = 600
gdata.dt = dL/5.0/sound_speed  # probably safe
gdata.dt_plot = gdata.max_time/4  # want some intermediate plots
print "low density is", pL/(287.1*300.0), "kg/m**3"


