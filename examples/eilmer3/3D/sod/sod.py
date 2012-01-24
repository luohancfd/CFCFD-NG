## \file sod.py
## \brief Test job-specification file for e3prep.py
## \author PJ, 08-Sep-2006 adapted from Tcl script to Python
##             11-Feb-2009 ported to Eilmer3 to demonstrate the use
##                         of user-supplied functions for geometry
##                         and flow conditions.
job_title = "One-dimensional shock tube with air driving air."
print job_title
gdata.dimensions = 3

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])

def tube_volume(r, s, t):
    """
    User-defined function for the parametric volume maps from
    parametric space to physical space.
    Note that a (Python) tuple of coordinates is returned.
    """
    # A simple hexahedron, one unit long in the i-direction.
    return (1.0*r, 0.1*s, 0.1*t)

def tube_gas(x, y, z):
    """
    User-defined function for the initial gas state
    works in physical space.
    Note that this function returns a dictionary of flow properties.
    """
    if x < 0.5:
        # Fill the left-half of the volume with high-pressure gas.
        p = 1.0e5; T = 348.4
    else:
        # and the right-half with low-pressure gas.
        p = 1.0e4; T = 278.8
    # We use the FlowCondition object to conveniently set all of
    # the relevant properties. 
    return FlowCondition(p=p, u=0.0, v=0.0, T=T, add_to_list=0).to_dict()

# Define a single block for the tube.
Block3D(PyFunctionVolume(tube_volume), 
        nni=100, nnj=2, nnk=2,
        fill_condition=tube_gas)

# We can set individual attributes of the global data object.
# These are often used to control the simulation process.
gdata.title = job_title
gdata.flux_calc = AUSMDV
gdata.max_time = 0.6e-3  # seconds
gdata.max_step = 600
gdata.dt = 1.0e-6

