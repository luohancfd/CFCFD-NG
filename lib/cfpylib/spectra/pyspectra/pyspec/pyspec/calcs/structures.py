import crystal

def Nickelate(x = 0):
    c = crystal.Crystal(3.822647, 3.822647, 12.72239)
    c.setSpaceGroup('I4/mmm')
    Sr = x / 2.
    La = (2 - x) / 2.
    c.addAtom('Ni', 'Ni', '2+', 0, 0, 0)
    c.addAtom('LaSr', 'La', '', 0, 0, 0.36180, occ = La)
    c.addAtom('LaSr', 'Sr', '', 0, 0, 0.36180, occ = Sr)
    c.addAtom('O1', 'O', '2-', 0, 0, 0.17372)
    c.addAtom('O2', 'O', '2-', 0.5, 0, 0)

    return c

def SrTiO3():
    c = crystal.Crystal(3.9034, 3.9304, 3.9304)
    c.setSpaceGroup('Pm3-m')
    c.addAtom('Ti1', 'Ti', '4+', 0.5, 0.5, 0.5)
    c.addAtom('Sr1', 'Sr', '', 0., 0., 0.)
    c.addAtom('O1', 'O', '2-', 0., 0.5, 0.5)
    
    return c

def BilayerManganite(x = 0, la = 0, sr = 0, thermal = False):
    """Return crystal object for a bilayer manganite"""

    if (la is 0) and (sr is 0):
        La = 2 - (2 * x)
        Sr = 1 + (2 * x)
    else:
        La = la
        Sr = sr

    print la, sr

    if thermal:
        U1 = 0.0021
        U2 = 0.0051
        U3 = 0.0098
        U4 = 0.0077
    else:
        U1 = 0.
        U2 = 0.
        U3 = 0.
        U4 = 0.

    c = crystal.Crystal(3.852, 3.852, 19.76)
    c.setSpaceGroup('I4/mmm')
    c.addAtom('Mn1', 'Mn', '3+', 0, 0, 0.09701, U = U1)
    c.addAtom('LaSr1', 'Sr', '', 0, 0, 0.5, occ = sr, U = U2)
    c.addAtom('LaSr1', 'La', '', 0, 0, 0.5, occ = la, U = U2)
    c.addAtom('LaSr2', 'La', '', 0, 0, 0.318198, occ = la, U = U2)
    c.addAtom('LaSr2', 'Sr', '', 0, 0, 0.318198, occ = sr, U = U2)
    
    c.addAtom('O1', 'O', '', 0, 0, 0, U = U3)
    c.addAtom('O2', 'O', '', 0, 0, 0.1951613, U = U3)
    c.addAtom('O3', 'O', '', 0.5, 0, 0.094909, U = U4)

    return c
