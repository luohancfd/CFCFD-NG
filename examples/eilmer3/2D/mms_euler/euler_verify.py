# euler_verify.py
from math import sin, cos, pi

R_air = 287.1

class EulerManufacturedSolution:
    def __init__(self, L,
                 rho0, rhox, rhoy, a_rhox, a_rhoy,
                 press0, pressx, pressy, a_pressx, a_pressy,
                 uvel0, uvelx, uvely, a_uvelx, a_uvely,
                 vvel0, vvelx, vvely, a_vvelx, a_vvely):
        self.L = L
        self.rho0 = rho0
        self.rhox = rhox
        self.rhoy = rhoy
        self.a_rhox = a_rhox
        self.a_rhoy = a_rhoy
        self.press0 = press0
        self.pressx = pressx
        self.pressy = pressy
        self.a_pressx = a_pressx
        self.a_pressy = a_pressy
        self.uvel0 = uvel0
        self.uvelx = uvelx
        self.uvely = uvely
        self.a_uvelx = a_uvelx
        self.a_uvely = a_uvely
        self.vvel0 = vvel0
        self.vvelx = vvelx
        self.vvely = vvely
        self.a_vvelx = a_vvelx
        self.a_vvely = a_vvely
        return

    def calculate_rho(self, x, y):
        rho = self.rho0
        rho += self.rhox*sin((self.a_rhox*pi*x)/self.L)
        rho += self.rhoy*cos((self.a_rhoy*pi*y)/self.L)
        return rho
    
    def calculate_p(self, x, y):
        p = self.press0
        p += self.pressx*cos((self.a_pressx*pi*x)/self.L)
        p += self.pressy*sin((self.a_pressy*pi*y)/self.L)
        return p

    def calculate_u(self, x, y):
        u = self.uvel0
        u += self.uvelx*sin((self.a_uvelx*pi*x)/self.L)
        u += self.uvely*cos((self.a_uvely*pi*y)/self.L)
        return u

    def calculate_v(self, x, y):
        v = self.vvel0
        v += self.vvelx*cos((self.a_vvelx*pi*x)/self.L)
        v += self.vvely*sin((self.a_vvely*pi*y)/self.L)
        return v
    

        
        
