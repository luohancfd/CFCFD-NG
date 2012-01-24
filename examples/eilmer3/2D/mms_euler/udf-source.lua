-- Lua script for the source terms
-- of a Manufactured Solution which
-- treats Euler flow.
--
-- Author: Rowan J. Gollan
-- Date: 04-Feb-2008

-- dummy functions to keep eilmer3 happy

function at_timestep_start(args) return nil end
function at_timestep_end(args) return nil end

M_PI = math.pi
cos = math.cos
sin = math.sin
pow = math.pow

L = 1.0
gam = 1.4

rho0 = 1.0
rhox = 0.15
rhoy = -0.1

uvel0 = 800.0
uvelx = 50.0
uvely = -30.0

vvel0 = 800.0
vvelx = -75.0
vvely = 40.0
wvel0 = 0.0

press0 = 1.0e5
pressx = 0.2e5
pressy = 0.5e5

function rho_source(x, y)
   f_m = (3*M_PI*uvelx*cos((3*M_PI*x)/(2.*L))*(rho0 + rhoy*cos((M_PI*y)/(2.*L)) +
rhox*sin((M_PI*x)/L)))/(2.*L) + (2*M_PI*vvely*cos((2*M_PI*y)/(3.*L))*(rho0 + rhoy
*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)/L)))/(3.*L) + (M_PI*rhox*cos((M_PI*x)/L)
*(uvel0 + uvely*cos((3*M_PI*y)/(5.*L)) + uvelx*sin((3*M_PI*x)/(2.*L))))/L - (M_PI
*rhoy*sin((M_PI*y)/(2.*L))*(vvel0 + vvelx*cos((M_PI*x)/(2.*L)) + vvely*sin((2*M_PI
*y)/(3.*L))))/(2.*L)
    return f_m
end

function xmom_source(x, y)
   f_x = (3*M_PI*uvelx*cos((3*M_PI*x)/(2.*L))*(rho0 + rhoy*cos((M_PI*y)/(2.*L)) +
rhox*sin((M_PI*x)/L))*(uvel0 + uvely*cos((3*M_PI*y)/(5.*L)) + uvelx*sin((3*M_PI*x)
/(2.*L))))/L + (2*M_PI*vvely*cos((2*M_PI*y)/(3.*L))*(rho0 + rhoy*cos((M_PI*y)/(2.
*L)) + rhox*sin((M_PI*x)/L))*(uvel0 + uvely*cos((3*M_PI*y)/(5.*L)) + uvelx*sin((3
*M_PI*x)/(2.*L))))/(3.*L) + (M_PI*rhox*cos((M_PI*x)/L)*pow(uvel0 + uvely*cos((3*
M_PI*y)/(5.*L)) + uvelx*sin((3*M_PI*x)/(2.*L)),2))/L - (2*M_PI*pressx*sin((2*M_PI
*x)/L))/L - (M_PI*rhoy*(uvel0 + uvely*cos((3*M_PI*y)/(5.*L)) + uvelx*sin((3*M_PI*x)
/(2.*L)))*sin((M_PI*y)/(2.*L))*(vvel0 + vvelx*cos((M_PI*x)/(2.*L)) + vvely*sin((2*
M_PI*y)/(3.*L))))/(2.*L) - (3*M_PI*uvely*(rho0 + rhoy*cos((M_PI*y)/(2.*L)) + rhox*
sin((M_PI*x)/L))*sin((3*M_PI*y)/(5.*L))*(vvel0 + vvelx*cos((M_PI*x)/(2.*L)) +
vvely*sin((2*M_PI*y)/(3.*L))))/(5.*L)
    return f_x
end

function ymom_source(x, y)
   f_y = (M_PI*pressy*cos((M_PI*y)/L))/L - (M_PI*vvelx*sin((M_PI*x)/(2.*L))*(rho0 +
rhoy*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)/L))*(uvel0 + uvely*cos((3*M_PI*y)/(5.
*L)) + uvelx*sin((3*M_PI*x)/(2.*L))))/(2.*L) + (3*M_PI*uvelx*cos((3*M_PI*x)/(2.*L))
*(rho0 + rhoy*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)/L))*(vvel0 + vvelx*cos((M_PI
*x)/(2.*L)) + vvely*sin((2*M_PI*y)/(3.*L))))/(2.*L) + (4*M_PI*vvely*cos((2*M_PI*y)/
(3.*L))*(rho0 + rhoy*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)/L))*(vvel0 + vvelx*
cos((M_PI*x)/(2.*L)) + vvely*sin((2*M_PI*y)/(3.*L))))/(3.*L) + (M_PI*rhox*cos((M_PI
*x)/L)*(uvel0 + uvely*cos((3*M_PI*y)/(5.*L)) + uvelx*sin((3*M_PI*x)/(2.*L)))*(vvel0
+ vvelx*cos((M_PI*x)/(2.*L)) + vvely*sin((2*M_PI*y)/(3.*L))))/L - (M_PI*rhoy*sin((
M_PI*y)/(2.*L))*pow(vvel0 + vvelx*cos((M_PI*x)/(2.*L)) + vvely*sin((2*M_PI*y)/(3.*L)
),2))/(2.*L)
    return f_y
end

function energy_source(x, y)
   f_e = (uvel0 + uvely*cos((3*M_PI*y)/(5.*L)) + uvelx*sin((3*M_PI*x)/(2.*L)))*((-2*
M_PI*pressx*sin((2*M_PI*x)/L))/L + (rho0 + rhoy*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI
*x)/L))*((-2*M_PI*pressx*sin((2*M_PI*x)/L))/((-1 + gam)*L*(rho0 + rhoy*cos((M_PI*y)/
(2.*L)) + rhox*sin((M_PI*x)/L))) + ((3*M_PI*uvelx*cos((3*M_PI*x)/(2.*L))*(uvel0 + 
uvely*cos((3*M_PI*y)/(5.*L)) + uvelx*sin((3*M_PI*x)/(2.*L))))/L - (M_PI*vvelx*sin((
M_PI*x)/(2.*L))*(vvel0 + vvelx*cos((M_PI*x)/(2.*L)) + vvely*sin((2*M_PI*y)/(3.*L))))
/L)/2. - (M_PI*rhox*cos((M_PI*x)/L)*(press0 + pressx*cos((2*M_PI*x)/L) + pressy*sin(
(M_PI*y)/L)))/((-1 + gam)*L*pow(rho0 + rhoy*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)
/L),2))) + (M_PI*rhox*cos((M_PI*x)/L)*((pow(wvel0,2) + pow(uvel0 + uvely*cos((3*M_PI
*y)/(5.*L)) + uvelx*sin((3*M_PI*x)/(2.*L)),2) + pow(vvel0 + vvelx*cos((M_PI*x)/(2.*L)
) + vvely*sin((2*M_PI*y)/(3.*L)),2))/2. + (press0 + pressx*cos((2*M_PI*x)/L) + pressy
*sin((M_PI*y)/L))/((-1 + gam)*(rho0 + rhoy*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)/L
)))))/L) + (3*M_PI*uvelx*cos((3*M_PI*x)/(2.*L))*(press0 + pressx*cos((2*M_PI*x)/L) +
pressy*sin((M_PI*y)/L) + (rho0 + rhoy*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)/L))*
((pow(wvel0,2) + pow(uvel0 + uvely*cos((3*M_PI*y)/(5.*L)) + uvelx*sin((3*M_PI*x)/(2.*
L)),2) + pow(vvel0 + vvelx*cos((M_PI*x)/(2.*L)) + vvely*sin((2*M_PI*y)/(3.*L)),2))/2.
+ (press0 + pressx*cos((2*M_PI*x)/L) + pressy*sin((M_PI*y)/L))/((-1 + gam)*(rho0 + 
rhoy*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)/L))))))/(2.*L) + (2*M_PI*vvely*cos((2*
M_PI*y)/(3.*L))*(press0 + pressx*cos((2*M_PI*x)/L) + pressy*sin((M_PI*y)/L) + (rho0 +
rhoy*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)/L))*((pow(wvel0,2) + pow(uvel0 + uvely*
cos((3*M_PI*y)/(5.*L)) + uvelx*sin((3*M_PI*x)/(2.*L)),2) + pow(vvel0 + vvelx*cos((
M_PI*x)/(2.*L)) + vvely*sin((2*M_PI*y)/(3.*L)),2))/2. + (press0 + pressx*cos((2*M_PI*
x)/L) + pressy*sin((M_PI*y)/L))/((-1 + gam)*(rho0 + rhoy*cos((M_PI*y)/(2.*L)) + rhox*
sin((M_PI*x)/L))))))/(3.*L) + (vvel0 + vvelx*cos((M_PI*x)/(2.*L)) + vvely*sin((2*M_PI*
y)/(3.*L)))*((M_PI*pressy*cos((M_PI*y)/L))/L - (M_PI*rhoy*sin((M_PI*y)/(2.*L))*((pow(
wvel0,2) + pow(uvel0 + uvely*cos((3*M_PI*y)/(5.*L)) + uvelx*sin((3*M_PI*x)/(2.*L)),2)
+ pow(vvel0 + vvelx*cos((M_PI*x)/(2.*L)) + vvely*sin((2*M_PI*y)/(3.*L)),2))/2. +
(press0 + pressx*cos((2*M_PI*x)/L) + pressy*sin((M_PI*y)/L))/((-1 + gam)*(rho0 + rhoy
*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)/L)))))/(2.*L) + (rho0 + rhoy*cos((M_PI*y)/
(2.*L)) + rhox*sin((M_PI*x)/L))*((M_PI*pressy*cos((M_PI*y)/L))/((-1 + gam)*L*(rho0 +
rhoy*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)/L))) + ((-6*M_PI*uvely*(uvel0 + uvely*
cos((3*M_PI*y)/(5.*L)) + uvelx*sin((3*M_PI*x)/(2.*L)))*sin((3*M_PI*y)/(5.*L)))/(5.*L)
+ (4*M_PI*vvely*cos((2*M_PI*y)/(3.*L))*(vvel0 + vvelx*cos((M_PI*x)/(2.*L)) + vvely*
sin((2*M_PI*y)/(3.*L))))/(3.*L))/2. + (M_PI*rhoy*sin((M_PI*y)/(2.*L))*(press0 + 
pressx*cos((2*M_PI*x)/L) + pressy*sin((M_PI*y)/L)))/(2.*(-1 + gam)*L*pow(rho0 + rhoy
*cos((M_PI*y)/(2.*L)) + rhox*sin((M_PI*x)/L),2))))
    return f_e
end

function source_vector(args, cell)
   src = {}
   src.mass = rho_source(cell.x, cell.y)
   src.momentum_x = xmom_source(cell.x, cell.y)
   src.momentum_y = ymom_source(cell.x, cell.y)
   src.momentum_z = 0.0
   src.total_energy = energy_source(cell.x, cell.y)
   src.species = {}
   src.species[0] = src.mass
   return src
end
