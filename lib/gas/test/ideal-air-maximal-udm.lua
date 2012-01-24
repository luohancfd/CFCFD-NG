-- Author: Rowan J. Gollan
-- Date: 08-Jul-2008
--
-- User-defined gas model
-- ----------------------
-- This is an example model which 
-- implements ideal gas behaviour
-- for a single component gas.
--

-- Mandatory, set nsp and nmodes
model = 'user-defined'
nsp = 1
nmodes = 1

-- Local parameters for model
local R0 = 8.31451
local R = 287.1
local gamma = 1.4
local C_v = R / (gamma - 1)
local C_p = R + C_v

local mu0 = 1.716e-5
local T0_v = 273.0
local S_v = 111.0

local k0 = 0.0241
local T0_k = 273.0
local S_k = 194.0

-- Local helper functions
local sqrt, pow = math.sqrt, math.pow
local function sound_speed(gamma, R, T)
   return sqrt(gamma*R*T)
end

local function Sutherland_viscosity(T)
   return mu0 * pow(T/T0_v, 3/2) * (T0_v + S_v)/(T + S_v)
end

local function Sutherland_thermal_conductivity(T)
   return k0 * pow(T/T0_k, 3/2) * (T0_k + S_k)/(T + S_k)
end

-- Mandatory function:
function eval_thermo_state_rhoe(Q)
   -- Assume rho and e[1] are given, compute the
   -- remaining thermodynamic variables.
   -- Remember: we need to access the temperature
   -- and energy as the first value in an array
   -- of possible energies/temperatures.
   Q.T[1] = Q.e[1]/C_v
   Q.p = Q.rho*R*Q.T[1]
   Q.a = sound_speed(gamma, R, Q.T[1])
   -- Pass back the updated table
   return Q
end

-- Optional functions:
-- If not supplied, an iterative scheme based on
-- the rhoe() function will be used.  It is more
-- computationally efficient to provided these
-- functions when it is easy to do so.

function eval_thermo_state_pT(Q)
   -- Assume p and T are given, compute the rest.
   Q.e[1] = C_v*Q.T[1]
   Q.rho = Q.p/(R*Q.T[1])
   Q.a = sound_speed(gamma, R, Q.T[1])
   return Q
end

function eval_thermo_state_rhoT(Q)
   -- Assume rho and T[1] are given, compute the rest.
   Q.e[1] = C_v*Q.T[1]
   Q.p = Q.rho*R*Q.T[1]
   Q.a = sound_speed(gamma, R, Q.T[1])
   return Q
end

function eval_thermo_state_rhop(Q)
   -- Assume rho and p are given, compute the rest.
   Q.T[1] = Q.p/(R*Q.rho)
   Q.e[1] = C_v*Q.T[1]
   Q.a = sound_speed(gamma, R, Q.T[1])
   return Q
end

function dTdp_const_rho(Q)
   -- Assume T, rho are given
   return 1.0/(Q.rho*R)
end

function dTdrho_const_p(Q)
   -- Assume T, p are given
   return (-1.0*Q.p)/(R*Q.rho^2)
end

function dpdrho_const_T(Q)
   -- Assume T, rho are given
   return R*Q.T[1]
end

function dedT_const_v(Q)
   return C_v
end

function dhdT_const_p(Q)
   return C_p
end

function molecular_weight(isp)
   return R0/R
end

function eval_transport_coefficients(Q)
   -- Assume that all pertinent values in Q are
   -- at the correct state.  In this particular
   -- model, viscosity and thermal conductivity
   -- are only dependent on temperature, ie. Q.T[1]
   Q.mu = Sutherland_viscosity(Q.T[1])
   Q.k[1] = Sutherland_thermal_conductivity(Q.T[1])
   return Q
end

function eval_diffusion_coeficients(Q)
   -- PJ added July 2010
   Q.D_AB[1][1] = 0.0
   return Q
end
