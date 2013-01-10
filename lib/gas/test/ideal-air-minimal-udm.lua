-- Author: Rowan J. Gollan
-- Date: 08-Jul-2008
--
-- User-defined gas model
-- ----------------------
-- This is an example model which 
-- implements ideal gas behaviour
-- for a single component gas.
--
-- This is a minimal implementation:
-- numerical techniques are used to 
-- give the rest of the functionality.
--
-- Notes:
--   20-Nov-2012 : Updated to compute thermal conductivity
--                 from Prandtl number
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

local Pr = 0.72

-- Local helper functions
local sqrt, pow = math.sqrt, math.pow
local function sound_speed(gamma, R, T)
   return sqrt(gamma*R*T)
end

local function Sutherland_viscosity(T)
   return mu0 * pow(T/T0_v, 3/2) * (T0_v + S_v)/(T + S_v)
end

local function thermal_conductivity(T)
   local mu = Sutherland_viscosity(T)
   local k = C_p*mu/Pr
   return 
end

-- Mandatory function:
function eval_thermo_state_rhoe(Q)
   -- Assume rho and e[1] are given, compute the
   -- remaining thermodynamic variables.
   -- Remember: we need to access the temperature
   -- and energy as the first value in an array
   -- of possible energies/temperatures.
   Q.T[0] = Q.e[0]/C_v
   Q.p = Q.rho*R*Q.T[0]
   Q.a = sound_speed(gamma, R, Q.T[0])
   -- Pass back the updated table
   return Q
end

function eval_transport_coefficients(Q)
   -- Assume that all pertinent values in Q are
   -- at the correct state.  In this particular
   -- model, viscosity and thermal conductivity
   -- are only dependent on temperature, ie. Q.T[1]
   Q.mu = Sutherland_viscosity(Q.T[0])
   Q.k[0] = thermal_conductivity(Q.T[0])
   return Q
end

function molecular_weight(isp)
   -- PJ added July 2010
   return R0/R
end

function eval_diffusion_coeficients(Q)
   -- PJ added July 2010
   Q.D_AB[0][0] = 0.0
   return Q
end
