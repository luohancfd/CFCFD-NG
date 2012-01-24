-- very-viscous-air.lua
--
-- User-defined gas model adapted from Rowan's example
-- PJ, 08-Jun-2011

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

local mu0 = 1.0e1
local Pr = 1.0
local k0 = mu0 * C_p / Pr

-- Local helper functions
local sqrt, pow = math.sqrt, math.pow
local function sound_speed(gamma, R, T)
   return sqrt(gamma*R*T)
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

function eval_transport_coefficients(Q)
   -- Assume that all pertinent values in Q are
   -- at the correct state.  In this particular
   -- model, viscosity and thermal conductivity
   -- are constants.
   Q.mu = mu0
   Q.k[1] = k0
   return Q
end

function molecular_weight(isp)
   -- PJ added July 2010
   return R0/R
end

function eval_diffusion_coeficients(Q)
   -- PJ added July 2010
   Q.D_AB[1][1] = 0.0
   return Q
end
