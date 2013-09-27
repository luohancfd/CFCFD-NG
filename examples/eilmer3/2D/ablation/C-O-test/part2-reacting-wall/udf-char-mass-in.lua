-- udf-massflux-in.lua
-- User-specified mass-flux at an inlet boundary, defined
-- by reactions at the gas-surface interface between the
-- flowfield and char layer in an ablating solid.
-- EJF, 17/9/2013

-- mass fractions to be provided as a table
-- We take these data from a preliminary run of the preparation script.
nsp = 9
ntm = 2
massf = {}
for isp=0, nsp-1 do
   massf[isp] = 0.0
end

-- Species:
--['C','O2','O','CO','C2','C3','C_plus','O_plus','e_minus']

N_A = 6.023e23		-- Avogadro's number (1/mol)
k_B = 1.381e-23		-- Boltzmann's constant (m^2 kg s^-1 K^-1)

-- species molecular mass (kg/mol)
M = {}
for isp=0, nsp-1 do
   M[isp] = 0.0
end
M[0] = 12.0107000e-3
M[2] = 15.9994e-3
M[3] = 28.010100e-3

-- species molecular mass
m = {}
for isp=0, nsp-1 do
   m[isp] = M[isp]/N_A
end

function reaction_velocity(rho,T)

   T_w = T[0]
   rho1 = rho

   -- reaction rates (numbers match index for produced species)
   k_f = {}
   for isp=0, nsp-1 do
      k_f[isp] = 0.0
   end
   --k_f[2] = 0.63*math.exp(-1160.0/T_w)   -- reaction term for CO formation
   k_f[3] = 1.0
   v_reac = {}
   for isp=0, nsp-1 do
      v_reac[isp] = 0.0
   end
   v_reac[3] = k_f[3] * math.sqrt(k_B*T_w/(2.0*math.pi*m[2]))    -- m_O = m[1]

   return v_reac
end

function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   -- print("Hello from function ghost_cell.")
   -- Sample the flow field at the current cell 
   -- which is beside the boundary.
   cell = sample_flow(block_id, args.i, args.j, args.k)
   ghost = {}
   ghost.p = cell.p -- pressure, Pa
   ghost.T = {} -- temperatures, K (as a table)
   ghost.T[0] = cell.T[0]
   ghost.T[1] = cell.T[1]
   ghost.u = -cell.u  -- x-velocity, m/s
   ghost.v = cell.v     -- y-velocity, m/s
   ghost.w = 0.0
   ghost.massf = {} -- mass fractions to be provided as a table
   for isp=0,nsp-1 do
      ghost.massf[isp] = cell.massf[isp]
   end
   return ghost, ghost
end

function interface(args)
--   print('hello from interface')
   face = {}
   face = sample_flow(block_id, args.i, args.j, args.k)
   return face
end

function convective_flux(args)
   -- Function that returns the fluxes of conserved quantities.
   -- For use in the inviscid flux calculations.
   --print("Hello from function flux.")
   cell = {}
   cell = sample_flow(block_id, args.i, args.j, args.k)
   u = -cell.u       -- x-velocity, m/s
   v = cell.v       -- y-velocity, m/s
   w = 0.0
   p = cell.p
   rho = cell.rho
   T = {}
   T[0] = cell.T[0]
   T[1] = cell.T[1]

   Q = create_empty_gas_table(nsp,ntm)   
   Q.T[0] = T[0]
   Q.T[1] = T[1]
   for isp=0,nsp-1 do 
      Q.massf[isp] = cell.massf[isp]
   end
   Cv = eval_Cv(Q)
   v_reac = reaction_velocity(rho,T)
   massf_t = 0.0
   for isp=0,nsp-1 do
      massf[isp] = cell.massf[isp]
   end
   -- Assemble flux vector
   F = {}
   -- species mass flux
   F.species = {}
   for isp=0,nsp-1 do
      F.species[isp] = 0.0
   end
   F.species[3] = -rho * massf[2] * v_reac[3]   -- created due to amount of O present
   F.species[2] = rho * massf[2] * v_reac[3] -- amount of O is depleted.
   F.mass = 0.0
   for isp=0,nsp-1 do
      F.mass = F.mass + F.species[isp]
   end
   F.momentum_x = p * args.csX + u * F.mass
   F.momentum_y = p * args.csY + v * F.mass
   F.momentum_z = 0.0
   F.renergies = {}
   F.renergies[0] = F.mass * (Cv*T[0])
   F.renergies[1] = F.mass * (Cv*T[1])
   F.total_energy = F.renergies[0] + F.renergies[1] + F.mass * (0.5*(u*u+v*v) + p/rho)
   return F
end
