-- udf-massflux-in.lua
-- User-specified mass-flux at an inlet boundary, defined
-- by reactions at the gas-surface interface between the
-- flowfield and char layer in an ablating solid.
-- EJF, 17/9/2013

-- mass fractions to be provided as a table
-- We take these data from a preliminary run of the preparation script.
nsp = 3
ntm = 1
massf = {}
for isp=0, nsp-1 do
   massf[isp] = 0.0
end
--massf[0] = 1.0

-- Order of species:
--  0: C     1: O     2: CO

--  0: C     1: O     2: N      3: H      4: CO 
--  5: C2    6: N2    7: CN     8: NO     9: O2
-- 10: H2   11: C3   12: C2H   13: C+    14: O+ 
-- 15: H+   16: N+   17: NO+   18: N2+   19: e-

N_A = 6.023e23		-- Avogadro's number (1/mol)
k_B = 1.381e-23		-- Boltzmann's constant (m^2 kg s^-1 K^-1)

-- species molecular mass (kg/mol)
M = {}
for isp=0, nsp-1 do
   M[isp] = 0.0
end
M[0] = 12.0107000e-3
M[1] = 15.9994e-3
M[2] = 28.010100e-3

-- species molecular mass
m = {}
for isp=0, nsp-1 do
   m[isp] = M[isp]/N_A
end

function reaction_velocity(rho)

   T_w = 3500.0	-- if it's pre-set! or take from flow props...
   rho1 = rho

   -- reaction rates (numbers match index for produced species)
   k_f = {}
   for isp=0, nsp-1 do
      k_f[isp] = 0.0
   end
   --k_f[2] = 0.63*math.exp(-1160.0/T_w)   -- reaction term for CO formation
   k_f[2] = 1.0
   v_reac = {}
   for isp=0, nsp-1 do
      v_reac[isp] = 0.0
   end
   v_reac[2] = k_f[2] * math.sqrt(k_B*T_w/(2.0*math.pi*m[1]))    -- m_O = m[1]

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
   ghost.T[0] = cell.T
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
   face.T[0] = 3500.0
   --face.T[1] = T_wall             -- just like in the python script, need to set T[0] and T[1] for the two-temp model.
   --face.massf = {}
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
   rhoxu = rho*(u*args.csX + v*args.csY)
   Cv = 100.0   		--FIX-ME
   T = 3500.0		--FIX-ME
   v_reac = reaction_velocity(rho)
   massf_t = 0.0
   for isp=0,nsp-1 do		-- massf created at top of lua file.
      massf[isp] = cell.massf[isp]
   end
   --print(massf[0],massf[1],massf[2])
   -- Assemble flux vector
   F = {}
   -- species mass flux
   F.species = {}
   for isp=0,nsp-1 do
      F.species[isp] = 0.0
   end
   F.species[2] = -rho * massf[1] * v_reac[2]   -- created due to amount of O present
   F.species[1] =  rho * massf[1] * v_reac[2] -- amount of O is depleted.
   --print(F.species[1],F.species[2])
   F.mass = 0.0
   for isp=0,nsp-1 do
      F.mass = F.mass + F.species[isp]
   end
   --F.mass = rho * (u*args.csX + v*args.csY) -- kg/s/m**2
   F.momentum_x = p * args.csX + u * F.mass
   F.momentum_y = p * args.csY + v * F.mass
   F.momentum_z = 0.0
   F.total_energy = F.mass * (Cv*T + 0.5*(u*u+v*v) + p/rho)
   F.renergies = {}
   F.renergies[0] = F.mass * (Cv*T)
   return F
end
