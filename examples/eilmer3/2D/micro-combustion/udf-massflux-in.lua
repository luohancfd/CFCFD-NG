-- udf-massflux-in.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedGhostCell BC.


function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.

--------------------------!!!Input parameters------------------------------------
   dt_plot = 1e-4 -- timestep for saving the ghost cell information, s
   mass = 1.122497365547*0.4 -- mass flow rate, kg/s/m2
   T0 = 300 --total temperature, K
   massf = {}
   for isp=0,(nsp-1) do
      massf[isp] = 0.000000e+00   
   end
   massf[3] = 2.201527e-01  --O2
   massf[10] = 5.518596e-02 --CH4
   massf[19] = 7.246613e-01 --N2
   -- mass fractions are indexed from 0 to nsp-1
---------------------------------------------------------------------------------

------------------for the very first timestep, set values in the ghost cell------
   if (args.t == 0) then
      cell0 = sample_flow(block_id, args.i, args.j, args.k)
      filename = "update-"..string.format("%04d", args.j)..".data"
      file = io.open(filename, "w")
      file:write(dt_plot,"\t",cell0.u,"\t",cell0.p,"\t",cell0.T[0],"\n")
      file:close()
   end

   if (args.t_step == 0 and args.j == 2) then
----- initialize the data records in each ghost cell along j axis
      step_previous = {}
      dt_save = {} 
      update_u = {}
      update_p = {}
      update_T = {}
   end

   if (args.t_step == 0) then
      step_previous[args.j] = args.t_step
      filename = "update-"..string.format("%04d", args.j)..".data"
      file = io.open(filename, "r")
      dt_save[args.j], update_u[args.j], update_p[args.j], update_T[args.j] \
      = file:read("*number","*number","*number","*number")
      file:close()
   end
---------------------------------------------------------------------------------

   Q = create_empty_gas_table()
   Q.p = update_p[args.j]
   Q.T[0] = update_T[args.j]
   for isp=0,(nsp-1) do
      Q.massf[isp] = massf[isp]    
   end   
   eval_thermo_state_pT(Q)
   eval_sound_speed(Q)
   a = Q.a
   Cp = eval_Cp(Q)
   gamma = eval_gamma(Q)
   R = eval_R(Q)
   rho = Q.rho
   u = update_u[args.j]
   p = update_p[args.j]
   M = math.abs(u/a)
   -- Sample the flow field from the inner cells near the boundary.
   cell1 = sample_flow(block_id, args.i, args.j, args.k)
   x1 = cell1.x
   u1 = cell1.u
   p1 = cell1.p
   cell2 = sample_flow(block_id, args.i+1, args.j, args.k)
   x2 = cell2.x
   u2 = cell2.u
   p2 = cell2.p
   cell3 = sample_flow(block_id, args.i+2, args.j, args.k)
   x3 = cell3.x
   u3 = cell3.u
   p3 = cell3.p
   cell4 = sample_flow(block_id, args.i+3, args.j, args.k)
   x4 = cell4.x
   u4 = cell4.u
   p4 = cell4.p


   if (args.t_step ~= step_previous[args.j]) then
      -- NSCBC Wave amplitude and LODI relations by T.J.Poinsot:
      dpdx = (-25/12*p+4*p1-3*p2+4/3*p3-1/4*p4)/(x2-x1)
      dudx = (-25/12*u+4*u1-3*u2+4/3*u3-1/4*u4)/(x2-x1)
      L1 = (u-a)*(dpdx-rho*a*dudx) -- sound wave at speed u-c
      L2 = (1-M)/(M+1/(gamma-1))*L1 -- entropy wave at speed u
      L5 = (M-1)*(M*(gamma-1)-1)/(M+1)/(M*(gamma-1)+1)*L1
      -- sound wave at speed u+c

      -- As total temperature and mass flow rate are specified,
      -- only continuity equation needs to be solved on the boundary:
      d1 = 1/a/a*(L2+0.5*(L1+L5))
      rho = rho-d1*args.dt
      step_previous[args.j] = args.t_step
   end

   update_u[args.j] = mass/rho  
   update_T[args.j] = T0-0.5/Cp*update_u[args.j]*update_u[args.j]
   update_p[args.j] = rho*R*update_T[args.j]


   -- update ghost cells
   ghost = {}
   ghost.T = {}  -- temperatures, K (as a table)
   ghost.T[0] = update_T[args.j]
   ghost.u = update_u[args.j] -- x-velocity, m/s
   ghost.v = 0.0 -- y-velocity, m/s
   ghost.w = 0.0 -- z-velocity, m/s
   ghost.p = update_p[args.j] -- pressure, Pa
   ghost.massf = massf -- mass fractions


------------------save ghost cell information every dt_plot----------------------
   if (args.t >= dt_save[args.j]) then
      dt_save[args.j] = dt_save[args.j] + dt_plot
      filename = "update-"..string.format("%04d", args.j)..".data"
      file = io.open(filename, "w")
      file:write(dt_save[args.j],"\t",update_u[args.j],"\t",update_p[args.j],\
      "\t",update_T[args.j],"\n")
      file:close()

      filename1 = "data-records.data"
      file = io.open(filename1, "a")
      file:write(args.t_step,"\t",args.t_level,"\t",args.dt,"\t",d1,"\t",rho,\
      "\t",ghost.u,"\t",ghost.p,"\t",ghost.T[0],"\n")
      file:close()
   end
---------------------------------------------------------------------------------


   return ghost, ghost
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   return sample_flow(block_id, args.i, args.j, args.k)
end


