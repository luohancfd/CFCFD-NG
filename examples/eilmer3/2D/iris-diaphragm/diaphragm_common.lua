-- diaphragm_common.lua
-- Lua script for modelling an iris-diaphragm.
--
-- It is used as a user-defined boundary condition 
-- associated with the AdjacentPlusUDFBC.
--
-- This file is only part of the user-defined boundary condition.
-- It is called up by each of the diaphragm_?.lua files that 
-- has set some control parameters.
--
-- Adapted from udf-slip-wall.lua, by Stefan Hess, 01-Jun-2009
-- Merged with code "diaphragm-test.lua" by Peter Jacobs
-- Gradual opening and hold time capabilities added 02-Mar-2011 
-- by David Gildfind, Fabian Zander, Peter Jacobs

function reflect_normal_velocity(ux, vy, cosX, cosY)
   un = ux * cosX + vy * cosY;     -- Normal velocity
   vt = -ux * cosY + vy * cosX;    -- Tangential velocity
   un = -un;                       -- Reflect normal component
   ux = un * cosX - vt * cosY;     -- Back to Cartesian coords
   vy = un * cosY + vt * cosX;
   return ux, vy
end

function ghost_cell(args)
   -- Function that returns the flow state for a ghost cell
   -- for use in the inviscid flux calculations.
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   -- cell0 is the cell in the block; cell1 and cell2 are the ghost cells.
   i = args.i; j = args.j; k = args.k
   cell0 = sample_flow(block_id, i, j, k)
   if args.which_boundary == NORTH then
      cell1 = sample_flow(block_id, i, j+1, k)   
      cell2 = sample_flow(block_id, i, j+2, k)
   elseif args.which_boundary == EAST then
      cell1 = sample_flow(block_id, i+1, j, k)   
      cell2 = sample_flow(block_id, i+2, j, k)
   elseif args.which_boundary == SOUTH then
      cell1 = sample_flow(block_id, i, j-1, k)   
      cell2 = sample_flow(block_id, i, j-2, k)
   elseif args.which_boundary == WEST then
      cell1 = sample_flow(block_id, i-1, j, k)   
      cell2 = sample_flow(block_id, i-2, j, k)
   end

   if is_burst == false then
      -- Diaphragm hasn't ruptured, so we can proceed to see if it is newly burst.
      if math.abs(cell0.p-cell1.p) >= p_burst then 
	 -- Pressure difference exceeds burst pressure to trigger opening.
	 t_trigger = args.t
	 is_burst = true
	 --print("Diaphragm burst at t=",args.t,", in block,i,j = ",block_id,i,j)
      end
      -- Initially maintain wall boundary.
      cell0.u, cell0.v = 0.0, 0.0
      return cell0, cell0
   else
      cell0.u, cell0.v = 0.0, 0.0 -- maintain wall boundary.
      t_chk1 = args.t - t_trigger; -- amount of time after trigger occurs.
      t_chk2 = args.t - t_trigger - t_hold; -- amount of time after trigger AND hold time has occurred.
      if t_chk1<t_hold then
	 -- Diaphragm is held closed for time t_hold.
	 return cell0, cell0
      elseif t_chk2 < dt_burst then
	 -- Diaphragm is still opening.
	 y_current = math.pow(((args.t - t_trigger - t_hold)/dt_burst*y_max^2),0.5);
	 -- print("t = ",args.t," s, r = ",y_current," m")
	 if args.y > y_current then
 	    -- We are above open section of diaphragm;
	    -- maintain the solid boundary by filling the ghost cells with reflected data.
	    return cell0, cell0
	 else
  	    -- Open up the boundary by not overwriting the exchange data.
	    return nil, nil
	 end
      else 			
	 -- Diaphragm is fully open; always use the exchange data.
	 return nil, nil
      end
   end
end

function zero_normal_velocity(ux, vy, cosX, cosY)
   -- Just the interesting bits from reflect_normal_velocity().
   vt = -ux * cosY + vy * cosX;    -- Tangential velocity
   ux = -vt * cosY;                -- Back to Cartesian coords
   vy =  vt * cosX;                -- just tangential component
   return ux, vy
end

function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   face = sample_flow(block_id, args.i, args.j, args.k);
   -- note these are face (interface) values (i.e. between two cells).
   face.u, face.v = zero_normal_velocity(face.u, face.v, args.csX, args.csY);
   if is_burst==true then
      t_chk1 = args.t - t_trigger;
      t_chk2 = args.t - t_trigger - t_hold;
      if t_chk1 < t_hold then
	 --is diaphragm held closed for time t_hold?
	 return face
      elseif t_chk2 < dt_burst then
	 -- is diaphragm still opening.
	 -- assume diaphragm opens at constant rate of increase of area.
	 y_current = math.pow(((args.t - t_trigger - t_hold)/dt_burst*y_max^2),0.5);
	 --print("t = ",args.t," s, r = ",y_current," m")
	 if args.y > y_current then
	    return face
	 else
	    return nil
	 end
      else
	 --diaphragm is open
	 return nil	
      end
   else
      --diaphragm is closed
      return face
   end
end





