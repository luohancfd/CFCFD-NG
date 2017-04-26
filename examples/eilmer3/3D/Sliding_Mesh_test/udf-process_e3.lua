-- udf-process.lua
-- This file sets up functions that will be called
-- from the main time-stepping loop.

-- Version 3 is based on vtx locations and implements a conservative scheme. I.e. overalpping areas are considered to set correct fluxes.

print("Hello from the set-up stage of udf-process.")
print("nblks=", nblks)

Area_up = 0

function at_timestep_start(args)


   -- Rport the rotational angle at regular intervals. 
   if (args.step % 50) == 0 then
      --OMEGA = 0.
      OMEGA = -60.e3/60 * 2 * math.pi
      print("Angle = ", args.t*OMEGA)
   end



   -- The followign code does some pre-allocation work for the sliding mesh interfce.       
   if args.step == 0 then

      -- load DATA  
      DATA,err = table.load("CONFIG_tbl.lua")
      N_BLADE = DATA[0]
      OMEGA = DATA[1]
      UP_row_list = DATA[2]
      DOWN_row_list = DATA[3]
      LIST = DATA[4]
      theta_min = DATA[11]
      theta_max = DATA[12]

      -- check that blks are in current block list. 
      -- This ensures code is only executed by correct MPI process
      flag = 0
      for i=0,(nblks-1) do
         if blks[i].id == 0 then
            flag = 1
            break
         end
      end

    -- Must only be run by MPI process that has interface.
      if flag == 1 then

          -- check that all interface blocks are present on current Node
          A = {}
          for i=0,(nblks-1) do
              table.insert(A,blks[i].id)
          end
          for k,b in pairs(LIST) do
            for i=0,(nblks-1) do
                if blks[i].id == b then
                   break
                end
                if i == nblks-1 then
                    print('ERROR: Not all blocks that are part of the sliding interface appear on same MPI Node.')
                    print('Blocks in interface: ', LIST)
                    print('Blocks on current Node: ', A)
                end
            end
          end

       

          BLK_MAP = {}
          -- Create Mapping table that allows mapping from global index to local indices.
          for i=0,(nblks-1) do
             for k,ind in pairs(LIST) do
                if blks[i].id == ind then
                   BLK_MAP[ind] = i
                   print("Mapped Global Block ", ind, " to local:",i) 
                else
                   print("Global Block:",blks[i].id, "Not mapped")       
                end
             end
          end

          -- upstream face
          -- Example 1:
          -- Looking Radially outwards
          --    Col1  Col0 
          -- +------+------+
          -- | blk1 | blk0 | ^
          -- |      |      | | Z,k
          -- +------+------+ |
          --            <----+
          --              T,j
          -- Example 2:
          -- Looking Radially outwards
          --   Col4  Col3  Col2  Col1  Col0 
          -- +-----+-----+-----+-----+-----+
          -- |     |     |     |     |     |
          -- |  12 |  13 |  14 |  15 |  16 |  Row0
          -- |     |     |     |     |     |  ^
          -- +-----+-----+-----+-----+-----+  | Z,k
          --                           T <----+
          --                               j


          print("\n")
          print("Getting Z-positions of vtx on upstream face")
          Z_up_pos_low = {}
          Z_up_pos_high = {}
          Z_up_k = {}
          Z_up_row = {}
          Z_up_area = {}
          indx = 1
          row = -1
          for k, ROW in pairs(UP_row_list) do -- get Z position
             global_id = ROW[1]  
             blk_id = BLK_MAP[global_id]
             print("Block currently being explored. blk_id =",blk_id)
             row = row + 1 
             imin = blks[blk_id].imin; imax = blks[blk_id].imax
             jmin = blks[blk_id].jmin; jmax = blks[blk_id].jmax
             kmin = blks[blk_id].kmin; kmax = blks[blk_id].kmax 
             -- print(imin,imax,jmin,jamx,kmin,kmax) 
             i = imin -- Summing up across WEST face
             j = jmin -- only need to follow one line
             for k=kmin,kmax do  
                face = sample_j_face(blk_id, i, j, k) 
                -- vtx3 = sample_vtx(blk_id, i, j+1, k)
                vtx0 = sample_vtx(blk_id, i, j, k)
                vtx4 = sample_vtx(blk_id, i, j, k+1)
                -- vtx7 = sample_vtx(blk_id, i, j+1, k+1)
                Z_up_pos_low[indx] = vtx0.z
                Z_up_pos_high[indx] = vtx4.z
                Z_up_k[indx] = k
                Z_up_row[indx] = row
                Z_up_area[indx] = face.area
                print("check creation of Z_up list: pos_low, pos_high, k, row:",Z_up_pos_low[indx],Z_up_pos_high[indx],Z_up_k[indx],Z_up_row[indx])
                indx = indx + 1
             end

          end 
          Z_up_length = indx-1 -- subtract +1 from last loop
          --store number of items in list for later access  

          print("\n")
          print("Getting Theta-positions of vtx on upstream face")
          T_up_pos_low={}
          T_up_pos_high={}
          T_up_j={}
          T_up_col={}
          T_up_area={} 
          indx = 1
          col = -1
          List = UP_row_list[1]
          for i=1,2 do -- get Z position
             blk_id = List[i]
             print("Block currently being explored. blk_id =",blk_id)
             col = col+1 
             imin = blks[blk_id].imin; imax = blks[blk_id].imax
             jmin = blks[blk_id].jmin; jmax = blks[blk_id].jmax
             kmin = blks[blk_id].kmin; kmax = blks[blk_id].kmax
             i = imin -- Summing up across WEST face
             k = kmin
             for j=jmin,jmax do  
                face = sample_j_face(blk_id, i, j, k) 
                vtx3 = sample_vtx(blk_id, i, j+1, k)
                vtx0 = sample_vtx(blk_id, i, j, k)
                -- vtx4 = sample_vtx(blk_id, i, j, k+1)
                -- vtx7 = sample_vtx(blk_id, i, j+1, k+1)
                T_up_pos_low[indx] = math.atan(vtx0.y/vtx0.x)
                T_up_pos_high[indx] = math.atan(vtx3.y/vtx3.x)
                T_up_j[indx] = j
                T_up_col[indx] = col
                T_up_area[indx] = face.area
                print("check creation of T_up list: angle_low, angle_high, j, col:",T_up_pos_low[indx],T_up_pos_high[indx],T_up_j[indx],T_up_col[indx])
                indx = indx + 1  
             end 
          end
          T_up_length = indx - 1 -- subtract +1 from last loop


          -- downstream face
          -- Example 1:
          -- Looking Radially outwards
          --      Col0 
          -- +-------------+
          -- |    blk3     |  
          -- |             |  
          -- +-------------+     
          -- |    blk2     |  ^
          -- |             |  | Z,k
          -- +-------------+  |   
          --            <-----+ 
          --              T,i
          --
          -- Example 2:
          -- Looking Radially outwards
          --   Col1  Col0 
          -- +-----+-----+
          -- |     |     |
          -- | IT0 | IT1 |  Row2
          -- |  20 |  21 |  
          -- +-----+-----+
          -- |     |     |
          -- | IC0 | IC1 |  Row1
          -- |  23 |  24 |  
          -- +-----+-----+
          -- |     |     |
          -- | IB0 | IB1 |  Row0
          -- |  26 |  27 |  ^
          -- +-----+-----+  | Z,k
          --         T <----+
          --              i  

          print("\n")
          print("Getting Z-positions of vtx on downstream face")
          Z_down_pos_low = {}
          Z_down_pos_high = {}
          Z_down_k = {}
          Z_down_row = {}
          Z_down_area = {}
          indx = 1
          row = - 1
          for k, ROW in pairs(DOWN_row_list) do -- get Z position
             global_id = ROW[1] 
             blk_id = BLK_MAP[global_id]
             print("Block currently being explored. blk_id =",blk_id)
             row = row + 1 
             imin = blks[blk_id].imin; imax = blks[blk_id].imax
             jmin = blks[blk_id].jmin; jmax = blks[blk_id].jmax
             kmin = blks[blk_id].kmin; kmax = blks[blk_id].kmax  
             i = imin -- Summing up across WEST face
             j = jmin -- only need to follow one line
             for k=kmin,kmax do  
                face = sample_j_face(blk_id, i, j, k) 
                vtx0 = sample_vtx(blk_id, i, j, k)
                -- vtx1 = sample_vtx(blk_id, i+1, j, k)
                -- vtx5 = sample_vtx(blk_id, i+1, j, k+1)
                vtx4 = sample_vtx(blk_id, i, j, k+1)
                Z_down_pos_low[indx] = vtx0.z
                Z_down_pos_high[indx] = vtx4.z
                Z_down_k[indx] = k
                Z_down_row[indx] = row
                Z_down_area[indx] = face.area
                print("check creation of Zdown list: pos_low, pos_high, k, row:",Z_down_pos_low[indx],Z_down_pos_high[indx],Z_down_k[indx],Z_down_row[indx])
                indx = indx + 1
             end
          end 
          Z_down_length = indx-1 -- subtract +1 from last loop
          --store nuber of items in list for later access  

          print("\n")
          print("Getting Theta-positions of vtx on downstream face")
          T_down_pos_low = {}
          T_down_pos_high = {}
          T_down_i = {}
          T_down_col = {}
          T_down_area = {}
          indx = 1
          col = -1
          List = DOWN_row_list[1]
          for i=1,1 do -- get T position
             blk_id = List[i]
             print("Block currently being explored. blk_id =",blk_id)
             col = col+1 
             imin = blks[blk_id].imin; imax = blks[blk_id].imax
             jmin = blks[blk_id].jmin; jmax = blks[blk_id].jmax
             kmin = blks[blk_id].kmin; kmax = blks[blk_id].kmax
             i = imin -- Summing up across WEST face
             k = kmin
             for i=imin,imax do  
                face = sample_j_face(blk_id, i, j, k) 
                vtx0 = sample_vtx(blk_id, i, j, k)
                vtx1 = sample_vtx(blk_id, i+1, j, k)
                -- vtx5 = sample_vtx(blk_id, i+1, j, k+1)
                -- vtx4 = sample_vtx(blk_id, i, j, k+1)
                T_down_pos_low[indx] = math.atan(vtx0.y/vtx0.x)
                T_down_pos_high[indx] = math.atan(vtx1.y/vtx1.x)
                T_down_i[indx] = i
                T_down_col[indx] = col
                T_down_area[indx] = face.area
                print("check creation of Tdown list: ang_low, ang_high, j, col:",T_down_pos_low[indx],T_down_pos_high[indx],T_down_i[indx],T_down_col[indx])
                indx = indx + 1    
             end 
          end
          T_down_length = indx-1 -- subtract +1 from last loop
          --store number of items in list for later access 


          print("\n")
          print("UP UP UP")
          -- create map for z direction (for mapping DOWN --> UP)
          UP_Nfaces = {}
          UP_face_ind = {}
          UP_row_ind = {}
          UP_weight = {}
          UP_area = {}
          for i=1,Z_up_length do
             U_low_vtx = Z_up_pos_low[i]
             U_high_vtx = Z_up_pos_high[i]
             count = 0
             faces = {}
             rows = {}
             weights = {}
             areas = {}
             for j = 1,Z_down_length do 
                D_low_vtx = Z_down_pos_low[j]    
                D_high_vtx = Z_down_pos_high[j]
                if ( (D_high_vtx <= U_high_vtx) and (D_low_vtx >= U_low_vtx) ) then -- face fully inside
                    count = count+1
                    faces[count] = Z_down_k[j]
                    rows[count] = Z_down_row[j]
                    weights[count] = 1.0
                    areas[count] = Z_down_area[j]
                    --print("fully inside")
                elseif ( (D_high_vtx > U_high_vtx) and (D_low_vtx < U_high_vtx) ) then -- only overlap at the top
                    count = count+1
                    faces[count] = Z_down_k[j]
                    rows[count] = Z_down_row[j]
                    weights[count] = (U_high_vtx-D_low_vtx) / (D_high_vtx-D_low_vtx)
                    areas[count] = Z_down_area[j]
                    --print("top overlap")
                elseif ( (D_high_vtx > U_low_vtx) and (D_low_vtx < U_low_vtx) ) then -- only overlap at the bottom
                    count = count+1
                    faces[count] = Z_down_k[j]
                    rows[count] = Z_down_row[j]
                    weights[count] = (D_high_vtx - U_low_vtx) / (D_high_vtx-D_low_vtx)
                    areas[count] = Z_down_area[j]
                    --print("bottom overlap")
                elseif ( (D_high_vtx > U_high_vtx) and (D_low_vtx < U_low_vtx) ) then -- overlap at both side
                    count = count+1
                    faces[count] = Z_down_k[j]
                    rows[count] = Z_down_row[j]
                    weights[count] = (U_high_vtx - U_low_vtx) / (D_high_vtx-D_low_vtx)
                    areas[count] = Z_down_area[j]
                    --print("both sides")
                else -- other cases
                    -- No overlap so no action required. 
                end
             end
             UP_Nfaces[i] = count
             UP_face_ind[i] = faces
             UP_row_ind[i] = rows
             UP_weight[i] = weights 
             UP_area[i] = areas        
          end

          for i=1,Z_up_length do
             print("Number of downstream faces overlaping with upstream face ",  i," is :",UP_Nfaces[i])
             for j=1,UP_Nfaces[i] do
                print("Face_ind, Row_ind, weight:", UP_face_ind[i][j],UP_row_ind[i][j],UP_weight[i][j])
             end
          end

          UP = {}
          UP[0] = UP_Nfaces
          UP[1] = UP_face_ind
          UP[2] = UP_row_ind
          UP[3] = UP_weight
          UP[4] = UP_area
          UP[5] = T_down_pos_low  
          UP[6] = T_down_pos_high
          UP[7] = T_down_col
          UP[8] = T_down_i
          UP[9] = T_down_length

          --print(T_down_i[0],T_down_i[1],T_down_i[2],T_down_i[3],T_down_i[4],T_down_i[5])
          
          print("\n")
          print("DOWN DOWN DOWN")
          -- create map for z direction (for mapping UP --> DOWN)
          DOWN_Nfaces = {}
          DOWN_face_ind = {}
          DOWN_row_ind = {}
          DOWN_weight = {}
          DOWN_area = {}
          for i=1,Z_down_length do
             D_low_vtx = Z_down_pos_low[i]
             D_high_vtx = Z_down_pos_high[i]
             count = 0
             faces = {}
             rows = {}
             weights = {}
             areas = {} 
             for j = 1,Z_up_length do 
                U_low_vtx = Z_up_pos_low[j]    
                U_high_vtx = Z_up_pos_high[j]
                if ( (U_high_vtx <= D_high_vtx) and (U_low_vtx >= D_low_vtx) ) then -- face fully inside
                    count = count+1
                    faces[count] = Z_up_k[j]
                    rows[count] = Z_up_row[j]
                    weights[count] = 1.0
                    areas[count] = Z_up_area[j]
                    --print("fully inside")
                elseif ( (U_high_vtx > D_high_vtx) and (U_low_vtx < D_high_vtx) ) then -- only overlap at the top
                    count = count+1
                    faces[count] = Z_up_k[j]
                    rows[count] = Z_up_row[j]
                    weights[count] = (D_high_vtx-U_low_vtx) / (U_high_vtx-U_low_vtx)
                    areas[count] = Z_up_area[j]
                    --print("top overlap")
                elseif ( (U_high_vtx > D_low_vtx) and (U_low_vtx < D_low_vtx) ) then -- only overlap at the bottom
                    count = count+1
                    faces[count] = Z_up_k[j]
                    rows[count] = Z_up_row[j]
                    weights[count] = (U_high_vtx - D_low_vtx) / (U_high_vtx-U_low_vtx)
                    areas[count] = Z_up_area[j]
                    --print("bottom overlap")
                elseif ( (U_high_vtx > D_high_vtx) and (U_low_vtx < D_low_vtx) ) then -- overlap at both side
                    count = count+1
                    faces[count] = Z_up_k[j]
                    rows[count] = Z_up_row[j]
                    weights[count] = (D_high_vtx - D_low_vtx) / (U_high_vtx-U_low_vtx)
                    areas[count] = Z_up_area[j]
                    --print("both sides")
                else -- other cases
                    -- No overlap so no action required. 
                end
             end
             DOWN_Nfaces[i] = count
             DOWN_face_ind[i] = faces
             DOWN_row_ind[i] = rows
             DOWN_weight[i] = weights 
             DOWN_area[i] = areas        
          end

          for i=1,Z_down_length do
             print("Number of upstream faces overlaping with downstream face ",  i," is :",DOWN_Nfaces[i])
             for j=1,DOWN_Nfaces[i] do
                print("Face_ind, Row_ind, weight:", DOWN_face_ind[i][j],DOWN_row_ind[i][j],DOWN_weight[i][j])
             end
          end


          DOWN = {}
          DOWN[0] = DOWN_Nfaces
          DOWN[1] = DOWN_face_ind
          DOWN[2] = DOWN_row_ind
          DOWN[3] = DOWN_weight
          DOWN[4] = DOWN_area
          DOWN[5] = T_up_pos_low 
          DOWN[6] = T_up_pos_high
          DOWN[7] = T_up_col
          DOWN[8] = T_up_j
          DOWN[9] = T_up_length
          --print(T_up_col[1],T_up_col[2],T_up_col[3])
     
          -- write mapping tables to files  
          assert( table.save( UP, "UP_tbl.lua" ) == nil )
          -- write mapping tables to files  
          assert( table.save( DOWN, "DOWN_tbl.lua" ) == nil )
     
       end
   end

   return
end

function at_timestep_end(args)
   --if (args.step % 100) == 0 then
   --   -- Run every 100 timesteps
   --   print("At end of timestep ", args.step, " t=", args.t)
   --   sum_total_mass(args)
   --end
   --print("At end of timestep")


   return
end

function sum_total_mass(args)
   --print("running mass summing fucntion")
   mass = 0.0
   for ib=0,(nblks-1) do
      imin = blks[ib].imin; imax = blks[ib].imax
      jmin = blks[ib].jmin; jmax = blks[ib].jmax
      kmin = blks[ib].kmin; kmax = blks[ib].kmax
      --print(imin, imax,jmin,jmax,kmin,kmax)
      blk_id = blks[ib].id
      for k=kmin,kmax do
         for j=jmin,jmax do
            for i=imin,imax do
               cell = sample_flow(blk_id, i, j, k)
               --print( blk_id, i, j, k)
               -- We are only given p and T
               -- so need to compute density
               -- using gas model
               Q = create_empty_gas_table()
               Q.p = cell.p
               Q.T = cell.T
               for isp=0,(nsp-1) do Q.massf[isp] = cell.massf[isp] end
               print("Here?",Q.p,Q.T[0])
               eval_thermo_state_pT(Q)
               rho = Q.rho
               -- Now we can compute mass in cell using volume of cell
               mass = mass + rho*cell.vol
            end
         end
      end
   end
   print("At timestep, ",args.step,", Mass (kg) of gas in domain: ", mass)
   return
end

function source_vector(args, cell_data)
   -- args contains t
   -- cell_data table contains most else
   src = {}
   src.mass = 0.0
   src.momemtum_x = 0.0
   src.momentum_y = 0.0
   src.momentum_z = 0.0
   if cell_data.x > 0.2 and cell_data.x < 0.4 and 
      cell_data.y > 0.2 and cell_data.y < 0.3 then
      src.total_energy = 100.0e+6  -- J/m**3
      src.energies = {[0]=100.0e6}
   else
      src.total_energy = 0.0
      src.energies = {[0]=0.0}
   end
   src.species = {[0]=0.0,}
   return src
end




--// test save to file
--assert( table.save( t, "test_tbl.lua" ) == nil )
   
-- load table from file
--local t2,err = table.load( "test_tbl.lua" )

--assert( err == nil )




-- declare local variables
--// exportstring( string )
--// returns a "Lua" portable version of the string
local function exportstring( s )
   return string.format("%q", s)
end

--// The Save Function
   function table.save(  tbl,filename )
      local charS,charE = "   ","\n"
      local file,err = io.open( filename, "wb" )
      if err then return err end

      -- initiate variables for save procedure
      local tables,lookup = { tbl },{ [tbl] = 1 }
      file:write( "return {"..charE )

      for idx,t in ipairs( tables ) do
         file:write( "-- Table: {"..idx.."}"..charE )
         file:write( "{"..charE )
         local thandled = {}

         for i,v in ipairs( t ) do
            thandled[i] = true
            local stype = type( v )
            -- only handle value
            if stype == "table" then
               if not lookup[v] then
                  table.insert( tables, v )
                  lookup[v] = #tables
               end
               file:write( charS.."{"..lookup[v].."},"..charE )
            elseif stype == "string" then
               file:write(  charS..exportstring( v )..","..charE )
            elseif stype == "number" then
               file:write(  charS..tostring( v )..","..charE )
            end
         end

         for i,v in pairs( t ) do
            -- escape handled values
            if (not thandled[i]) then
            
               local str = ""
               local stype = type( i )
               -- handle index
               if stype == "table" then
                  if not lookup[i] then
                     table.insert( tables,i )
                     lookup[i] = #tables
                  end
                  str = charS.."[{"..lookup[i].."}]="
               elseif stype == "string" then
                  str = charS.."["..exportstring( i ).."]="
               elseif stype == "number" then
                  str = charS.."["..tostring( i ).."]="
               end
            
               if str ~= "" then
                  stype = type( v )
                  -- handle value
                  if stype == "table" then
                     if not lookup[v] then
                        table.insert( tables,v )
                        lookup[v] = #tables
                     end
                     file:write( str.."{"..lookup[v].."},"..charE )
                  elseif stype == "string" then
                     file:write( str..exportstring( v )..","..charE )
                  elseif stype == "number" then
                     file:write( str..tostring( v )..","..charE )
                  end
               end
            end
         end
         file:write( "},"..charE )
      end
      file:write( "}" )
      file:close()
   end
   
--// The Load Function
function table.load( sfile )
   local ftables,err = loadfile( sfile )
   if err then return _,err end
   local tables = ftables()
   for idx = 1,#tables do
      local tolinki = {}
      for i,v in pairs( tables[idx] ) do
         if type( v ) == "table" then
            tables[idx][i] = tables[v[1]]
         end
         if type( i ) == "table" and tables[i[1]] then
            table.insert( tolinki,{ i,tables[i[1]] } )
         end
      end
      -- link indices
      for _,v in ipairs( tolinki ) do
         tables[idx][v[2]],tables[idx][v[1]] =  tables[idx][v[1]],nil
      end
   end
   return tables[1]
end


