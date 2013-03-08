-- Author: Rowan J. Gollan
-- Date: 30-Oct-2008

-- diatomic iodine

I2 = {}
I2.M = { 
   value = 253.8089400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
I2.gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}

I2.atomic_constituents = {I=2}
I2.charge = 0
-- Presently these values are oxygen values.
I2.d = { 
   value = 3.433e-10,
   units = 'm',
   description = 'equivalent hard sphere diameter, based on L-J parameters',
   reference = 'Bird, Stewart and Lightfoot (2001), p. 864'
}
I2.viscosity = { 
   model = "Sutherland",
   parameters = { 
      mu_ref = 1.919e-05, T_ref = 273.0, S = 139.0,
      ref = "Table 1-2, White (2006)"
   }
}
I2.thermal_conductivity = { 
   model = "Sutherland",
   parameters = { 
      k_ref = 0.0244, T_ref = 273.0, S = 240.0,
      ref = "Table 1-3, White (2006)"
   }
}
