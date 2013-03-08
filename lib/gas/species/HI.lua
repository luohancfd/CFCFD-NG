-- Author: Rowan J. Gollan
-- Date: 30-Oct-2008

-- hydrogen iodine

HI = {}
HI.M = { 
   value = 127.9124100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
HI.gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
HI.atomic_constituents = {H=1, I=1}
HI.charge = 0

HI.e_zero = {
   value = 206070.7,
   units = 'J/kg',
   description = 'reference energy, from CEA (enthalpy of formation)'
}

-- Presently these values are oxygen values.
HI.d = { 
   value = 3.433e-10,
   units = 'm',
   description = 'equivalent hard sphere diameter, based on L-J parameters',
   reference = 'Bird, Stewart and Lightfoot (2001), p. 864'
}
HI.viscosity = { 
   model = "Sutherland",
   parameters = { 
      mu_ref = 1.919e-05, T_ref = 273.0, S = 139.0,
      ref = "Table 1-2, White (2006)"
   }
}
HI.thermal_conductivity = { 
   model = "Sutherland",
   parameters = { 
      k_ref = 0.0244, T_ref = 273.0, S = 240.0,
      ref = "Table 1-3, White (2006)"
   }
}
