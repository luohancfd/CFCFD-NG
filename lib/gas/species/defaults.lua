-- Author: Rowan J. Gollan
-- Date: 29-July-2008

default = {}
default.charge = 0
default.atomic_constituents = {}
default.e_zero = {
   value = 0.0,
   units = 'J/kg',
   description = 'reference energy'
}
default.q = { 
   value = 0.0,
   units = 'J/kg',
   description = 'heat release'
}
default.d = {
   value = 3.617e-10,
   units = 'm',
   description = 'equivalent hard-sphere diameter, sigma from L-J parameters',
   reference = 'value for air: Bird, Stewart and Lightfoot (2001), p. 864'
}
default.p_a = {
   value = 101325.0,
   units = 'Pa',
   description = 'reference pressure for air in simple gas model',
   reference = 'Ask Kan (Jason) Qin'
}
default.rho_a = {
   value = 10.0,
   units = 'kg/m^3',
   description = 'reference density for air in simple gas model',
   reference = 'Ask Kan (Jason) Qin'
}
default.k_s = {
   value = 1.0,
   units = '-',
   description = 'reference k_s for air in simple gas model',
   reference = 'Ask Kan (Jason) Qin'
}
-- transport properties for air (as default)
default.viscosity = { 
   model = "Sutherland",
   parameters = { 
      mu_ref = 1.716e-05, T_ref = 273.0, S = 111.0,
      ref = "Table 1-2, White (2006)"
   }
}

default.thermal_conductivity = { 
   model = "Sutherland",
   parameters = {
      k_ref = 0.0241, T_ref = 273.0, S = 194.0,
      ref = "Table 1-3, White (2006)"
   }
}
-- when Blottner viscosity is required, use defaults for N2.
-- default.viscosity = {
--    model = "Blottner",
--    parameters = { 
--       A_mu = 0.026814, B_mu = 0.317784, C_mu = -11.31555,
--       ref = "Table 4, ESA Radiation Test Case 8 (2014)"
--    }
-- }
-- these properties vary widely between species, using the default is probably not a good idea
-- these values result in effectively an ideal gas
default.T_c = {
   value = 0.0,
   units = 'K',
   description = 'critical temperature',
}
default.p_c = {
   value = 1.0,
   units = 'Pa',
   description = 'critical pressure',
}
