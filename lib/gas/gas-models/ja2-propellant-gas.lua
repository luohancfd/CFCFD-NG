-- JA2 propellant gas thermodynamic model.
-- 
-- Adapted from Rowan's ideal gas model for air by PJ.
--
-- Ref: Ian A. Johnston
--      The Noble-Abel Equation of State: 
--      Thermodynamic Derivations for Ballistics Modelling.
--      DST0-TN-0670 Defence Science and Technology Organisation
--      Australian Government, Department of Defence.
-- (just look it up with Google to get a copy)
--
-- The transport properties are left with the values for ideal-air.

model = 'composite gas'
equation_of_state = 'Noble-Abel gas'
thermal_behaviour = 'constant specific heats'
mixing_rule = 'Wilke'
sound_speed = 'equilibrium'
ignore_mole_fraction = 1.0e-15
species = {'propellant', }

propellant = {}
propellant.M = {
  value = 0.02489,
  reference = "Johnston (2005) page 9 Table 1 has R=334 J/(kg.K)",
  description = "molecular mass",
  units = "kg/mol",
}
propellant.gamma = {
  value = 1.225,
  reference = "Johnston (2005) page 9 Table 1",
  description = "(ideal) ratio of specific heats at room temperature",
  units = "non-dimensional",
}
propellant.b = {
  value = 0.001,
   reference = "Johnston (2005) page 9 Table 1",
  description = "co-volume",
  units = "m**3/kg",
}
propellant.d = {
  value = 3.617e-10,
  reference = "Air value from Bird, Stewart and Lightfoot (2001), p. 864",
  description = "equivalent hard-sphere diameter, sigma from L-J parameters",
  units = "m",
}
propellant.e_zero = {
  value = 0,
  description = "reference energy",
  units = "J/kg",
}
propellant.q = {
  value = 0,
  description = "heat release",
  units = "J/kg",
}
propellant.viscosity = {
  parameters = {
    T_ref = 273,
    ref = "Air value from Table 1-2, White (2006)",
    S = 111,
    mu_ref = 1.716e-05,
  },
  model = "Sutherland",
}
propellant.thermal_conductivity = {
  parameters = {
    S = 194,
    ref = "Air value from Table 1-3, White (2006)",
    k_ref = 0.0241,
    T_ref = 273,
  },
  model = "Sutherland",
}
