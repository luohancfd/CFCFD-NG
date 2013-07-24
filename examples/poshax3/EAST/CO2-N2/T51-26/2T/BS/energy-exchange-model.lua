scheme_t = {
    update = "energy exchange ODE",
    temperature_limits = {
        lower = 20.0,
        upper = 100000.0
    },
    error_tolerance = 0.000001
}

ode_t = {
    step_routine = 'rkf',
    max_step_attempts = 4,
    max_increase_factor = 1.15,
    max_decrease_factor = 0.01,
    decrease_factor = 0.333
}

-- all VT exchange mechanisms identified by Park (1994)
-- all ET exchange mechanisms from Gnoffo (1989)
-- using Park (1994) a and b values for VT relaxation

mechanism{
   'N2 ~~ ( N, N2, CO, CO2 ) : V-T',
   rt={'Millikan-White:HTC', HTCS = { type = "Park", sigma_dash = 3.0e-17 } }
}

mechanism{
   'N2 ~~ O : V-T',
   rt={'Millikan-White:HTC', a = 72.400, b = 0.015, HTCS = { type = "Park", sigma_dash = 3.0e-17 } }
}

mechanism{
   'N2 ~~ C : V-T',
   rt={'Millikan-White:HTC', a = 72.400, b = 0.015, HTCS = { type = "Park", sigma_dash = 3.0e-17 } }
}

mechanism{
   'CO ~~ ( N2, CO, CO2 ) : V-T',
   rt={'Millikan-White:HTC', HTCS = { type = "Park", sigma_dash = 3.0e-18 } }
}

mechanism{
   'CO ~~ N : V-T',
   rt={'Millikan-White:HTC', a = 47.700, b = 0.050, HTCS = { type = "Park", sigma_dash = 3.0e-18 } }
}

mechanism{
   'CO ~~ O : V-T',
   rt={'Millikan-White:HTC', a = 47.700, b = 0.050, HTCS = { type = "Park", sigma_dash = 3.0e-18 } }
}

mechanism{
   'CO ~~ C : V-T',
   rt={'Millikan-White:HTC', a = 47.700, b = 0.050, HTCS = { type = "Park", sigma_dash = 3.0e-18 } }
}

mechanism{
   'CO2 ~~ ( N, O, C, N2, CO ) : V-T',
   rt={'Millikan-White:HTC', HTCS = { type = "Park", sigma_dash = 1.0e-16 } }
}

mechanism{
   'CO2 ~~ CO2 : V-T',
   rt={'Millikan-White:HTC', a = 36.5000, b = -0.0193, HTCS = { type = "Park", sigma_dash = 1.0e-16 } }
}

mechanism{
   'e- ~~ C+ : E-T',
   rt={'Appleton-Bray:Ion'}
}

mechanism{
   'e- ~~ N+ : E-T',
   rt={'Appleton-Bray:Ion'}
}

mechanism{
   'e- ~~ O+ : E-T',
   rt={'Appleton-Bray:Ion'}
}

mechanism{
   'e- ~~ C : E-T',
   rt={'Appleton-Bray:Neutral', sigma_coefficients = { 1.2e-20, 1.7e-24,  -2.0e-29 } }
}

mechanism{
   'e- ~~ N : E-T',
   rt={'Appleton-Bray:Neutral', sigma_coefficients = { 5.0e-20, 0.0, 0.0 } }
}

mechanism{
   'e- ~~ O : E-T',
   rt={'Appleton-Bray:Neutral', sigma_coefficients = { 1.2e-20, 1.7e-24,  -2.0e-29 } }
}

mechanism{
   'e- ~~ CO : E-T',
   rt={'Appleton-Bray:Neutral', sigma_coefficients = { 2.0e-20, 6.0e-24,  0.0 } }
}

mechanism{
   'e- ~~ N2 : E-T',
   rt={'Appleton-Bray:Neutral', sigma_coefficients = { 7.5e-20, 5.5e-24, -1.0e-28 } }
}

mechanism{
   'e- ~~ NO : E-T',
   rt={'Appleton-Bray:Neutral', sigma_coefficients = { 1.0e-19, 0.0,      0.0 } }
}

mechanism{
   'e- ~~ O2 : E-T',
   rt={'Appleton-Bray:Neutral', sigma_coefficients = { 2.0e-20, 6.0e-24,  0.0 } }
}

