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

-- all VT exchange mechanisms identified by Park (1993) [except NO~(N,0)], but
-- using thivet curve fit model
-- all ET exchange mechanisms from Gnoffo (1989)

mechanism{
   'N2 ~~ ( N2, NO, O2 ) : V-T',
   rt={'SSH-VT'}
}

mechanism{
   'N2 ~~ O : V-T',
   rt={'Thivet-cf', B=70.3, C=-24.35}
}

mechanism{
  'N2 ~~ N : V-T',
   rt={'Thivet-cf', B=70.3, C=-24.35}
}

mechanism{
   'NO ~~ ( N2, NO, O2 ) : V-T',
   rt={'SSH-VT'}
}

mechanism{
   'O2 ~~ ( N2, NO, O2 ) : V-T',
   rt={'SSH-VT'}
}

mechanism{
   'O2 ~~ O : V-T',
   rt={'Thivet-cf', B=55.2, C=-26.95}
}

 -- NOTE: B and C assumed to be same as for O2 ~ O
mechanism{
   'O2 ~~ N : V-T',
   rt={'Thivet-cf', B=55.2, C=-26.95}
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
   'e- ~~ N : E-T',
   rt={'Appleton-Bray:Neutral', sigma_coefficients = { 5.0e-20, 0.0, 0.0 } }
}

mechanism{
   'e- ~~ O : E-T',
   rt={'Appleton-Bray:Neutral', sigma_coefficients = { 1.2e-20, 1.7e-24,  -2.0e-29 } }
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
