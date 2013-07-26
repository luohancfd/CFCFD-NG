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

-- First attempt using standard models and a guess at sigma_dash, the limiting cross-section

mechanism{
   'H2 ~~ ( H2, H, He ) : V-T',
   rt={'Millikan-White:HTC', HTCS = { type = "Park", sigma_dash = 3.0e-18 } }
}

mechanism{
   'e- ~~ H+ : E-T',
   rt={'Appleton-Bray:Ion'}
}

