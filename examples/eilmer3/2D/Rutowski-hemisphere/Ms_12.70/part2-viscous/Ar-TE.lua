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

-- ET exchange via the Appleton model

mechanism{
   'e- ~~ Ar : E-T',
   rt={'Appleton-Bray:TwoRangeNeutral',
        T_switch=10000.0,
	sigma_low_T={  3.9e-21, -5.51e-25, 5.95e-29},
	sigma_high_T={-3.5e-21,  7.75e-25, 0.0}
   }
}

mechanism{
   'e- ~~ Ar+ : E-T',
   rt={'Appleton-Bray:Ion'}
}

