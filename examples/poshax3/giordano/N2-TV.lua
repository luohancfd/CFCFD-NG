scheme{
    update = "energy exchange ODE",
    temperature_limits = {
        lower = 20.0,
        upper = 100000.0
    },
    error_tolerance = 0.000001
}

ode_solver{
    step_routine = 'rkf',
    max_step_attempts = 4,
    max_increase_factor = 1.15,
    max_decrease_factor = 0.01,
    decrease_factor = 0.333
}

mechanism{
   'N2 ~~ N2 : V-T',
   rt={'Landau-Teller-cf', A=7.12e-9, B = 124.07, C = 0.0}
}

