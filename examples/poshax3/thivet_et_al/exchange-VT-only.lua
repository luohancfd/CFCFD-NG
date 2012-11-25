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
   'O2 ~~ (O2, N2) : V-T',
   rt={'SSH-VT'}
}

mechanism{
   'N2 ~~ (O2, N2) : V-T',
   rt={'SSH-VT'}
}

