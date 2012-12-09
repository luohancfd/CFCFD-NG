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

rates = {
    {
        mechanisms = {
            {
                type = 'ET_exchange',
                relaxation_time = {
                    type = 'ET_AppletonBray_Ion',
                    c_name = 'Ar_plus'
                }
            },
            {
                type = 'ET_exchange',
                relaxation_time = {
                    type = 'ET_AppletonBray_TwoRangeNeutral',
                    c_name = 'Ar',
                    LT_sigma_coefficients = {  3.9e-19, -5.51e-23, 5.95e-27 },
                    HT_sigma_coefficients = { -3.5e-19,  7.75e-23, 0.0 },
                    T_switch = 10000.0
                }
            }  
        }
    }
}

equilibriation_mechanisms = {}
