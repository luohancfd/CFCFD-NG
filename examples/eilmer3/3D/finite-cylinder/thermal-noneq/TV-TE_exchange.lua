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

-- all VT exchange mechanisms identified by Park (1993)
-- all ET exchange mechanisms from Gnoffo (1989)

rates = {
    {
        mechanisms = {
            {
                type = 'VT_exchange',
                p_name = 'N2',
                relaxation_time = {
                    type = 'VT_MillikanWhite_HTC',
                    HTCS_model = {
                        type = 'Park',
                        sigma_dash = 3.0e-17
                    },
                    p_name = 'N2',
                    q_names = { 'N2', 'N' },
                    a_values = {  -1,   -1 },
                    b_values = {  -1,   -1 }
                }
            },
            {
                type = 'ET_exchange',
                relaxation_time = {
                    type = 'ET_AppletonBray',
                    ions = {
                        { c_name = 'N_plus', },
                    },
                    neutrals = {
                        { c_name = 'N', sigma_coefficients = { 5.0e-20, 0.0, 0.0 } },
                        { c_name = 'N2', sigma_coefficients = { 7.5e-20, 5.5e-24, -1.0e-28 } },
                    }
                }
            }   
        }
    }
}

equilibriation_mechanisms = {}
