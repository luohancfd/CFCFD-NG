from exact_binary_diffusion import calc_rho, calc_f_N2

def ref_function(x, y, z, t):
    rho = calc_rho(x, t)
    f_N2 = calc_f_N2(rho)
    return {"rho":rho, "massf[0]":f_N2, "massf[1]":1-f_N2}

