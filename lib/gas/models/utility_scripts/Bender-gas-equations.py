#!/usr/bin/env python
"""
Bender-gas-equations.py

Program to generate the equations for the Bender real gas model using the
symbolic math module, sympy, to perform the analytical derivatives and integrals.

It is reasonably straight forward to perform the calculus by hand, but this
script is used to alleviate maintenance and transcription issues. The generated
output equations may not be fully factorised or simplified but are simply
pasted into the corresponding Bender gas model functions.

Refer to Reynolds, WC (1970). Thermodynamic Properties in SI.
for a description of the equations and calculus performed in this module.

.. Author: Peter Blyton
.. Version: 29 March 2012
"""

from sympy.utilities.codegen import codegen
from sympy import symbols, diff, integrate, exp, latex
from sympy.abc import rho, R, T
Tinv = T**(-1)
T0 = symbols("T0")

def bender_EOS():
    """Returns a sympy expression for the Bender p-rho-T equation of state."""
    A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20 \
    = symbols("A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20")
    return rho*R*T \
        + rho**2*(A1*T + A2 + Tinv*(A3 + Tinv*(A4 + Tinv*A5))) \
        + rho**3*(A6*T + A7 + Tinv*A8) \
        + rho**4*(A9*T + A10) \
        + rho**5*(A11*T + A12) \
        + rho**6*A13 \
        + rho**3*Tinv*Tinv*(A14 + Tinv*(A15 + Tinv*A16))*exp(-A20*rho*rho) \
        + rho**5*Tinv*Tinv*(A17 + Tinv*(A18 + Tinv*A19))*exp(-A20*rho*rho)

def Cv0_polynomial():
    """Returns a sympy expression for the classic fourth order ideal gas
    heat capacity polynomial."""
    G1,G2,G3,G4,G5,G6 = symbols("G1,G2,G3,G4,G5,G6")
    return G1*Tinv + G2 + G3*T + G4*T**2 + G5*T**3 + G6*T**4

def create_code_str(expression):
    """Turn a sympy expression into a string in C code syntax."""
    [(c_name, c_code), (h_name, c_header)] = \
        codegen(("f", expression), "C", "test", header=False, empty=False)
    # Codegen generates an entire C function, just grab the equation from it.
    return c_code.splitlines(True)[3].lstrip()

def create_cpp_syntax(code_str, coeff):
    """Turn the equation coefficients contained within the string into C++ vector syntax."""
    # Replace from large to small numbered coefficient so that A12 is not changed to A[1]2
    for i in range(code_str.count(coeff), 0, -1):
        code_str = code_str.replace(str(coeff) + str(i), str(coeff) + "[" + str(i) + "]")
    return code_str

def create_latex_str(expression):
    """"Turn the sympy expression into a latex formula in a dmath
    environment, which auto-breaks long latex equations."""
    latex_eqn = latex(expression, mode="equation")
    latex_eqn = latex_eqn.replace("equation", "dmath")
    latex_eqn = latex_eqn.replace("\operatorname{log}", "\log")
    return latex_eqn

if __name__ == "__main__":
    print "Creating Bender gas model equations..."
    P = bender_EOS()
    Cv0 = Cv0_polynomial()

    print "Differentiating and integrating equations..."
    dPdT = diff(P, T)
    dPdrho = diff(P, rho)
    u_Cv0 = integrate(Cv0, (T, T0, T))
    s_Cv0 = integrate(Cv0/T, (T, T0, T))
    u_dPdT = integrate(rho**(-2)*(P - T*dPdT), (rho, 0, rho))
    s_dPdT = integrate(rho**(-2)*(rho*R - dPdT), (rho, 0, rho))

    print "Writing output to file..."
    f = open("Bender-gas-equations.txt", "w")
    f.write("**************************************************\n")
    f.write("Sympy generated equations for the Bender gas model\n")
    f.write("**************************************************\n\n")
    f.write("Equations in C++ format\n")
    f.write("=======================\n\n")
    f.write("The Bender p-rho-T equation\n")
    f.write("---------------------------\n")
    f.write(create_cpp_syntax(create_code_str(P), "A") + "\n")
    f.write("dPdT\n")
    f.write("----\n")
    f.write(create_cpp_syntax(create_code_str(dPdT), "A") + "\n")
    f.write("dPdrho\n")
    f.write("------\n")
    f.write(create_cpp_syntax(create_code_str(dPdrho), "A") + "\n")
    f.write("Integral of Cv0\n")
    f.write("---------------\n")
    f.write(create_cpp_syntax(create_code_str(u_Cv0), "G") + "\n")
    f.write("Integral of Cv0/T\n")
    f.write("-----------------\n")
    f.write(create_cpp_syntax(create_code_str(s_Cv0), "G") + "\n")
    f.write("Integral of (P - T*dPdT)/rho**2\n")
    f.write("-------------------------------\n")
    f.write(create_cpp_syntax(create_code_str(u_dPdT), "A") + "\n")
    f.write("Integral of (R*T - dPdT)/rho**2\n")
    f.write("-------------------------------\n")
    f.write(create_cpp_syntax(create_code_str(s_dPdT), "A") + "\n")
    f.write("Equations in LaTeX format\n")
    f.write("=========================\n\n")
    f.write("The Bender p-rho-T equation\n")
    f.write("---------------------------\n")
    f.write(create_latex_str(P) + "\n\n")
    f.write("dPdT\n")
    f.write("----\n")
    f.write(create_latex_str(dPdT) + "\n\n")
    f.write("dPdrho\n")
    f.write("------\n")
    f.write(create_latex_str(dPdrho) + "\n\n")
    f.write("Integral of Cv0\n")
    f.write("---------------\n")
    f.write(create_latex_str(u_Cv0) + "\n\n")
    f.write("Integral of Cv0/T\n")
    f.write("-----------------\n")
    f.write(create_latex_str(s_Cv0) + "\n\n")
    f.write("Integral of (P - T*dPdT)/rho**2\n")
    f.write("-------------------------------\n")
    f.write(create_latex_str(u_dPdT) + "\n\n")
    f.write("Integral of (R*T - dPdT)/rho**2\n")
    f.write("-------------------------------\n")
    f.write(create_latex_str(s_dPdT) + "\n\n")
    f.close()
    print "Done."
