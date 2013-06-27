#! /usr/bin/python
"""
Created on Thur Apr 18 2013

authors: Nikhil Banerji
         Jeremy Mora-Monteros
"""

# The equations used in this file have been taken from
# "Molecular Theory of Gases and Liquids",
# Hirschfelder, Curtiss and Bird, in 1954
# pp 32, 527

#for Stockmayer Potential/Lennard-Jones Potential

import numpy as np
import math

# The function calc_rm_r computes the reduced distance of closest approach rm_r
# This is done by solving "1-4*((1/rm_r)^12-(1/rm_r)^6)/g_r^2-b_r^2/rm_r^2=0"
# The goal is to find the maximum limit in "int_ki[i]" in function "ki_r"
# Because we're integrating from this value.

def calc_rm_r(g_r, b_r):
    coeff = [1, 0, -b_r**2, 0, 0, 0, 4/g_r**2, 0, 0, 0, 0, 0, -4/g_r**2]
    rm_r = np.roots(coeff)
    for i in range(0,len(rm_r)):
        if np.isreal(rm_r[i])==False:
            rm_r[i] = 0
        elif rm_r[i]<0:
            rm_r[i] = 0    
    rm_r = max(np.real(rm_r))
    return rm_r + 0.0001
# Here 0.0001 has been added to the solution to avoid the integration of an
# infinite value which would give a "nan" error
# 0.0001 represents a small value compared to r_step which the integration step

# The function phi_r simply calculates the reduced Stockmayer potential

def phi_r(r_ri,delta):
    phi_r = 4*((1/r_ri)**12-(1/r_ri)**6-delta*(1/r_ri)**3)
    return phi_r

# The function ki_r calculates the reduced angle of deflection
# by integrating on r_r a function of phi_r, b_r(i),g_r(i)
# from rm_r to infinity

def ki_r(g_ri,b_ri,r_step,r_max,delta):
    r_min = calc_rm_r(g_ri,b_ri)
    r_r = np.arange(r_min,r_max+r_step,r_step)
    int_ki = [0]
    for i in range(0,len(r_r)):
        int_ki[i] = 1/((r_r[i]**2)*(1-(b_ri**2/r_r[i]**2)-(phi_r(r_r[i],delta)/g_ri**2))**0.5)
        if i < len(r_r)-1:
            int_ki.append(0)
    ki_r = math.pi-2*b_ri*np.trapz(int_ki,r_r)
    return ki_r

# The function Ql_r calculates the reduced cross-section
# by integrating on b_r a function of ki_r
# from 0 to infinity

def Ql_r(l,g_ri,b_min,b_step,b_max,r_step,r_max,delta):
    b_r = np.arange(b_min,b_max+b_step,b_step)
    int_Q = [0]
    for i in range(0,len(b_r)):
        int_Q[i] = (1-(math.cos(ki_r(g_ri,b_r[i],r_step,r_max,delta))**l))*b_r[i]
        if i < len(b_r)-1:
            int_Q.append(0)
    Ql_r = 2/(1-0.5*(1+(-1)**l)/(1+l))*np.trapz(int_Q, b_r)
    return Ql_r

# The function Omegals_r calculates the reduced collision integral
# by integrating on g_r a function of Ql_r and T_r
# from 0 to infinity

def Omegals_r(l,s,T_r,g_min,g_step,g_max,b_min,b_step,b_max,r_step,r_max,delta):
    g_r = np.arange(g_min,g_max+g_step,g_step)
    int_Omega = [0]
    for i in range(0,len(g_r)):
        int_Omega[i] = math.exp(-g_r[i]**2/T_r)*g_r[i]**(2*s+3)*Ql_r(l,g_r[i],b_min,b_step,b_max,r_step,r_max,delta)
        if i < len(g_r)-1:
            int_Omega.append(0)
    Omegals_r = 2/(math.factorial(s+1)*T_r**(s+2))*np.trapz(int_Omega, g_r)
    return Omegals_r