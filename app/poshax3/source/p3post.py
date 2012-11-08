#!/usr/bin/env python

from gaspy import *
from librad2 import *
import sys

def ndim_from_header( fname ):
    ifile = open( fname, "r" )
    lines = ifile.readlines()
    ifile.close()
    ntm = 0
    nsp = 0
    species = []
    for line in lines:
        if "T[" in line: ntm+=1
        if "massf" in line:
            i = line.find("-")
            species.append( line[i+1:-1].replace("\n","") )
            nsp+=1
    return nsp,ntm,species

def extract_flowstate( fname, x_star, gm=None ):
    if gm!=None:
        nsp = gm.get_number_of_species()
        ntm = gm.get_number_of_modes()
        species = []
        for isp in range(nsp):
            species.append( gm.species_name(isp) )
    else:
        nsp,ntm,species = ndim_from_header(fname)
    
    ifile = open( fname, "r" )
    lines = ifile.readlines()
    ifile.close()
    
    above = False
    
    for line in lines:
        if len(line)==0: continue
        if line[0]=="#": continue
        tks = line.split()
        x = float(tks[0])
        if x >= x_star: break
        
    Q = Gas_data(nsp,ntm)
    for itm in range(ntm):
        Q.T[itm] = float(tks[1+itm])
    Q.p = float(tks[1+ntm])
    Q.rho = float(tks[1+ntm+1])
    for isp in range(nsp):
        Q.massf[isp] = float(tks[1+ntm+3+isp])
    
    if "e_minus" in species:
        R_e = PC_k_SI / PC_m_SI
        Q.p_e = Q.rho * Q.massf[-1] * R_e * Q.T[-1]
    
    if gm!=None:
        gm.eval_thermo_rhoT(Q)
    
    u = float(tks[1+ntm+2])
    
    return Q,u

def extract_averaged_flowstate( fname, x_star, delta_x, gm=None ):
    if gm!=None:
        nsp = gm.get_number_of_species()
        ntm = gm.get_number_of_modes()
        species = []
        for isp in range(nsp):
            species.append( gm.species_name(isp) )
    else:
        nsp,ntm,species = ndim_from_header(fname)
    
    ifile = open( fname, "r" )
    lines = ifile.readlines()
    ifile.close()
    
    x_list = []
    Q_list = []
    u_list = []
    
    for line in lines:
        if len(line)==0: continue
        if line[0]=="#": continue
        tks = line.split()
        x = float(tks[0])
        if x >= ( x_star - delta_x/2.0 ) and x <= ( x_star + delta_x/2.0 ):
            Q = Gas_data(nsp,ntm)
            for itm in range(ntm):
                Q.T[itm] = float(tks[1+itm])
            Q.p = float(tks[1+ntm])
            Q.rho = float(tks[1+ntm+1])
            for isp in range(nsp):
                Q.massf[isp] = float(tks[1+ntm+3+isp])
            if "e_minus" in species:
                R_e = PC_k_SI / PC_m_SI
                Q.p_e = Q.rho * Q.massf[-1] * R_e * Q.T[-1]
            if gm!=None:
                gm.eval_thermo_rhoT(Q)
            u = float(tks[1+ntm+2])
            # append to lists
            x_list.append(x)
            Q_list.append(Q)
            u_list.append(u)
            
    print x_list
               
                
    # do averaging
    Q = Gas_data(nsp,ntm); u = 0.0; dx_sum = 0.0
    for i in range(len(x_list)):
        if i==0: dx = x_list[1] - x_list[0]
        elif i==len(x_list)-1: dx = x_list[-1] - x_list[-2]
        else:  dx = 0.5 * ( x_list[i+1] - x_list[i-1] )
        for itm in range(ntm):
            Q.T[itm] += Q_list[i].T[itm]*dx
        Q.p += Q_list[i].p*dx
        Q.rho += Q_list[i].rho*dx
        for isp in range(nsp):
            Q.massf[isp] += Q_list[i].massf[isp]*dx
        u += u_list[i]*dx
        dx_sum += dx    
        
    print "dx_sum = ", dx_sum
    
    for itm in range(ntm):
        Q.T[itm] /= dx_sum
    Q.p /= dx_sum
    Q.rho /= dx_sum
    massf_sum = 0.0
    for isp in range(nsp):
        Q.massf[isp] /= dx_sum
        massf_sum += Q.massf[isp]
    for isp in range(nsp):
        Q.massf[isp] /= massf_sum
    u /= dx_sum
    if "e_minus" in species:
        R_e = PC_k_SI / PC_m_SI
        Q.p_e = Q.rho * Q.massf[-1] * R_e * Q.T[-1]
    if gm!=None:
        gm.eval_thermo_rhoT(Q)
    
    return Q,u
    
def spectra_at_position( fname, x_star, rsm, gm=None ):
    Q,u = extract_flowstate(fname,x_star,gm)
    X = CoeffSpectra(rsm)
    rsm.radiative_spectra_for_gas_state(Q,X)
    return X