#!/usr/bin/env python

from gaspy import *
from librad2 import *

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
    
def spectra_at_position( fname, x_star, rsm, gm=None ):
    Q,u = extract_flowstate(fname,x_star,gm)
    X = CoeffSpectra(rsm)
    rsm.radiative_spectra_for_gas_state(Q,X)
    return X