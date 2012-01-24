/*Copyright (C) Joseph Tang, 2007
                                                                                                        
    This file is part of OctVCE.
                                                                                                        
    OctVCE is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
                                                                                                        
    OctVCE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
                                                                                                        
    You should have received a copy of the GNU General Public License
    along with OctVCE.  If not, see <http://www.gnu.org/licenses/>.
*/

/**\file Source file for thermodynamic evaluations e.g. sound speed, mixture quantities etc*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ov_kernel.h"
#include "ov_thermo.h"

extern Gas_dat gas_data;
extern Jwlb_coeff JWLB_coeffs;

/*------------------------------------------------------------------*/
/**\brief Get sound speed of mixture*/

double get_mixture_soundspd(State_vector *sv)
{
  short int i;
  double R, R1, v, prod_m_frac, L, T, soundspd, e;
  double rho_amb, Cv_mix, R_mix, R_p, R_a, dTdrp, dTdra, dedrp, dedra, dPdT_p, dPdT_a, amb_m_frac, dTdr, dedr, dTdP, dedP;
  double f = 0;
  double q = 0;
  double dqdr = 0;
  double dLdr = 0;

  double vw = 0;
  double expRv[5] = {0};
  double expRLv[5] = {0};
  double expRL[5] = {0};
  soundspd = 0;
  prod_m_frac = (sv->rhof_prod)/(sv->Rhof);

  if(gas_data.eos_prod == PERFECT_GAS) {
    R = (gas_data.Cv_prod)*(gas_data.gam_prod - 1.0); /*Explosion products gas R*/
    R1 = (gas_data.Cv_amb)*(gas_data.gam_amb - 1.0);
    
    R = prod_m_frac*R + (1.0-prod_m_frac)*R1; /*Mixture gas R*/
    soundspd = SQRT((1.0+R/((gas_data.Cv_prod)*prod_m_frac+(gas_data.Cv_amb)*(1.0-prod_m_frac)))*(sv->Pres)/(sv->Rhof));
  }
  else if(gas_data.eos_prod == JWLB) {

    if(sv -> rhof_prod > JWLB_RHO_CUTOFF) { /*No sense in using JWLB if no products in cell*/
      
      if((gas_data.jwlb_cutoff == 'y') && (treat_products_as_air(sv->rhof_prod,sv->Rhof) == TRUE)) 
	soundspd = SQRT((gas_data.gam_amb)*(sv->Pres)/(sv->Rhof));
      else {
	v = (1.0/(sv->rhof_prod))/(JWLB_coeffs.v0); /*Relative specific volume*/

	amb_m_frac = 1.0-prod_m_frac;	
	rho_amb = amb_m_frac*(sv->Rhof);
	Cv_mix = prod_m_frac*(gas_data.Cv_prod) + amb_m_frac*(gas_data.Cv_amb);
      
	/*Precompute quantities involving expensive power operators*/
	
	R = (gas_data.Cv_amb)*(gas_data.gam_amb - 1.0); /*Ambient gas R*/

	for(i = 0; i < 5; i++) {
	  if(JWLB_coeffs.R[i] > 0)
	    expRv[i] = exp(-(JWLB_coeffs.R[i])*v);

	  if(JWLB_coeffs.RL[i] > 0) 
	    expRLv[i] = exp(-(JWLB_coeffs.RL[i])*v);
	}

	/*Get Gruneisen coefficient*/
	L = JWLB_coeffs.omega;
	if(JWLB_coeffs.RL[0] > 0)
	  for(i = 0; i < 5; i++) 
	    if(JWLB_coeffs.RL[i] > 0)
	      L += (JWLB_coeffs.AL[i]*v + JWLB_coeffs.BL[i])*expRLv[i];    
  
	R_a = rho_amb*R; /*Species R multiplied by corresponding species density*/
	R_p = L*(sv->rhof_prod)*(gas_data.Cv_prod);
	R_mix = R_a + R_p;
	
	if(JWLB_coeffs.RL[0] > 0) {
	  vw = pow(1.0/v, JWLB_coeffs.omega);

	  for(i = 0; i < 5; i++)
	    if(JWLB_coeffs.RL[i] > 0)
	      dLdr += (v/(sv->rhof_prod))*expRLv[i]*((JWLB_coeffs.RL[i])*(JWLB_coeffs.AL[i]*v + JWLB_coeffs.BL[i]) - JWLB_coeffs.AL[i]);

	  q = JWLB_q(v, expRLv, expRL);
	  dqdr = JWLB_dqdr(v, q, expRLv, expRL);

	  f = ((JWLB_coeffs.C)*(JWLB_coeffs.v0)*(1.0-q)/(JWLB_coeffs.omega))*vw;
	}

	/*Now compute temperature*/

	e = ((sv->RhoE) - 0.5*(SQR(sv->RhoU)+SQR(sv->RhoV)+SQR(sv->RhoW))/(sv->Rhof))/(sv->Rhof); /*Internal energy of mixture*/

	for(i = 0; i < 5; i++)
	  if(JWLB_coeffs.R[i] > 0)
	    f += ((JWLB_coeffs.A[i])*(JWLB_coeffs.v0)/(JWLB_coeffs.R[i]))*expRv[i];

	T = (e - prod_m_frac*f)/Cv_mix; /*Temperature of mixture*/
	
	/*Now obtain dT/d(rhof_prod)*/
	f = (gas_data.Cv_prod)*T*(L + (sv->rhof_prod)*dLdr); /*Reset temp variable f*/

	if(JWLB_coeffs.RL[0] > 0)
	  f += (JWLB_coeffs.C)*((JWLB_coeffs.v0)*((JWLB_coeffs.omega)+1.0)*vw*(1.0 - L*q/(JWLB_coeffs.omega)) - 
				(vw/(v*(JWLB_coeffs.omega)))*(q*dLdr + L*dqdr));
	for(i = 0; i < 5; i++)
	  if(JWLB_coeffs.R[i] > 0)
	    f += ((JWLB_coeffs.A[i])*(JWLB_coeffs.R[i])*v/(sv->rhof_prod))*expRv[i];

	dTdrp = -f/R_mix;

	/*Now obtain dT/d(rhof_amb)*/
	dTdra = -T*R/R_mix;

	/*Now obtain de/d(rhof_prod)*/
	if(JWLB_coeffs.RL[0] > 0)
	  f = ((JWLB_coeffs.C)*(JWLB_coeffs.v0)/(JWLB_coeffs.omega))*vw*v*((JWLB_coeffs.v0)*(JWLB_coeffs.omega)*(1.0-q) - dqdr/v);
	else f = 0;

	for(i = 0; i < 5; i++)
	  if(JWLB_coeffs.R[i] > 0)
	    f += ((JWLB_coeffs.A[i])/SQR(sv->rhof_prod))*expRv[i];

	dedrp = f;

	/*Now obtain de/d(rhof_amb)*/
	dedra = 0;

	/*Now obtain dP(prod)/dT*/
	dPdT_p = R_p;

	/*Now obtain dP(amb)/dT*/
	dPdT_a = R_a;

	/*Now get dT/dr for mixture (at const. pressure)*/
	dTdr = prod_m_frac*dTdrp + amb_m_frac*dTdra;

	/*Now get de/dr for mixture (at const. pressure)*/
	dedr = SQR(prod_m_frac)*dedrp + Cv_mix*dTdr; /*dedra = 0*/

	/*Now get dTdP for mixture*/
	dTdP = 1.0/(dPdT_p + dPdT_a);

	/*Now get dedP*/
	dedP = Cv_mix*dTdP;

	/*Finally get sqaured sound speed*/
	soundspd = ((sv->Pres) - SQR(sv->Rhof)*dedr)/(SQR(sv->Rhof)*dedP);

#if 1	
	  if(soundspd <= 0) {
	    printf("Imaginary soundspd!  soundspd^2 %e\n",soundspd);
	    printf("v0 %e, v %e, L %e, dL/dr %e\n",JWLB_coeffs.v0,v,L,dLdr);
	    printf("sv - r %e, rp %e, ru %e, rv %e, rw %e, re %e, p %e, e %e\n",sv->Rhof,sv->rhof_prod,sv->RhoU,sv->RhoV,sv->RhoW,
		   sv->RhoE,sv->Pres,e);
	    
	    abort();
	  }
	  else if(soundspd > 0) {
	  }
	  else {
	    printf("Soundspd stuffed!  soundspd^2 %e\n",soundspd);
	    printf("v0 %e, v %e, L %e, dL/dr %e\n",JWLB_coeffs.v0,v,L,dLdr);
	    printf("sv - r %e, rp %e, ru %e, rv %e, rw %e, re %e, p %e, e %e\n",sv->Rhof,sv->rhof_prod,sv->RhoU,sv->RhoV,sv->RhoW,
		   sv->RhoE,sv->Pres,e);
	    
	    abort();
	  }
#endif

	soundspd = SQRT(soundspd);
      }
    }
    else soundspd = SQRT((gas_data.gam_amb)*(sv->Pres)/(sv->Rhof)); /*Ideal gas sound speed for ambient gas*/
  }
  
  return(soundspd); /*Finally can get the sound speed*/	
}

/*------------------------------------------------------------------*/
/**\brief Get total intensive energy of mixture.  Density and pressure must thus be known*/

double get_mixture_RhoE(State_vector *sv, int Ton)
{
  short int i;
  double v, prod_m_frac, L, term1, T, R, R1, RhoE, rho_amb, Cv_mix;
  double q = 0;
  double vw = 0;
  double f = 0;
  double term2 = 0;
  
  double expRv[5] = {0};
  double expRLv[5] = {0};
  double expRL[5] = {0};
  T = 0; RhoE = 0;
  R = (gas_data.Cv_amb)*((gas_data.gam_amb) - 1.0); /*Ambient gas R*/
  prod_m_frac = (sv->rhof_prod)/(sv->Rhof); 

  if(gas_data.eos_prod == PERFECT_GAS) {
    R1 = (gas_data.Cv_prod)*((gas_data.gam_prod) - 1.0); /*Explosion products gas R*/
    R = R1*prod_m_frac + (1.0-prod_m_frac)*R; /*Mixture R*/
    
    T = (sv->Pres)/(R*(sv->Rhof));
    RhoE = (sv->Rhof)*(prod_m_frac*(gas_data.Cv_prod) + (1.0-prod_m_frac)*(gas_data.Cv_amb))*T + 
      0.5*(SQR(sv->RhoU) + SQR(sv->RhoV) + SQR(sv->RhoW))/(sv->Rhof);
  }
  else if(gas_data.eos_prod == JWLB) {

    if(sv -> rhof_prod > JWLB_RHO_CUTOFF) {
      
      if((gas_data.jwlb_cutoff == 'y') && (treat_products_as_air(sv->rhof_prod,sv->Rhof) == TRUE)) 
	RhoE = (sv->Pres)/((gas_data.gam_amb)-1.0) + 0.5*(SQR(sv->RhoU) + SQR(sv->RhoV) + SQR(sv->RhoW))/(sv->Rhof);
      else {
	v = (1.0/(sv->rhof_prod))/(JWLB_coeffs.v0); /*Relative specific volume*/
  
	rho_amb = (1.0-prod_m_frac)*(sv->Rhof); /*Ambient gas density*/
	Cv_mix = prod_m_frac*(gas_data.Cv_prod) + (1.0-prod_m_frac)*(gas_data.Cv_amb);

	/*First precompute some quantities involving expensive power operators*/

	for(i = 0; i < 5; i++) {
	  if(JWLB_coeffs.R[i] > 0)
	    expRv[i] = exp(-(JWLB_coeffs.R[i])*v);

	  if(JWLB_coeffs.RL[i] > 0) 
	    expRLv[i] = exp(-(JWLB_coeffs.RL[i])*v);
	}
		
	/*Now compute temperature*/

	L = JWLB_coeffs.omega;
	if(JWLB_coeffs.RL[0] > 0)
	  for(i = 0; i < 5; i++)
	    if(JWLB_coeffs.RL[i] > 0)
	      L += ((JWLB_coeffs.AL[i])*v + (JWLB_coeffs.BL[i]))*expRLv[i];

	term1 = 0;
	for(i = 0; i < 5; i++)
	  if(JWLB_coeffs.R[i] > 0)
	    term1 += (JWLB_coeffs.A[i])*expRv[i];

	if(JWLB_coeffs.RL[0] > 0) { /*Test if JWLB or JWL form used*/
	  q = JWLB_q(v, expRLv, expRL);
	  vw = pow(1.0/v, JWLB_coeffs.omega);
	  term2 = (JWLB_coeffs.C)*(1.0 - L*q/(JWLB_coeffs.omega))*vw/v;
	  f = ((JWLB_coeffs.C)*(JWLB_coeffs.v0)*(1.0 - q)/(JWLB_coeffs.omega))*vw;
	}
	
	T = ((sv->Pres)-term1-term2)/(L*(gas_data.Cv_prod)*(sv->rhof_prod) + rho_amb*R); /*Mixture temperature*/
	
	/*Use this to compute energy*/
		
	for(i = 0; i < 5; i++)
	  if(JWLB_coeffs.R[i] > 0)
	    f += ((JWLB_coeffs.A[i])*(JWLB_coeffs.v0)/(JWLB_coeffs.R[i]))*expRv[i];

	RhoE =  (sv->Rhof)*(Cv_mix*T + prod_m_frac*f) +
	  0.5*(SQR(sv->RhoU) + SQR(sv->RhoV) + SQR(sv->RhoW))/(sv->Rhof); /*Must remember to add kinetic energy*/
      }
    }
    else RhoE = (sv->Pres)/((gas_data.gam_amb)-1.0) + 0.5*(SQR(sv->RhoU) + SQR(sv->RhoV) + SQR(sv->RhoW))/(sv->Rhof);
  }

  if(Ton == 1) {
    return(T);
  }

  return(RhoE); 
}

/*------------------------------------------------------------------*/
/**\brief Get total static pressure of mixture.  Density and temperature/emergy must thus be known*/

void get_mixture_Pres(State_vector *sv)
{
  short int i;
  double v, prod_m_frac, L, term1, R, R1, e, rho_amb, Cv_mix;
  double term2 = 0;
  double f = 0;
  double q = 0;
  double vw = 0;

  double expRv[5] = {0};
  double expRLv[5] = {0};
  double expRL[5] = {0};
  
  R = (gas_data.Cv_amb)*((gas_data.gam_amb) - 1.0); /*Ambient gas R*/
  prod_m_frac = (sv->rhof_prod)/(sv->Rhof); 

  if(gas_data.eos_prod == PERFECT_GAS) {
    R1 = (gas_data.Cv_prod)*((gas_data.gam_prod) - 1.0);
    R = prod_m_frac*R1 + (1.0-prod_m_frac)*R;
    sv -> Pres = ((sv->RhoE) - 0.5*(SQR(sv->RhoU) + SQR(sv->RhoV) + SQR(sv->RhoW))/(sv->Rhof))*R/
      (prod_m_frac*(gas_data.Cv_prod) + (1.0-prod_m_frac)*(gas_data.Cv_amb));
    sv -> T = (sv -> Pres)/(R*(sv -> Rhof));
  }
  else if(gas_data.eos_prod == JWLB) {
    
    if(sv -> rhof_prod > JWLB_RHO_CUTOFF) {
      
      if((gas_data.jwlb_cutoff == 'y') && (treat_products_as_air(sv->rhof_prod,sv->Rhof) == TRUE)) {
	sv -> Pres = ((sv->RhoE) - 0.5*(SQR(sv->RhoU) + SQR(sv->RhoV) + SQR(sv->RhoW))/(sv->Rhof))*((gas_data.gam_amb) - 1.0);
	sv -> T = (sv -> Pres)/(R*(sv -> Rhof));
      }
      else {
	v = (1.0/(sv->rhof_prod))/(JWLB_coeffs.v0); /*Relative specific volume*/

	rho_amb = (1.0-prod_m_frac)*(sv->Rhof); /*Ambient gas density*/
	Cv_mix = prod_m_frac*(gas_data.Cv_prod) + (1.0-prod_m_frac)*(gas_data.Cv_amb);

	/*Precompute quantities involving expensive power operators*/
  
	for(i = 0; i < 5; i++) {
	  if(JWLB_coeffs.R[i] > 0)
	    expRv[i] = exp(-(JWLB_coeffs.R[i])*v);

	  if(JWLB_coeffs.RL[i] > 0) 
	    expRLv[i] = exp(-(JWLB_coeffs.RL[i])*v);
	}
	
	/*Now compute temperature*/

	L = JWLB_coeffs.omega;
	if(JWLB_coeffs.RL[0] > 0)
	  for(i = 0; i < 5; i++)
	    if(JWLB_coeffs.RL[i] > 0)
	      L += ((JWLB_coeffs.AL[i])*v + (JWLB_coeffs.BL[i]))*expRLv[i];

	term1 = 0;
	for(i = 0; i < 5; i++)
	  if(JWLB_coeffs.R[i] > 0)
	    term1 += ((JWLB_coeffs.A[i])*(JWLB_coeffs.v0)/(JWLB_coeffs.R[i]))*expRv[i];

	if(JWLB_coeffs.RL[0] > 0) { /*I.e. if JWLB form, not JWL form, is actually used*/
	  vw = pow(1.0/v, JWLB_coeffs.omega);
	  q = JWLB_q(v, expRLv, expRL);
	  term2 = ((JWLB_coeffs.C)*(JWLB_coeffs.v0)*(1.0 - q)/(JWLB_coeffs.omega))*vw;
	  f = (JWLB_coeffs.C)*(1.0 - L*q/(JWLB_coeffs.omega))*vw/v;
	}
	
	e = ((sv->RhoE) - 0.5*(SQR(sv->RhoU) + SQR(sv->RhoV) + SQR(sv->RhoW))/(sv->Rhof))/(sv->Rhof); /*Internal energy e*/
	/*Must remember to subtract kinetic energy*/

	sv -> T = (e - prod_m_frac*(term1+term2))/Cv_mix; /*Temperature of mixture*/
	
	for(i = 0; i < 5; i++)
	  if(JWLB_coeffs.R[i] > 0)
	    f += (JWLB_coeffs.A[i])*expRv[i];

	/*Using pressure-dependent form for mixtures*/
	sv -> Pres = (rho_amb*R + L*(gas_data.Cv_prod)*(sv -> rhof_prod))*(sv -> T) + f;

#if 0
	if(sv->T < 0) {
	  printf("In get_mix_Pres() is -ve temp\n");
	  printf("v %g\n",v);
	  printf("P %9.8e, term1, %9.8e, term2 %9.8e, L %g, rhof_prod %9.8e, rho_amb %g, R %9.8e\n",sv->Pres,term1,term2,L,sv->rhof_prod,rho_amb,R);
	  printf("expRv[0] %g, expRv[1] %g\n", expRv[0], expRv[1]);
	  printf("m frac %g\n",prod_m_frac);
	  if(sv->rhof_prod > sv->Rhof)
	    printf("prod m is too much\n");
	  abort();
	}
#endif
      }
    }
    else {
      sv -> Pres = ((sv->RhoE) - 0.5*(SQR(sv->RhoU) + SQR(sv->RhoV) + SQR(sv->RhoW))/(sv->Rhof))*((gas_data.gam_amb) - 1.0);
      sv -> T = (sv -> Pres)/(R*(sv -> Rhof));
    }
  }
}

/*------------------------------------------------------------------*/
/**\brief Get temperature of mixture (I really only use this currently when computing the EFM flux)*/

double get_mixture_T(double Rhof, double rhof_prod, double RhoE, double Pres, double U, double V, double W)
{
  short int i;
  double v, prod_m_frac, L, term1, T, R, R1, e, rho_amb, Cv_mix;
  double term2 = 0;
  double q = 0;
  double vw = 0;

  double expRv[5] = {0.0};
  double expRLv[5] = {0.0};
  double expRL[5] = {0.0};
  T = 0;
  R = (gas_data.Cv_amb)*((gas_data.gam_amb) - 1.0); /*Ambient gas R*/
  prod_m_frac = (rhof_prod)/(Rhof); 

  if(gas_data.eos_prod == PERFECT_GAS) {
    R1 = (gas_data.Cv_prod)*((gas_data.gam_prod) - 1.0);
    R = prod_m_frac*R1 + (1.0-prod_m_frac)*R;
    T = Pres/(R*Rhof);
  }
  else if(gas_data.eos_prod == JWLB) {

    if(rhof_prod > JWLB_RHO_CUTOFF) {

      if((gas_data.jwlb_cutoff == 'y') && (treat_products_as_air(rhof_prod,Rhof) == TRUE)) 
	T = Pres/(R*Rhof); /*Just treat as air (prod_m_frac too small to be of significance)*/
      else {
	v = (1.0/rhof_prod)/(JWLB_coeffs.v0); /*Relative specific volume*/

	rho_amb = (1.0-prod_m_frac)*(Rhof); /*Ambient gas density*/
	Cv_mix = prod_m_frac*(gas_data.Cv_prod) + (1.0-prod_m_frac)*(gas_data.Cv_amb);

	/*Precompute quantities involving expensive power operators*/
  
	for(i = 0; i < 5; i++) {
	  if(JWLB_coeffs.R[i] > 0)
	    expRv[i] = exp(-(JWLB_coeffs.R[i])*v);
      
	  if(JWLB_coeffs.RL[i] > 0) 
	    expRLv[i] = exp(-(JWLB_coeffs.RL[i])*v);
	}
	
	/*Now compute temperature*/
    
	L = JWLB_coeffs.omega;
	if(JWLB_coeffs.RL[0] > 0)
	  for(i = 0; i < 5; i++)
	    if(JWLB_coeffs.RL[i] > 0)
	      L += ((JWLB_coeffs.AL[i])*v + (JWLB_coeffs.BL[i]))*expRLv[i];

	term1 = 0;
	for(i = 0; i < 5; i++)
	  if(JWLB_coeffs.R[i] > 0)
	    term1 += ((JWLB_coeffs.A[i])*(JWLB_coeffs.v0)/(JWLB_coeffs.R[i]))*expRv[i];
    
	if(JWLB_coeffs.RL[0] > 0) { /*I.e. if JWLB form, not JWL form, is actually used*/
	  vw = pow(1.0/v, JWLB_coeffs.omega);
	  q = JWLB_q(v, expRLv, expRL);
	  term2 = ((JWLB_coeffs.C)*(JWLB_coeffs.v0)*(1.0 - q)/(JWLB_coeffs.omega))*vw;
	}
	
	e = (RhoE/Rhof) - 0.5*(SQR(U) + SQR(V) + SQR(W)); /*Internal energy e*/
	/*Must remember to subtract kinetic energy*/

	T = (e - prod_m_frac*(term1+term2))/Cv_mix; /*Temperature of mixture*/
      }
    }
    else T = Pres/(R*Rhof); /*Just treat as air - prod_m_frac way too small to be of use*/
  }

  return(T);
}

/*------------------------------------------------------------------*/
/**\brief Compute the 'q' term appearing in JWLB EOS when temperature exists*/

double JWLB_q(double V, double expRLv[], double expRL[]) /**< V is normalized specific volume v/v0*/ 
{
  short int i;
  double q, RL, AL, BL;
  double mu = (1.0/V)-1.0;

  q = 0;
  for(i = 0; i < 5; i++) {
    AL = JWLB_coeffs.AL[i];
    BL = JWLB_coeffs.BL[i];
    RL = JWLB_coeffs.RL[i];

    if(RL > 0) {
      expRL[i] = exp(-RL);
      q += (AL/RL)*expRLv[i] + BL*mu*expRL[i]*JWLB_K(i,mu);
    }
  }

  q = exp(q);  

  return(q);
}

/*------------------------------------------------------------------*/
/**\brief Precompute terms in JWLB 'dqdr' expression*/

double JWLB_dqdr(double V, double q, double expRLv[], double expRL[])
{
  short int i;
  double AL, BL, dqdr;
  double mu = (1.0/V)-1.0;

  dqdr = 0;
  for(i = 0; i < 5; i++) {
    AL = JWLB_coeffs.AL[i];
    BL = JWLB_coeffs.BL[i];
    if(JWLB_coeffs.RL[i] > 0) 
      dqdr += AL*SQR(V)*expRLv[i] + BL*expRL[i]*(JWLB_K(i,mu) + JWLB_DKDM(i,mu)); 
  }

  /*Currently dqdr is dqdm/q*/
  dqdr = (JWLB_coeffs.v0)*q*dqdr;

  return(dqdr);
}

/*------------------------------------------------------------------*/
/**\brief Precompute terms appearing in JWLB 'q' expression*/

void precompute_JRL(void)
{
  short int i;
  double RL;
  
  for(i = 0; i < 5; i++) {
    RL = JWLB_coeffs.RL[i];
    if(JWLB_coeffs.RL[i] > 0) {
      JWLB_coeffs.JRL[i][0] = RL - 1.0;
      JWLB_coeffs.JRL[i][1] = SQR(RL) - 4.0*(RL) + 2.0;
      JWLB_coeffs.JRL[i][2] = CUBE(RL) - 9.0*SQR(RL) + 18.0*RL - 6.0;
      JWLB_coeffs.JRL[i][3] = QUART(RL) - 16.0*CUBE(RL) + 72.0*SQR(RL) - 96.0*RL + 24.0;
    }
  }  
}

/*------------------------------------------------------------------*/
/**\brief Test if explosion products can be treated as air (for faster computation)*/

short int treat_products_as_air(double rho_p, double rho)
{
  double Rmin, RLmin, A, AL, V;

#if 0
  if(rho_p/rho < AIR_THRESHOLD)
    return(FALSE); /*Maybe this is an even simpler test - but this seems redundant as I've got already a JWLB_RHOF_CUTOFF variable.
		     Nonetheless this centralizes the test so I don't need to repeat it everywhere*/
#endif

  UNUSED_VARIABLE(rho);

  Rmin = JWLB_coeffs.R[(JWLB_coeffs.smallest_R_index)];
  RLmin = JWLB_coeffs.RL[(JWLB_coeffs.smallest_RL_index)];
  A = JWLB_coeffs.A[(JWLB_coeffs.smallest_R_index)];
  AL = JWLB_coeffs.AL[(JWLB_coeffs.smallest_RL_index)];

  V = (1.0/rho_p)/(JWLB_coeffs.v0); /*Normalized specific volume*/
  
  /*Perhaps at the very initial stages of the blast (say < 100 timesteps) we still use JWLB?*/

  if(((A*exp(-Rmin*V)) < AIR_THRESHOLD) && ((AL*V*exp(-RLmin*V)) < AIR_THRESHOLD)) 
    return(TRUE); /*Exponential and Gruneisen terms are small enough - can treat products as air*/
  
  return(FALSE);
}

/*------------------------------------------------------------------*/

/**\brief Find smallest values of R and RL for future reference*/

void get_smallest_JWLB_terms(void)
{
  short int i;
  double R, RL;

  JWLB_coeffs.smallest_R_index = 0; /*Initialize*/
  JWLB_coeffs.smallest_RL_index = 0;

  R = JWLB_coeffs.R[0];
  for(i = 1; i < 5; i++) {
    if((JWLB_coeffs.R[i] < R) && (JWLB_coeffs.R[i] > 0)) {
      R = JWLB_coeffs.R[i];
      JWLB_coeffs.smallest_R_index = i;
    }
  }

  RL = JWLB_coeffs.RL[0];
  for(i = 1; i < 5; i++) {
    if((JWLB_coeffs.RL[i] < RL) && (JWLB_coeffs.RL[i] > 0)) {
      RL = JWLB_coeffs.RL[i];
      JWLB_coeffs.smallest_RL_index = i;
    }
  }  
}

/*------------------------------------------------------------------*/
