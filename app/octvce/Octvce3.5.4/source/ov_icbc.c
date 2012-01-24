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

/**\file Source file for setting ICs and BCs for octree cartesian cells.   
   Typically BCs of root are freestream/ambient (but may need to be ambient non-reflecting); blast wave simulations also need an 
   initial condition represented as region of high temperature and/or pressure gas*/
  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ov_kernel.h"
#include "ov_icbc.h"
#include "ov_thermo.h"

#define NOTINCLUDE 0

#if NOTINCLUDE
extern State_vector IC_inside_state;
extern State_vector IC_outside_state;
#else
extern State_vector IC_ambient;
extern State_vector **IC_states;
extern Deton_pts Deton;
#endif

extern double root_scale;
extern Gas_dat gas_data;
extern Jwlb_coeff JWLB_coeffs;
extern BC_data BC_est;
extern BC_data BC_wst;
extern BC_data BC_nth;
extern BC_data BC_sth;
extern BC_data BC_upr;
extern BC_data BC_lwr;

#if SUPERSONIC_VORTEX
extern double M_in;
#endif

/*------------------------------------------------------------------*/

/**\brief Sets state of flow opposite a boundary interface to enforce required boundary conditions.  For the moment
 only BCs are those corresponding to border conditions or solid walls; special embedded BCs may require location of cell on octree*/

State_vector set_boundary_state(State_vector *sve, State_vector *svc, double norm[], short int direction, short int BC_type, 
				short int stage, double time, double dx, Grad *Grads, double lim)
{
  /**< BC_type specifies whether wall or border BCs employed, sve is extrapolated state vector to cell face, svc is
     cell-centered values right beside face.  I don't recommend (for the moment) using two species and a specific
     or transient BC (though I do try and account for this).  Basically, explosion problems within domains should be regarded as
     different from shock waves entering domains.*/
    
  State_vector sv;
  double U_dot_norm, U[3];

  if(BC_type == WALL) 
    {
      U[0] = (sve->RhoU)/(sve->Rhof);
      U[1] = (sve->RhoV)/(sve->Rhof);
      U[2] = (sve->RhoW)/(sve->Rhof);
      U_dot_norm = U[0]*norm[0] + U[1]*norm[1] + U[2]*norm[2];

      sv.Rhof = sve -> Rhof;
      sv.RhoU = (sv.Rhof)*(U[0] - 2.0*U_dot_norm*norm[0]);
      sv.RhoV = (sv.Rhof)*(U[1] - 2.0*U_dot_norm*norm[1]);
      sv.RhoW = (sv.Rhof)*(U[2] - 2.0*U_dot_norm*norm[2]);
      sv.Pres = sve -> Pres;
      sv.RhoE = sve -> RhoE;
      
      if(gas_data.num_species > 1) {
	sv.rhof_prod = sve -> rhof_prod;
#if AFTERBURN
	sv.rhof_prod_ub = sve -> rhof_prod_ub;
#endif
      }
    }
  else if(BC_type == BORDER) /*Implement border BCs*/
    {      
      Prim_vector *pvb;
      short int i, bnorm;
      double Uxtr, Vxtr, Wxtr, M, dPdx, soundspd, R, R1, prod_m_frac, v, term1, term3, q, vel;
      double T = 0;
      double expRv[5] = {0};
      double expRLv[5] = {0};
      double expRL[5] = {0};

      BC_data *BC_now;
      pvb = NULL; Uxtr = 0; Vxtr = 0; Wxtr = 0; dPdx = 0; R = 0; 
      prod_m_frac = 0; v = 0; term1 = 0; term3 = 0; q = 0; BC_now = NULL;
      switch(direction)
	{
	case EST:
	  BC_now = &BC_est; break;
	case WST:
	  BC_now = &BC_wst; break;
	case NTH:
	  BC_now = &BC_nth; break;
	case STH:
	  BC_now = &BC_sth; break;
	case UPR:
	  BC_now = &BC_upr; break;
	case LWR:
	  BC_now = &BC_lwr; break;
	}

      if(BC_now -> type == TRANSIENT)
	{
	  if(stage == 1)
	    pvb = &(BC_now -> Pv1);
	  else if(stage == 2)
	    pvb = &(BC_now -> Pv2);

	  /*First do extrapolated values*/
	  sv.Rhof = sve -> Rhof; /*'Right' state would just have extrapolated values here as no ghost cells actually used*/
	  sv.Pres = sve -> Pres;
	  sv.RhoE = sve -> RhoE;
	  sv.RhoU = sve -> RhoU;
	  sv.RhoV = sve -> RhoV;
	  sv.RhoW = sve -> RhoW;
	      	
	  Uxtr = (sv.RhoU)/(sv.Rhof); /*Extrapolated values of velocity*/
	  Vxtr = (sv.RhoV)/(sv.Rhof);
	  Wxtr = (sv.RhoW)/(sv.Rhof);

	  if((strcmp(BC_now -> stat[BC_now -> step][1], "Extrapolated") != 0) && 
	     ((BC_now -> step >= (BC_now -> num_pts-1)) || (strcmp(BC_now -> stat[(BC_now->step)+1][1], "Extrapolated") != 0) ||
#if FLOAT_EQ_STABLE
	      EQ(strtod(BC_now->stat[BC_now->step][0], NULL), time)
#else
	      (strtod(BC_now->stat[BC_now->step][0], NULL) == time)
#endif
	      )) 
	    { /*First variable is fixed - can only be fixed from table of transient values if time coincides with particular table
		entry, or its value lies between 2 non-extrapolated transient entries*/

	      if(BC_now -> pres_spec == 1) /*Pressure is fixed*/
		sv.Pres = pvb -> P;
	      else if(BC_now -> rhof_spec == 1)
		sv.Rhof = pvb -> R;
	      else if(BC_now -> temp_spec == 1)
		T = pvb -> T;
	    }

	  if((strcmp(BC_now -> stat[BC_now -> step][2], "Extrapolated") != 0) &&
	     (((BC_now -> step >= (BC_now -> num_pts-1)) || (strcmp(BC_now -> stat[(BC_now->step)+1][2], "Extrapolated") != 0)) ||
#if FLOAT_EQ_STABLE
	      EQ(strtod(BC_now->stat[BC_now->step][0], NULL), time)
#else
	      (strtod(BC_now->stat[BC_now->step][0], NULL) == time)
#endif
	      )) 
	    { /*2nd variable is fixed*/

	      if(BC_now -> pres_spec == 2)
		sv.Pres = pvb -> P;
	      else if(BC_now -> rhof_spec == 2)
		sv.Rhof = pvb -> R;
	      else if(BC_now -> temp_spec == 2)
		T = pvb -> T;
	    }

	  R1 = (gas_data.Cv_amb)*(gas_data.gam_amb - 1.0); /*Ambient gas R*/

	  if(gas_data.num_species > 1)
	    {
	      prod_m_frac = (sve -> rhof_prod)/(sve -> Rhof); /*Mass fraction constant throughout cell*/
	      sv.rhof_prod = prod_m_frac*(sv.Rhof);
#if AFTERBURN
	      sv.rhof_prod_ub = ((sve->rhof_prod_ub)/(sve->rhof_prod))*(sv.rhof_prod);
#endif
	    }

	  if(T != 0) /*Temperature isn't extrapolated i.e. it's fixed*/
	    {
	      if(gas_data.num_species > 1)
		{		  		  
		  if(gas_data.eos_prod == PERFECT_GAS) 
		    {
		      R = (gas_data.Cv_prod)*(gas_data.gam_prod - 1.0); /*Explosion products gas R*/
		      R = prod_m_frac*R + (1.0-prod_m_frac)*R1; /*Mixture R*/    

		      if(BC_now -> rhof_spec == 0) 
			sv.Rhof = (sv.Pres)/(R*T);
		      else if(BC_now -> pres_spec == 0)
			sv.Pres = (sv.Rhof)*R*T;
		    }
		  else if(gas_data.eos_prod == JWLB)
		    {
		      /*As rho = rho(p, T) difficult to find, dependent variable is pressure.  Get pressure from partial pressures
			(using density mass fractions)*/

		      if(sv.rhof_prod > JWLB_RHO_CUTOFF) { /*No sense in implementing JWLB EOS if no explosion products in cell*/
			
			if((gas_data.jwlb_cutoff == 'y') && (treat_products_as_air(sv.rhof_prod,sv.Rhof) == TRUE))  /*Can treat products as air*/
			  sv.Pres = (sv.rhof_prod)*R1*T;
			else {
			  v = (1.0/(sv.rhof_prod))/(JWLB_coeffs.v0); /*Normalized specific volume*/

			  /*Precompute some quantities*/
			  for(i = 0; i < 5; i++) {
			    if(JWLB_coeffs.R[i] > 0)
			      expRv[i] = exp(-(JWLB_coeffs.R[i])*v);

			    if(JWLB_coeffs.RL[i] > 0) 
			      expRLv[i] = exp(-(JWLB_coeffs.RL[i])*v);
			  }
			  q = JWLB_q(v, expRLv, expRL);

			  R = JWLB_coeffs.omega;
			  for(i = 0; i < 5; i++)
			    if(JWLB_coeffs.RL[i] > 0)
			      R += ((JWLB_coeffs.AL[i])*v + (JWLB_coeffs.BL[i]))*expRLv[i];

			  term1 = 0;
			  for(i = 0; i < 5; i++)
			    if(JWLB_coeffs.R[i] > 0)
			      term1 += (JWLB_coeffs.A[i])*expRv[i];
			
			  term3 = (JWLB_coeffs.C)*(1.0 - R*q/(JWLB_coeffs.omega))*pow(v,-(1.0+(JWLB_coeffs.omega)));
		    
			  sv.Pres = term1 + R*(gas_data.Cv_prod)*T/(v*(JWLB_coeffs.v0)) + term3;
			} 
		      }
		      else sv.Pres = 0;

		      /*Partial pressure of explosion products*/

		      (sv.Pres) += (1.0-prod_m_frac)*(sv.Rhof)*R1*T; /*Total static pressure of mixture*/ 
		    }		  
		}
	      else
		{
		  if(BC_now -> rhof_spec == 0) /*Density was dependent*/
		    sv.Rhof = (sv.Pres)/(R1*T);
		  else if(BC_now -> pres_spec == 0)
		    sv.Pres = (sv.Rhof)*R1*T;
		}
	    } /*If temperature extrapolated, then dependent state variable automatically extrapolated too*/

	  if(gas_data.num_species > 1) {
	    sv.rhof_prod = ((sve -> rhof_prod)/(sve -> Rhof))*(sv.Rhof);
#if AFTERBURN
	    sv.rhof_prod_ub = ((sve->rhof_prod_ub)/(sve->rhof_prod))*(sv.rhof_prod);
#endif
	  }

	  if((strcmp(BC_now -> stat[BC_now -> step][3], "Extrapolated") != 0) &&
	     ((BC_now -> step >= (BC_now -> num_pts-1)) || (strcmp(BC_now -> stat[(BC_now->step)+1][3], "Extrapolated") != 0) ||
#if FLOAT_EQ_STABLE
	      EQ(strtod(BC_now->stat[BC_now->step][0], NULL), time)
#else
	      (strtod(BC_now->stat[BC_now->step][0], NULL) == time)
#endif	      
	      )) 
	    sv.RhoU = (pvb -> U)*(sv.Rhof); /*Get fixed value of u velocity*/
	  else sv.RhoU = (sv.Rhof)*Uxtr;

	  if((strcmp(BC_now -> stat[BC_now -> step][4], "Extrapolated") != 0) &&
	     ((BC_now -> step >= (BC_now -> num_pts-1)) || (strcmp(BC_now -> stat[(BC_now->step)+1][4], "Extrapolated") != 0) ||
#if FLOAT_EQ_STABLE
	      EQ(strtod(BC_now->stat[BC_now->step][0], NULL), time)
#else	
	      (strtod(BC_now->stat[BC_now->step][0], NULL) == time)
#endif
)) 
	    sv.RhoV = (pvb -> V)*(sv.Rhof);
	  else sv.RhoV = (sv.Rhof)*Vxtr;

	  if((strcmp(BC_now -> stat[BC_now -> step][5], "Extrapolated") != 0) &&
	     ((BC_now -> step >= (BC_now -> num_pts-1)) || (strcmp(BC_now -> stat[(BC_now->step)+1][5], "Extrapolated") != 0) ||
#if FLOAT_EQ_STABLE
	      EQ(strtod(BC_now->stat[BC_now->step][0], NULL), time)
#else 
	      (strtod(BC_now->stat[BC_now->step][0], NULL) == time)
#endif
	      ))
	    sv.RhoW = (pvb -> W)*(sv.Rhof);
	  else sv.RhoW = (sv.Rhof)*Wxtr;

	  if(gas_data.num_species > 1)
	    {
	      if(gas_data.eos_prod == PERFECT_GAS)
		{
		  sv.RhoE = (sv.Pres)/(R/(prod_m_frac*(gas_data.Cv_prod) + (1.0-prod_m_frac)*(gas_data.Cv_amb))) +		    
		    0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof); 
		}
	      else if(gas_data.eos_prod == JWLB)
		{
		  if(sv.rhof_prod > JWLB_RHO_CUTOFF) {

		    if((gas_data.jwlb_cutoff == 'y') && (treat_products_as_air(sv.rhof_prod,sv.Rhof) == TRUE)) 
		      sv.RhoE = (sv.Pres)/((gas_data.gam_amb)-1.0) + 0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof);
		    else {
		      T = ((sv.Pres) - (term1 + term3))/((1.0-prod_m_frac)*(sv.Rhof)*R1 + R*(gas_data.Cv_prod)/((JWLB_coeffs.v0)*v));

		      term1 = 0;
		      for(i = 0; i < 5; i++)
			if(JWLB_coeffs.R[i] > 0)
			  term1 += ((JWLB_coeffs.A[i])*(JWLB_coeffs.v0)/(JWLB_coeffs.R[i]))*exp(-(JWLB_coeffs.R[i])*v);

		      term3 = ((JWLB_coeffs.C)*(JWLB_coeffs.v0)*(1.0 - q)/(JWLB_coeffs.omega))*pow(v, -(JWLB_coeffs.omega));
		    
		      sv.RhoE = (sv.Rhof)*((prod_m_frac*(gas_data.Cv_prod) + (1.0-prod_m_frac)*(gas_data.Cv_amb))*T + 
					   prod_m_frac*(term1+term3)) + 
			0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof);
		    }
		  }
		  else sv.RhoE = (sv.Pres)/((gas_data.gam_amb)-1.0) + 0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof);
		}	      
	    }
	  else sv.RhoE = (sv.Pres)/((gas_data.gam_amb)-1.0) + 0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof);
	}
      else if(BC_now -> type == SPECIFIC) 
	{
	  if(
#if FLOAT_EQ_STABLE
	     EQ(BC_now->rhof_stat[0],EXTRAPOLATED) || EQ(BC_now->pres_stat[0],EXTRAPOLATED) || EQ(BC_now->T_stat[0],EXTRAPOLATED) ||
	     EQ(BC_now->U_stat[0],EXTRAPOLATED) || EQ(BC_now->V_stat[0],EXTRAPOLATED) || EQ(BC_now->W_stat[0],EXTRAPOLATED) 
#else
	     (BC_now->rhof_stat[0] == EXTRAPOLATED) || (BC_now->pres_stat[0] == EXTRAPOLATED) || (BC_now->T_stat[0] == EXTRAPOLATED) ||
	     (BC_now->U_stat[0] == EXTRAPOLATED) || (BC_now->V_stat[0] == EXTRAPOLATED) || (BC_now->W_stat[0] == EXTRAPOLATED)
#endif
	     )
	    {
#if 0
	      sv.Rhof = sve -> Rhof; /*'Right' state would just have extrapolated values here as no ghost cells actually used*/
	      sv.Pres = sve -> Pres;
	      sv.RhoE = sve -> RhoE;
	      sv.RhoU = sve -> RhoU;
	      sv.RhoV = sve -> RhoV;
	      sv.RhoW = sve -> RhoW;
#else
	      sv.Rhof = svc -> Rhof;
              sv.Pres = svc -> Pres;
              sv.RhoE = svc -> RhoE;
              sv.RhoU = svc -> RhoU;
              sv.RhoV = svc -> RhoV;
              sv.RhoW = svc -> RhoW;

#endif	
	      
	      Uxtr = (sv.RhoU)/(sv.Rhof); /*Extrapolated values of velocity*/
	      Vxtr = (sv.RhoV)/(sv.Rhof);
	      Wxtr = (sv.RhoW)/(sv.Rhof);
	    }
	  
	  /*Set values of state variables first - if not extrapolated already, then fixed*/

	  if(
#if FLOAT_EQ_STABLE
	     EQ(BC_now -> rhof_stat[0], FIXED)
#else
	     BC_now -> rhof_stat[0] == FIXED
#endif
	     )
	    sv.Rhof = BC_now -> rhof_stat[1];

	  if(
#if FLOAT_EQ_STABLE
	     EQ(BC_now -> pres_stat[0], FIXED)
#else
	     BC_now -> pres_stat[0] == FIXED
#endif
	     )
	    sv.Pres = BC_now -> pres_stat[1];

	  R1 = (gas_data.Cv_amb)*(gas_data.gam_amb - 1.0); 

	  if(gas_data.num_species > 1)
	    {
	      prod_m_frac = (sve -> rhof_prod)/(sve -> Rhof);
	      sv.rhof_prod = prod_m_frac*(sv.Rhof);
#if AFTERBURN
	      sv.rhof_prod_ub = ((sve->rhof_prod_ub)/(sve->rhof_prod))*(sv.rhof_prod);
#endif 

	      R = (gas_data.Cv_prod)*(gas_data.gam_prod - 1.0); /*Explosion products gas R*/
	      R = prod_m_frac*R + (1.0-prod_m_frac)*R1; /*Mixture R*/    
	    }

	  if(
#if FLOAT_EQ_STABLE
	     EQ(BC_now -> T_stat[0], FIXED)
#else
	     BC_now -> T_stat[0] == FIXED
#endif
	     ) 
	    {
	      if(gas_data.num_species > 1)
		{		  		  
		  if(gas_data.eos_prod == PERFECT_GAS) 
		    {
		      if(
#if FLOAT_EQ_STABLE
			 EQ(BC_now -> rhof_stat[0], UNKNOWN)
#else
			 BC_now -> rhof_stat[0] == UNKNOWN
#endif
			 ) 
			sv.Rhof = (sv.Pres)/(R*(BC_now->T_stat[1]));
		      else if(
#if FLOAT_EQ_STABLE
			      EQ(BC_now -> pres_stat[0], UNKNOWN)
#else
			      BC_now -> pres_stat[0] == UNKNOWN
#endif
			      )
			sv.Pres = (sv.Rhof)*R*(BC_now->T_stat[1]);
		    }
		  else if(gas_data.eos_prod == JWLB)
		    {
		      if(sv.rhof_prod > JWLB_RHO_CUTOFF) {

			if((gas_data.jwlb_cutoff == 'y') && (treat_products_as_air(sv.rhof_prod,sv.Rhof) == TRUE))  
			  sv.Pres = (sv.rhof_prod)*R1*(BC_now->T_stat[1]);
			else {
			  v = (1.0/(sv.rhof_prod))/(JWLB_coeffs.v0); 

			  /*Precompute some quantities*/
			  for(i = 0; i < 5; i++) {
			    if(JWLB_coeffs.R[i] > 0)
			      expRv[i] = exp(-(JWLB_coeffs.R[i])*v);

			    if(JWLB_coeffs.RL[i] > 0) 
			      expRLv[i] = exp(-(JWLB_coeffs.RL[i])*v);
			  }
			  q = JWLB_q(v, expRLv, expRL);

			  R = JWLB_coeffs.omega;
			  for(i = 0; i < 5; i++)
			    if(JWLB_coeffs.RL[i] > 0)
			      R += ((JWLB_coeffs.AL[i])*v + (JWLB_coeffs.BL[i]))*expRLv[i];

			  term1 = 0;
			  for(i = 0; i < 5; i++)
			    if(JWLB_coeffs.R[i] > 0)
			      term1 += (JWLB_coeffs.A[i])*expRv[i];

			  term3 = (JWLB_coeffs.C)*(1.0 - R*q/(JWLB_coeffs.omega))*pow(v,-(1.0+(JWLB_coeffs.omega)));
						
			  sv.Pres = term1 + R*(gas_data.Cv_prod)*(BC_now->T_stat[1])/((JWLB_coeffs.v0)*v) + term3;
			}
		      }
		      else sv.Pres = 0;

		      (sv.Pres) += (1.0-prod_m_frac)*(sv.Rhof)*R1*(BC_now->T_stat[1]); /*Total static pressure of mixture*/ 
		    }		  
		}
	      else
		{
		  if(BC_now -> rhof_spec == 0) /*Density was dependent*/
		    sv.Rhof = (sv.Pres)/(R1*(BC_now->T_stat[1]));
		  else if(BC_now -> pres_spec == 0)
		    sv.Pres = (sv.Rhof)*R1*(BC_now->T_stat[1]);
		}
	    } /*If temperature extrapolated, then dependent state variable automatically extrapolated too*/

	  if(gas_data.num_species > 1) {
	    sv.rhof_prod = ((sve -> rhof_prod)/(sve -> Rhof))*(sv.Rhof);
#if AFTERBURN
	    sv.rhof_prod_ub = ((sve->rhof_prod_ub)/(sve->rhof_prod))*(sv.rhof_prod);
#endif
	  }

	  if(
#if FLOAT_EQ_STABLE
	     EQ(BC_now -> U_stat[0], FIXED)
#else
	     BC_now -> U_stat[0] == FIXED
#endif
	     )
	    sv.RhoU = (BC_now -> U_stat[1])*(sv.Rhof);
	  else sv.RhoU = (sv.Rhof)*Uxtr;

	  if(
#if FLOAT_EQ_STABLE
	     EQ(BC_now -> V_stat[0], FIXED)
#else
	     BC_now -> V_stat[0] == FIXED
#endif
	     )
	    sv.RhoV = (BC_now -> V_stat[1])*(sv.Rhof);
	  else sv.RhoV = (sv.Rhof)*Vxtr;

	  if(
#if FLOAT_EQ_STABLE
	     EQ(BC_now -> W_stat[0], FIXED)
#else
	     BC_now -> W_stat[0] == FIXED
#endif
	     )
	    sv.RhoW = (BC_now -> W_stat[1])*(sv.Rhof);
	  else sv.RhoW = (sv.Rhof)*Wxtr;

#if SUPERSONIC_VORTEX
	  if(direction == WST) {
	    sv.Rhof = sv_dens(1.0, M_in, 1.0, (sve->r));
	    sv.Pres = (1.0/(gas_data.gam_amb))*pow((sv.Rhof),gas_data.gam_amb);
	    sv.RhoU = (sv.Rhof)*(M_in/(sve->r));
	    sv.RhoV = 0;
	    sv.RhoW = 0;

	    /*printf("r %g, rho %g, p %g, ru %g\n",sve->r,sv.Rhof,sv.Pres,sv.RhoU);*/
	  }
#endif

	  if(gas_data.num_species > 1)
	    {
	      if(gas_data.eos_prod == PERFECT_GAS)
		{
		  sv.RhoE = (sv.Pres)/(R/(prod_m_frac*(gas_data.Cv_prod) + (1.0-prod_m_frac)*(gas_data.Cv_amb)))+		    
		    0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof); 
		}
	      else if(gas_data.eos_prod == JWLB)
		{
		  if(sv.rhof_prod > JWLB_RHO_CUTOFF) {

		    if((gas_data.jwlb_cutoff == 'y') && (treat_products_as_air(sv.rhof_prod,sv.Rhof) == TRUE)) 
		      sv.RhoE = (sv.Pres)/((gas_data.gam_amb)-1.0) + 0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof);
		    else {
		      T = ((sv.Pres) - (term1 + term3))/((1.0-prod_m_frac)*(sv.Rhof)*R1 + R*(gas_data.Cv_prod)/((JWLB_coeffs.v0)*v));

		      term1 = 0;
		      for(i = 0; i < 5; i++)
			if(JWLB_coeffs.R[i] > 0)
			  term1 += ((JWLB_coeffs.A[i])*(JWLB_coeffs.v0)/(JWLB_coeffs.R[i]))*exp(-(JWLB_coeffs.R[i])*v);

		      term3 = ((JWLB_coeffs.C)*(JWLB_coeffs.v0)*(1.0 - q)/(JWLB_coeffs.omega))*pow(v, -(JWLB_coeffs.omega));
		    
		      sv.RhoE = (sv.Rhof)*((prod_m_frac*(gas_data.Cv_prod) + (1.0-prod_m_frac)*(gas_data.Cv_amb))*T + 
					   prod_m_frac*(term1+term3)) + 
			0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof);
		    }
		  }
		  else sv.RhoE = (sv.Pres)/((gas_data.gam_amb)-1.0) + 0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof);
		}	      
	    }
	  else sv.RhoE = (sv.Pres)/((gas_data.gam_amb)-1.0) + 0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof);
	}
      else if(BC_now -> type == NONREFLECTING)
	{
	  if(gas_data.num_species > 1)
	    soundspd = get_mixture_soundspd(sve);
	  else soundspd = SQRT((gas_data.gam_amb)*(sve->Pres)/(sve->Rhof));

	  if((direction == EST) || (direction == WST))
	    {
	      i = 0;
	      vel = (sve->RhoU)/(sve->Rhof);
	    }
	  else if((direction == NTH) || (direction == STH))
	    {
	      i = 1;
	      vel = (sve->RhoV)/(sve->Rhof);
	    }
	  else 
	    {
	      i = 2;
	      vel = (sve->RhoW)/(sve->Rhof);
	    }
	    
	  if((direction == NTH) || (direction == EST) || (direction == UPR))
	    bnorm = 1;
	  else bnorm = -1;

	  M = vel*bnorm/soundspd;
	  
	  if(M > 1) /*Supersonic flow - can just extrapolate*/
	    {
	      sv.Rhof = svc -> Rhof; 
	      sv.RhoU = svc -> RhoU;
	      sv.RhoV = svc -> RhoV;
	      sv.RhoW = svc -> RhoW;
	      sv.RhoE = svc -> RhoE;
	      sv.Pres = svc -> Pres;

	      if(gas_data.num_species > 1) {
		sv.rhof_prod = svc -> rhof_prod;
#if AFTERBURN
		sv.rhof_prod_ub = svc -> rhof_prod_ub;
#endif
	      }		
	    }
	  else if(M < -1) /*Supersonic inflow (though this really shouldn't happen) - use ambient conditions*/ 
	    {
	      sv.Rhof = IC_ambient.Rhof;
	      sv.RhoU = IC_ambient.RhoU;
	      sv.RhoV = IC_ambient.RhoV;
	      sv.RhoW = IC_ambient.RhoW;
	      sv.RhoE = IC_ambient.RhoE;
	      sv.Pres = IC_ambient.Pres;

	      if(gas_data.num_species > 1) {
		sv.rhof_prod = IC_ambient.rhof_prod;
#if AFTERBURN
		sv.rhof_prod_ub = IC_ambient.rhof_prod_ub;
#endif
	      }		  
	    }
	  else if((M >= 0) && (M <= 1)) /*Subsonic outflow - use Thompson outflow non-reflecting BCs*/
	    {
#if !NOTINCLUDE
	      /*Characteristic analysis indicates 4 flow variables can be extrapolated whilst the other is fixed by the Thompson BC.
		So we determine the pressure based on velocity gradient parallel to axis*/

	      sv.Rhof = sve -> Rhof; 
	      sv.RhoU = sve -> RhoU;
	      sv.RhoV = sve -> RhoV;
	      sv.RhoW = sve -> RhoW;
	      
	      if(gas_data.num_species > 1)
		{
		  sv.rhof_prod = sve -> rhof_prod;
		  soundspd = get_mixture_soundspd(svc);
#if AFTERBURN
		  sv.rhof_prod_ub = sve -> rhof_prod_ub;
#endif
		}
	      else soundspd = SQRT((gas_data.gam_amb)*(svc->Pres)/(svc->Rhof));

	      if(i == 0)
		dPdx = (svc->Rhof)*soundspd*(Grads -> grad_u[0])*lim;
	      else if(i == 1)
		dPdx = (svc->Rhof)*soundspd*(Grads -> grad_v[1])*lim;
	      else if(i == 2)
		dPdx = (svc->Rhof)*soundspd*(Grads -> grad_w[2])*lim;
				
	      sv.Pres = (svc->Pres) + dPdx*dx; /*The Thompson BC*/

	      if(gas_data.num_species > 1)
		sv.RhoE = get_mixture_RhoE(&sv,0);
	      else sv.RhoE = (sv.Pres)/((gas_data.gam_amb)-1.0) + 0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof);

#if 0	      
	      /*To further enforce the boundary condition, it might be good also to modify the fluid cell's pressure
		gradient and extrapolated pressure to match the ghost condition.  The left and right states are thus identical*/

	      sve -> Pres = sv.Pres;
	      sve -> RhoE = sv.RhoE;

	      /*We no longer store pressure gradient*/
	      Grads -> grad_pres[i] = dPdx;
#endif

#else
	      /*For the moment just use ambient conditions*/

	      sv.Rhof = IC_ambient.Rhof;
	      sv.RhoU = IC_ambient.RhoU;
	      sv.RhoV = IC_ambient.RhoV;
	      sv.RhoW = IC_ambient.RhoW;
	      sv.RhoE = IC_ambient.RhoE;
	      sv.Pres = IC_ambient.Pres;
	      sv.rhof_prod = IC_ambient.rhof_prod;
#if AFTERBURN
	      sv.rhof_prod_ub = IC_ambient.rhof_prod_ub;
#endif	      
#endif	      
	    }
	  else if((M >= -1) && (M < 0)) /*Subsonic inflow - use Thompson inflow non-reflecting BCs*/
	    {
#if !NOTINCLUDE
	      if(gas_data.num_species > 1)
		soundspd = get_mixture_soundspd(svc);
	      else soundspd = SQRT((gas_data.gam_amb)*(svc->Pres)/(svc->Rhof));

	      if(i == 0)
		dPdx = (svc->Rhof)*soundspd*(Grads -> grad_u[0])*lim;		  
	      else if(i == 1)
		dPdx = (svc->Rhof)*soundspd*(Grads -> grad_v[1])*lim;
	      else if(i == 2)
		dPdx = (svc->Rhof)*soundspd*(Grads -> grad_w[2])*lim;
		
	      sv.Pres = (svc->Pres) + dPdx*dx;
	      sv.Rhof = (svc->Rhof) + (dx/SQR(soundspd))*dPdx;

	      if(gas_data.num_species > 1) {
		sv.rhof_prod = ((svc->rhof_prod)/(svc->Rhof))*(sv.Rhof); /*Mass fraction constant in cell through density reconstructed*/
#if AFTERBURN
		sv.rhof_prod_ub = ((svc->rhof_prod_ub)/(svc->rhof_prod))*(sv.rhof_prod);
#endif
	      }
		
	      if(i == 0)
		{
		  sv.RhoU = (sv.Rhof)*((sve->RhoU)/(sve->Rhof)); /*Extrapolate U velocity here, but use enforced density*/
		  sv.RhoV = (sv.Rhof)*((svc->RhoV)/(svc->Rhof)); /*Traverse velocities are constant*/
		  sv.RhoW = (sv.Rhof)*((svc->RhoW)/(svc->Rhof));

		  Grads -> grad_v[0] = 0; Grads -> grad_w[0] = 0;
		}
	      else if(i == 1)
		{
		  sv.RhoU = (sv.Rhof)*((svc->RhoU)/(svc->Rhof));
		  sv.RhoV = (sv.Rhof)*((sve->RhoV)/(sve->Rhof)); /*Extrapolated V velocity*/
		  sv.RhoW = (sv.Rhof)*((svc->RhoW)/(svc->Rhof));

		  Grads -> grad_u[1] = 0; Grads -> grad_w[1] = 0;
		}
	      else if(i == 2)
		{
		  sv.RhoU = (sv.Rhof)*((svc->RhoU)/(svc->Rhof));
		  sv.RhoV = (sv.Rhof)*((svc->RhoV)/(svc->Rhof));
		  sv.RhoW = (sv.Rhof)*((sve->RhoW)/(sve->Rhof));

		  Grads -> grad_u[2] = 0; Grads -> grad_v[2] = 0;
		}

	      if(gas_data.num_species > 1)
		sv.RhoE = get_mixture_RhoE(&sv,0);
	      else sv.RhoE = sv.RhoE = (sv.Pres)/((gas_data.gam_amb)-1.0) + 0.5*(SQR(sv.RhoU) + SQR(sv.RhoV) + SQR(sv.RhoW))/(sv.Rhof);

#if 0	      
	      /*Let left and right states be identical to enforce boundary condition and modify other gradients*/
	      sve -> Pres = sv.Pres;
	      sve -> Rhof = sv.Rhof;
	      sve -> RhoU = sv.RhoU;
	      sve -> RhoV = sv.RhoV;
	      sve -> RhoW = sv.RhoW;
	      sve -> RhoE = sv.RhoE;

	      /*As pressure gradient no longer stored, this seems unnecessary*/
	      Grads -> grad_pres[i] = dPdx;
	      Grads -> grad_rhof[i] = (1.0/SQR(soundspd))*dPdx;
#endif
	      		
#else
	      sv.Rhof = IC_ambient.Rhof;
	      sv.RhoU = IC_ambient.RhoU;
	      sv.RhoV = IC_ambient.RhoV;
	      sv.RhoW = IC_ambient.RhoW;
	      sv.RhoE = IC_ambient.RhoE;
	      sv.Pres = IC_ambient.Pres;
	      sv.rhof_prod = IC_ambient.rhof_prod;
#if AFTERBURN
	      sv.rhof_prod_ub = IC_ambient.rhof_prod_ub;
#endif
#endif	      
	    }
	}
    }

  return(sv);
}

/*------------------------------------------------------------------*/

/**\brief Interpolate heat addition when it is added to the initial explosion products*/

void calc_trans_IC(double time, IC_heat *ICh, short int stage)
{
  int i, step;
  double *f, *e, *T;
  f = NULL; e = NULL; T = NULL;
  if(stage == 1) {
    f = &(ICh -> f1);
    e = &(ICh -> e1);
    T = &(ICh -> T1);
  }
  else if(stage == 2) {
    f = &(ICh -> f2);
    e = &(ICh -> e2);
    T = &(ICh -> T2);
  }

  for(i = ICh -> step; i < (ICh -> num_pts); i++) { /*Search for next written time step larger than the current time*/
    if(i != (ICh->num_pts-1)) {      
      if((ICh -> heat[i+1][0]) > time)
	break;
      else (ICh -> step)++;
    }	
  }
    
  step = ICh -> step;

  if(step >= (ICh->num_pts-1))  /*Reached end of table, assume no more heat addition*/
    {
      *f = 0;
      *e = 0;
      ICh -> num_pts = 0; /*No need to execute this routine anymore*/
    }
  else {
    /*Use linear interpolation to get corresponding values*/

    if(
#if FLOAT_EQ_STABLE
       EQ(ICh -> heat[step][0], time)
#else
       (ICh -> heat[step][0]) == time
#endif
       ) {
      *f = ICh -> heat[step][1];
      *e = ICh -> heat[step][2];
      *T = ICh -> heat[step][3];
    }
    else {
      *f = (ICh -> heat[step][1]) + ((ICh->heat[step+1][1]) - (ICh->heat[step][1]))/((ICh->heat[step+1][0]) - (ICh->heat[step][0]));
      *e = (ICh -> heat[step][2]) + ((ICh->heat[step+1][2]) - (ICh->heat[step][2]))/((ICh->heat[step+1][0]) - (ICh->heat[step][0]));
      *T = (ICh -> heat[step][3]) + ((ICh->heat[step+1][3]) - (ICh->heat[step][3]))/((ICh->heat[step+1][0]) - (ICh->heat[step][0]));
    }    
  }
}

/*------------------------------------------------------------------*/

/**\brief Interpolate flow variables for transient BCs given in a table.  All variables which aren't
 extrapolated will be converted.  Interpolation isn't possible if at t1 a variable is fixed/extrapolated and t2 it's
 extrapolated/fixed, so any time t1 > t >= t2 the variable will be extrapolated*/

void calc_trans_BC(double time, BC_data *BC_dat, short int stage) 
{
  Prim_vector *pv;
  double u, v, w, q;
  int i, step;
  pv = NULL;
  if(stage == 1)
    pv = &(BC_dat -> Pv1);
  else if(stage == 2)
    pv = &(BC_dat -> Pv2);

#if NOTINCLUDE
  printf("BC_dat -> step is %d, num_pts is %d\n",BC_dat->step,BC_dat->num_pts);
#endif

  for(i = BC_dat -> step; i < (BC_dat -> num_pts); i++) { /*Search for next written time step larger than the current time*/
    
    if(i != (BC_dat->num_pts-1)) {
      
#if NOTINCLUDE
      printf("1st col is %s, step %d, time %g\n",BC_dat->stat[i+1][0],BC_dat->step,time);
      printf("The converted time point is %g\n", strtod(BC_dat -> stat[i+1][0], NULL));
#endif
      
      if(strtod(BC_dat -> stat[i+1][0], NULL) > time)
	break;
      else (BC_dat -> step)++;
    }	
  }    
    
  step = BC_dat -> step; 
  
  if(step >= (BC_dat->num_pts-1)) { /*Reached end of table, just use these end values*/

    q = strtod(BC_dat -> stat[BC_dat->num_pts-1][1], NULL);

    if(BC_dat -> rhof_spec == 1) /*Density is in the 1st column*/
      pv -> R = q;
    else if(BC_dat -> pres_spec == 1)
      pv -> P = q;
    else if(BC_dat -> temp_spec == 1)
      pv -> T = q;

    q = strtod(BC_dat -> stat[BC_dat->num_pts-1][2], NULL);

    if(BC_dat -> rhof_spec == 2)
      pv -> R = q;
    else if(BC_dat -> pres_spec == 2)
      pv -> P = q;
    else if(BC_dat -> temp_spec == 2)
      pv -> T = q;

    pv -> U = strtod((BC_dat -> stat[BC_dat->num_pts-1][3]), NULL);
    pv -> V = strtod((BC_dat -> stat[BC_dat->num_pts-1][4]), NULL);
    pv -> W = strtod((BC_dat -> stat[BC_dat->num_pts-1][5]), NULL);
  }
  else {
    /*Use linear interpolation to get corresponding values (where variables between timesteps are indeed fixed)*/

    if(
#if FLOAT_EQ_STABLE
       EQ(strtod(BC_dat -> stat[step][0], NULL), time)
#else
       strtod(BC_dat -> stat[step][0], NULL) == time
#endif
       ) { 
      q = strtod(BC_dat -> stat[step][1], NULL);
      if(BC_dat -> rhof_spec == 1)
	pv -> R = q;
      else if(BC_dat -> pres_spec == 1)
	pv -> P = q;
      else if(BC_dat -> temp_spec == 1)
	pv -> T = q;

      q = strtod(BC_dat->stat[step][2], NULL);
      if(BC_dat -> rhof_spec == 2)
	pv -> R = q;
      else if(BC_dat -> pres_spec == 2)
	pv -> P = q;
      else if(BC_dat -> temp_spec == 2)
	pv -> T = q;

      u = strtod(BC_dat->stat[step][3], NULL);
      v = strtod(BC_dat->stat[step][4], NULL);
      w = strtod(BC_dat->stat[step][5], NULL);
    }
    else {
      q = strtod(BC_dat->stat[step][1], NULL) + 
	((strtod(BC_dat->stat[step+1][1], NULL) - strtod(BC_dat->stat[step][1], NULL))/
	 (strtod(BC_dat->stat[step+1][0], NULL) - strtod(BC_dat->stat[step][0], NULL)))*(time - strtod(BC_dat->stat[step][0], NULL));

      if(BC_dat -> rhof_spec == 1)
	pv -> R = q;
      else if(BC_dat -> pres_spec == 1)
	pv -> P = q;
      else if(BC_dat -> temp_spec == 1)
	pv -> T = q;

      q = strtod(BC_dat->stat[step][2], NULL) + 
	((strtod(BC_dat->stat[step+1][2], NULL) - strtod(BC_dat->stat[step][2], NULL))/
	 (strtod(BC_dat->stat[step+1][0], NULL) - strtod(BC_dat->stat[step][0], NULL)))*(time - strtod(BC_dat->stat[step][0], NULL));
    
      if(BC_dat -> rhof_spec == 2)
	pv -> R = q;
      else if(BC_dat -> pres_spec == 2)
	pv -> P = q;
      else if(BC_dat -> temp_spec == 2)
	pv -> T = q;

      u = strtod(BC_dat->stat[step][3], NULL) + 
	((strtod(BC_dat->stat[step+1][3], NULL) - strtod(BC_dat->stat[step][3], NULL))/
	 (strtod(BC_dat->stat[step+1][0], NULL) - strtod(BC_dat->stat[step][0], NULL)))*(time - strtod(BC_dat->stat[step][0], NULL));
    
      v = strtod(BC_dat->stat[step][4], NULL) + 
	((strtod(BC_dat->stat[step+1][4], NULL) - strtod(BC_dat->stat[step][4], NULL))/
	 (strtod(BC_dat->stat[step+1][0], NULL) - strtod(BC_dat->stat[step][0], NULL)))*(time - strtod(BC_dat->stat[step][0], NULL));

      w = strtod(BC_dat->stat[step][5], NULL) + 
	((strtod(BC_dat->stat[step+1][5], NULL) - strtod(BC_dat->stat[step][5], NULL))/
	 (strtod(BC_dat->stat[step+1][0], NULL) - strtod(BC_dat->stat[step][0], NULL)))*(time - strtod(BC_dat->stat[step][0], NULL));
    }
    
    pv -> U = u;
    pv -> V = v;
    pv -> W = w;    
  }
}

/*------------------------------------------------------------------*/

/**\brief Set flow IC values for cells which don't need further refinement*/

void set_cell_IC(Cart_cell C)
{ 
  int i, j;
  
  C->Flow_data.Grads.grad_rhof[0] = 0.0; C->Flow_data.Grads.grad_rhof[1] = 0.0; C->Flow_data.Grads.grad_rhof[2] = 0.0;
  C->Flow_data.Grads.grad_E[0] = 0.0; C->Flow_data.Grads.grad_E[1] = 0.0; C->Flow_data.Grads.grad_E[2] = 0.0;
  C->Flow_data.Grads.grad_u[0] = 0.0; C->Flow_data.Grads.grad_u[1] = 0.0; C->Flow_data.Grads.grad_u[2] = 0.0;
  C->Flow_data.Grads.grad_v[0] = 0.0; C->Flow_data.Grads.grad_v[1] = 0.0; C->Flow_data.Grads.grad_v[2] = 0.0;
  C->Flow_data.Grads.grad_w[0] = 0.0; C->Flow_data.Grads.grad_w[1] = 0.0; C->Flow_data.Grads.grad_w[2] = 0.0;

  C->Flow_data.Lim.lim = 0.0;
  C->Flow_data.Lim.rhof_lim = 0.0;
  C->Flow_data.Lim.rhof_prod_lim = 0.0;
  C->Flow_data.Lim.E_lim = 0.0;
  C->Flow_data.Lim.u_lim = 0.0;
  C->Flow_data.Lim.v_lim = 0.0;
  C->Flow_data.Lim.w_lim = 0.0;

  if(
#if FLOAT_EQ_STABLE
     EQ(C -> IC_flag, OUTSIDE) || EQ(C -> IC_flag, UNKNOWN)
#else
     (C -> IC_flag == OUTSIDE) || (C -> IC_flag == UNKNOWN)
#endif
     ) {
    C -> Flow_data.State_vec.Rhof = IC_ambient.Rhof;
    C -> Flow_data.State_vec.RhoU = IC_ambient.RhoU;
    C -> Flow_data.State_vec.RhoV = IC_ambient.RhoV;
    C -> Flow_data.State_vec.RhoW = IC_ambient.RhoW;
    C -> Flow_data.State_vec.RhoE = IC_ambient.RhoE;
    C -> Flow_data.State_vec.Pres = IC_ambient.Pres;
    C -> Flow_data.State_vec.T = IC_ambient.T;

    if(gas_data.num_species > 1) {
      C -> Flow_data.State_vec.rhof_prod = IC_ambient.rhof_prod;
#if AFTERBURN
      C -> Flow_data.State_vec.rhof_prod_ub = IC_ambient.rhof_prod_ub;
#endif
    } 
  }
  else { 
    i = (int) floor(C -> IC_flag); /*Get the IC region the cell is inside/intersected by*/
    
    C -> Flow_data.State_vec.Rhof = IC_states[i] -> Rhof; 
    C -> Flow_data.State_vec.RhoU = IC_states[i] -> RhoU;
    C -> Flow_data.State_vec.RhoV = IC_states[i] -> RhoV;
    C -> Flow_data.State_vec.RhoW = IC_states[i] -> RhoW;
    C -> Flow_data.State_vec.RhoE = IC_states[i] -> RhoE;
    C -> Flow_data.State_vec.Pres = IC_states[i] -> Pres;
    C -> Flow_data.State_vec.T = IC_states[i] -> T;

    if(gas_data.num_species > 1) {
      C -> Flow_data.State_vec.rhof_prod = IC_states[i] -> rhof_prod;
#if AFTERBURN
      C -> Flow_data.State_vec.rhof_prod_ub = IC_states[i] -> rhof_prod_ub;
#endif
    }

    if(Deton.num_deton_pts > 0) {
      C -> un_det = TRUE; /*Assume this cell is undetonated yet*/
      
      for(j = 0; j < Deton.num_deton_pts; j++) {
	if((C -> verticies[0] -> loc[0] <= Deton.dpts[j][0]) && (C -> verticies[0] -> loc[1] <= Deton.dpts[j][1]) &&
	   (C -> verticies[0] -> loc[2] <= Deton.dpts[j][2]) && (C -> verticies[7] -> loc[0] >= Deton.dpts[j][0]) &&
	   (C -> verticies[7] -> loc[1] >= Deton.dpts[j][1]) && (C -> verticies[7] -> loc[2] >= Deton.dpts[j][2])) 
	  C -> un_det = FALSE; /*Initiation detonation point lies within a cell (its bounding box) - its un_det should be FALSE*/	
      } 
    }
  }
}

/*------------------------------------------------------------------*/
