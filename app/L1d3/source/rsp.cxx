/** \file rsp.cxx
 * \ingroup l1d3
 * \brief Riemann Sub-Problem functions for L1d3
 *
 * \author DF Potter
 * \version 19-Nov-06, c++ version for cfcfd2 repository.
 */

/*-----------------------------------------------------------------*/

#include <vector>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../../../lib/util/source/useful.h"
#include "l_kernel.hh"
#include "l1d.hh"
#include "rsp.hh"

/*=================================================================*/

/* Miscellaneous functions for Riemann Sub-problem */

/// \brief Checks to see if two numbers have the same sign
int sign_check  ( double a,
		  double b )
{
    double c ;
    c=a/b;
    if (c<0) return 0;
    else return 1;
}

/// \brief Writes a flow state to the RSP structure
int fs2bs ( double * x_mid,
            double * time,
            double * u,
            struct L_flow_state * Q,
            struct ersp_solution * base_sol )
{
    base_sol->x_p = (*x_mid)+(*u)*(*time);
    base_sol->P_p = Q[0].gas->p;
    base_sol->u_p = Q[0].u;
    base_sol->T_p = Q[0].gas->T[0];
    base_sol->rho_p = Q[0].gas->rho;
    base_sol->e_p = Q[0].gas->e[0];

    return 0;
}

/*=================================================================*/

/* Firstly functions for solving the PG Riemann Problem based on the
 * method outlined by Gottlieb and Groth (1988)
 */

/// \brief Perfect gas shock jump relation based on internal velocity
///
/// Given initial state, internal velocity and direction compute the
/// shock processed state and shock speed.
double shock_wave( struct L_flow_state * QI,
                   double ustar,
                   struct pg_int_state * QIstar_pg,
	           int left )

{
    double pIstar, pIstar_d, aIstar;
    double pI, aI;
    double uI, wI, gI, cI;
    Gas_model *gmodel = get_gas_model_ptr();

    uI = QI->u;
    gI = gmodel->gamma(*(QI[0].gas));
    pI = QI[0].gas->p;
    aI = QI[0].gas->a;
    cI = gI*pI/aI;

    if (left==0) wI = (gI+1)/4*(ustar-uI)/aI+sqrt(1+pow((gI+1)/4*(ustar-uI)/aI,2));
    else wI = (gI+1)/4*(ustar-uI)/aI-sqrt(1+pow((gI+1)/4*(ustar-uI)/aI,2));

    pIstar = pI+cI*(ustar-uI)*wI;
    pIstar_d = 2*cI*pow(wI,3)/(1+pow(wI,2));
    aIstar = aI*sqrt(((gI+1)+(gI-1)*pIstar/pI)/((gI+1)+(gI-1)*pI/pIstar));

    QIstar_pg->pstar = pIstar;
    QIstar_pg->pstar_d = pIstar_d;
    QIstar_pg->astar = aIstar;

    return wI;
}

/// \brief Perfect gas rarefaction wave jump relation based on internal velocity
///
/// Given initial state, internal velocity and direction compute the
/// expanded state.
int rarefaction( struct L_flow_state * QI,
                  double ustar,
		  int left,
                  struct pg_int_state * QIstar_pg )

{
    double pIstar, pIstar_d, aIstar;
    double pI, aI, gI, uI;
    double sign=1.0;
    Gas_model *gmodel = get_gas_model_ptr();

    uI = QI[0].u;
    gI = gmodel->gamma(*(QI[0].gas));
    pI = QI[0].gas->p;
    aI = QI[0].gas->a;

    if (left==0)  sign=1.0;
    else if (left==1)  sign=-1.0;
    aIstar = aI+sign*(gI-1)/2*(ustar-uI);
    pIstar = pI*pow(aIstar/aI,2*gI/(gI-1));
    pIstar_d = sign*gI*pIstar/aIstar;

    QIstar_pg->pstar = pIstar;
    QIstar_pg->pstar_d = pIstar_d;
    QIstar_pg->astar = aIstar;
    return 0;
}

/// \brief Gottlieb and Groth's efficient PG Riemann Solver
///
/// Firstly the wave pattern is determined a-priori, then the internal
/// states found by iterating velocity using Newton's method.
int perfect_gas ( struct L_flow_state * QL,
                  struct L_flow_state * QR,
                  struct pg_int_state * S_pg )

{
    double WSL, WSR;
    double uR, pR, aR, gR;
    double uL, pL, aL, gL;
    double u_scn, u_ncr, u_ncs, u_rcn, u_rcvr;
    double uLbar, uRbar, sigma, z, u0star, epsilon, ustar=0.0, pstar;
    double vR=0.0, vL=0.0;
    int designation=3, count;
    double RL, RR, T_starL, rho_starL, T_starR, rho_starR;
    double tol=1.0e-10;
    struct pg_int_state QLstar_pg[1], QRstar_pg[1];
    Gas_model *gmodel = get_gas_model_ptr();

    uL = QL[0].u; pL = QL[0].gas->p;
    aL = QL[0].gas->a; gL = gmodel->gamma(*(QL[0].gas));

    uR = QR[0].u; pR = QR[0].gas->p;
    aR = QR[0].gas->a; gR = gmodel->gamma(*(QR[0].gas));

    /*
     * Determine wave pattern a-priori
    */

    u_rcvr = uL+2*aL/(gL-1)+2*aR/(gR-1);

    if (pR > pL)  {
        sigma = gR;
        u_scn = uL-(aL/gL)*(pR/pL-1)/sqrt((gL+1)/(2*gL)*pR/pL+(gL-1)/(2*gL));
        u_ncr = uL+2*aR/(gR-1)*(1-pow(pL/pR,(gR-1)/(2*gR)));
        if (uR < u_scn)  designation = 0;             /* SCS */
        else if (uR == u_scn)  designation = 99;      /* SCN */
        else if (uR < u_ncr)  designation = 1;        /* SCR */
        else if (uR == u_ncr)  designation = 98;      /* NCR */
        else if (uR < u_rcvr)  designation = 2;       /* RCR */
        else if (uR == u_rcvr) designation = 2;      /* RCVR (NA) */
        else designation = 2;                        /* RCVCR (NA) */
    }
    else  {  /* pR <= pL */
        sigma = gL;
        u_ncs = uL-(aR/gR)*(pL/pR-1)/sqrt((gR+1)/(2*gR)*pL/pR+(gR-1)/(2*gR));
        u_rcn = uL+2*aL/(gL-1)*(1-pow(pR/pL,(gL-1)/(2*gL)));
        if (uR < u_ncs)  designation = 0;             /* SCS */
        else if (uR == u_ncs)  designation = 96;      /* NCS */
        else if (uR < u_rcn)   designation = 3;       /* RCS */
        else if (uR == u_rcn)  designation = 95;      /* RCN */
        else if (uR < u_rcvr)  designation = 2;       /* RCR */
        else if (uR == u_rcvr) designation = 2;        /* RCVR (NA) */
        else designation = 2;                          /* RCVCR (NA) */ 
    }
    /*
     * Initial guesses from Gottlieb and Groth (1994)
     */

    uLbar = uL+2.0/(gL-1.0)*aL;
    uRbar = uR-2.0/(gR-1.0)*aR;

    z=(gL-1)/(gR-1)*aR/aL*pow(pL/pR,(sigma-1)/(2*sigma));
    u0star = (uLbar*z+uRbar)/(1+z);

    /*
     * Unique iterative loop required for each wave pattern
     */

    count=1;
    epsilon=1;
    ustar=u0star;
    while (epsilon > tol)  {
        if ((designation==0)||(designation==1)||(designation==99))  {
            WSL = shock_wave( QL, ustar, QLstar_pg, 1 );
            vL=uL+aL*WSL;
        }
        if ((designation==96)||(designation==3)||(designation==0))  {
            WSR = shock_wave( QR, ustar, QRstar_pg, 0 );
            vR = uR+aR*WSR;
        }
        if ((designation==2)||(designation==3)||(designation==95))
            rarefaction( QL, ustar, 1, QLstar_pg );
        if ((designation==1)||(designation==2)||(designation==98))
            rarefaction( QR, ustar, 0, QRstar_pg );
        if ((designation==96)||(designation==98))  {
            QLstar_pg->pstar=QL[0].gas->p;
            QLstar_pg->pstar_d=0;
        }
        if ((designation==95)||(designation==99))  {
            QRstar_pg->pstar=QR[0].gas->p;
            QRstar_pg->pstar_d=0;
        }
        ustar -= (QLstar_pg->pstar-QRstar_pg->pstar)/
                 (QLstar_pg->pstar_d-QRstar_pg->pstar_d);
        epsilon=fabs(1.0-QLstar_pg->pstar/QRstar_pg->pstar);
        count++;
    }

    pstar=(QLstar_pg->pstar+QRstar_pg->pstar)/2.0;

    gL=gmodel->gamma(*(QL[0].gas)); RL=gmodel->R(*(QL[0].gas));
    gR=gmodel->gamma(*(QR[0].gas)); RR=gmodel->R(*(QR[0].gas));

    T_starL=pow(QLstar_pg->astar,2)/(gL*RL);
    rho_starL=(QLstar_pg->pstar)/(RL*T_starL);
    T_starR=pow(QRstar_pg->astar,2)/(gR*RR);
    rho_starR=(QRstar_pg->pstar)/(RR*T_starR);

    S_pg[0].pstar=pstar;
    S_pg[0].v_rwt = QL[0].u - QL[0].gas->a;
    S_pg[0].v_cs = ustar;
    S_pg[0].v_sw = vR;

    return 0;
}

/*=================================================================*/

/* Now an exact Riemann solver for an arbitrary EOS */

/// \brief Zero value expression for expected energy change across a shock
///
/// Calculates energy change based on Rankine-Hugoniot relations and 
/// computes normalised deviation from that expected from the thermodynamic
/// definitions.
double shock_energy ( double Us,
                      double v_i,
                      double P_i,
                      double rho_i,
                      double P_star,
                      struct L_flow_state * Qi )

{

    double V_i, V_star;
    double e, p, rho;
    double e_star, rho_star, delta_h, zero_val;
    L_flow_state Q[1];
    Gas_model *gmodel = get_gas_model_ptr();
    Q[0].gas = new Gas_data(gmodel);

    Q[0].gas->copy_values_from(*(Qi[0].gas));           /* makes mass fractions correct */
    V_i=v_i-Us;
    V_star=V_i/fabs(V_i)*fabs(((P_i-P_star)+rho_i*pow(V_i,2))/(rho_i*V_i));
    rho_star=rho_i*V_i/V_star;
    e=Qi[0].gas->e[0]; rho=Qi[0].gas->rho; p=Qi[0].gas->p;

    Q[0].gas->rho = rho_star;
    Q[0].gas->p = P_star;
    gmodel->eval_thermo_state_rhop(*(Q[0].gas));
    gmodel->eval_transport_coefficients(*(Q[0].gas));
    e_star=Q[0].gas->e[0];

    delta_h = - (e + p/rho) + ( e_star + P_star/rho_star ) ;
    zero_val=(0.5*(pow(V_star,2)-pow(V_i,2))+delta_h)/(delta_h);

    delete Q[0].gas;
    return zero_val;
}

/// \brief Arbitrary EOS shock jump routine
///
/// Iterates shock speed to solve the zero-value expression contained
/// in shock_energy.  The shock processed state and shock speed are 
/// returned.
double shock_jump ( double P_star,
                 struct L_flow_state * QI,
                 int shock_right,
		 struct L_flow_state * QIstar )
{
    double rhoI, TI, uI, PI, gI, C_vI;
    double dT_dpI, dT_drhoI;
    double Us_PG, Us_APRX, Us_n, Us_nm1, Us_np1, Us_I;
    double VI, V_starI, rho_starI, u_starI;
    double errors, tol=1.0e-10;
    double se_n, se_nm1, abs_err, first_err;
    int count, sign, ig_shift_flag=0;
    int status;
    L_flow_state Q1[1];
    Gas_model *gmodel = get_gas_model_ptr();
    Q1[0].gas = new Gas_data(gmodel);

    /* usually will be a right facing shock */
    if (shock_right==0) sign=-1;
    else sign=1;

    rhoI=QI[0].gas->rho; TI=QI[0].gas->T[0]; uI=QI[0].u;
    PI=QI[0].gas->p; gI=gmodel->gamma(*(QI[0].gas));
    C_vI=gmodel->Cv(*(QI[0].gas));
    Q1[0].gas=QI[0].gas;
    dT_drhoI = gmodel->dTdrho_const_p(*(Q1[0].gas), status);
    Q1[0].gas=QI[0].gas;
    dT_dpI = gmodel->dTdp_const_rho(*(Q1[0].gas), status);

    Us_PG=sign*sqrt((PI/(2.0*rhoI))*(P_star/PI*(gI+1)+gI-1))+uI;
    Us_APRX=Us_PG*1.02;
    count=0;
    errors=shock_energy(Us_PG,uI,PI,rhoI,P_star,QI);
    abs_err=fabs(errors);
    first_err=abs_err;

    Us_n=Us_PG;
    Us_nm1=Us_APRX;
    Us_np1=Us_PG;

    while (abs_err>tol)  {
        /* a check for divergence */
        if (count==4 && abs_err>first_err)  {
            ig_shift_flag++;
            if (ig_shift_flag==500)  {
                printf("WARNING: bad initial guess in shock_jump\n");
                Us_np1 = Us_PG;
                goto shortcut;
            }
            else  {
                count=0;
                Us_nm1=Us_PG*pow(1.001,pow(-1.0,ig_shift_flag)*ig_shift_flag);
                Us_n=Us_PG;
            }
        }
        else if (count!=0)  {
            Us_nm1=Us_n;
            Us_n=Us_np1;
        }
        se_n=shock_energy(Us_n,uI,PI,rhoI,P_star,QI);
        se_nm1=shock_energy(Us_nm1,uI,PI,rhoI,P_star,QI);
        Us_np1=Us_n-se_n*(Us_n-Us_nm1)/(se_n-se_nm1);
        errors=shock_energy(Us_np1,uI,PI,rhoI,P_star,QI);
        abs_err=fabs(errors);
        count++;
    }

    shortcut:

    Us_I=Us_np1;
    VI=uI-Us_np1;
    V_starI=VI/fabs(VI)*fabs(((PI-P_star)+rhoI*pow(VI,2))/(rhoI*VI)); 
    u_starI=V_starI+Us_I;
    rho_starI=rhoI*VI/V_starI;

    QIstar[0].gas=QI[0].gas;
    QIstar[0].u=u_starI;
    QIstar[0].gas->rho=rho_starI;
    QIstar[0].gas->p=P_star;
    gmodel->eval_thermo_state_rhop(*(QIstar[0].gas));
    gmodel->eval_transport_coefficients(*(QIstar[0].gas));

    delete Q1[0].gas;
    return Us_I;
}

/// \brief Arbitrary EOS rarefaction wave jump routine
///
/// Isentropically expands an initial state to a final pressure and computes
/// the expanded velocity using the Riemann invariant for an unsteady expansion.
int rare_jump ( double P_star,
                struct L_flow_state * QI,
                int rare_left,
	        struct L_flow_state * QIstar)
{
    double u_starI=QI[0].u, rho_scale=1.0e-3, rho_f=0.0;
    double T=0.0, rho=0.0, C_v, dT_dp, dT_drho_s;
    double f, T_ib, dT_drho_ib, dT_dp_ib,df_drho, f_new, d_rho;
    int run=1, sign, dens_too_high=1;
    int status;
    double fn_abs, f_abs, rho_limit=0.1, drho_sign=-1.0;
    L_flow_state Q1[1], Q[1];
    Gas_model *gmodel = get_gas_model_ptr();
    Q1[0].gas = new Gas_data(gmodel);
    Q[0].gas = new Gas_data(gmodel);

    if (rare_left==0) sign=1;
    else sign=-1;

    QIstar[0].gas->copy_values_from(*(QI[0].gas));
    Q1[0].gas->copy_values_from(*(QI[0].gas));
    Q[0].gas->copy_values_from(*(QI[0].gas));               /* Ensures correct mass fraction */

    /* Step along the isentrope until P_star is reached */
    Q[0].gas->T[0]=QI[0].gas->T[0]; Q[0].gas->rho=QI[0].gas->rho; C_v=gmodel->Cv(*(QI[0].gas));
    rho=QI[0].gas->rho; T=QI[0].gas->T[0]; u_starI = QI[0].u;
    gmodel->eval_thermo_state_rhoT(*(Q[0].gas));
    gmodel->eval_transport_coefficients(*(Q[0].gas));
    while (run==1)  {
        d_rho=drho_sign*rho*rho_scale;
        /* Isentrope diffs */
        C_v=gmodel->Cv(*(Q[0].gas));
        Q1[0].gas=Q[0].gas;
        dT_dp = gmodel->dTdp_const_rho(*(Q1[0].gas), status);
        dT_drho_s = T/(C_v*pow(rho,2)*dT_dp);
        /* Calculate step for isobar only if a reasonable density exists */
        if (dens_too_high==0)  {
            Q[0].gas->rho=rho; Q[0].gas->p=P_star;
            gmodel->eval_thermo_state_rhop(*(Q[0].gas)); T_ib=Q[0].gas->T[0];
            Q1[0].gas->copy_values_from(*(Q[0].gas)); dT_drho_ib = gmodel->dTdrho_const_p(*(Q1[0].gas), status);
            Q1[0].gas->copy_values_from(*(Q[0].gas)); dT_dp_ib = gmodel->dTdp_const_rho(*(Q1[0].gas), status);
            f=T-T_ib;
            df_drho=dT_drho_s-dT_drho_ib;               /* combine diffs */
            T+=dT_drho_s*d_rho;                         /* update both temperatures */ 
            T_ib+=dT_drho_ib*d_rho;
            rho+=d_rho;                                 /* step density forward */
            f_new=T-T_ib;
            if (sign_check(f,f_new)==0)  {
                rho_f=(rho-d_rho)-f/df_drho;            /* linear interpolation */
                T-=dT_drho_s*(rho_f-(rho-d_rho));
                d_rho = rho_f - (rho-d_rho);	        /* Fix d_rho for integration */ 
                run=0;                                  /* terminate run */
            } else {
                /* check for right direction if stepping is to continue */
                fn_abs=fabs(f_new);
                f_abs=fabs(f);
                if (fn_abs>f_abs) drho_sign*=-1.0;      /* reverse density step */
            }
        }
        else if (dens_too_high==1)  {
            T+=dT_drho_s*d_rho; rho+=d_rho;
            /* so that bad temperature warnings do not occur */
            rho_limit=P_star/(gmodel->R(*(QI[0].gas))*200.0);	
            if (rho<rho_limit)  dens_too_high=0;
        }
        /* Calculate expanded velocity */
        Q[0].gas->rho=rho; Q[0].gas->T[0]=T;
        gmodel->eval_thermo_state_rhoT(*(Q[0].gas));
        u_starI += sign*(Q[0].gas->a / Q[0].gas->rho) * d_rho;
    }
    QIstar[0].gas->rho=rho_f;
    QIstar[0].gas->T[0]=T;
    QIstar[0].u=u_starI;
    QIstar[0].gas->p=P_star;
    gmodel->eval_thermo_state_rhop(*(QIstar[0].gas));
    gmodel->eval_transport_coefficients(*(QIstar[0].gas));

    delete Q1[0].gas;
    delete Q[0].gas;
    return 0;
}

/// \brief Arb. EOS Riemann zero-value expression for internal velocity
///
/// Given initial LHS & RHS states and a guess for internal pressure
/// calculate the error in internal velocity across the CS
double f_zero ( double P_star,
                struct L_flow_state * QL,
                struct L_flow_state * QR )
{
    double PL, PR, uL, uR;
    double u_starL, u_starR;
    double zero_val;
    L_flow_state QLstar[1], QRstar[1];
    Gas_model *gmodel = get_gas_model_ptr();
    QLstar[0].gas = new Gas_data(gmodel);
    QRstar[0].gas = new Gas_data(gmodel);

    PL=QL[0].gas->p; uL=QL[0].u;
    PR=QR[0].gas->p; uR=QR[0].u;

/*
 * Determine wave pattern for current internal pressure guess
 */

    if (P_star/PR>1)  {
        shock_jump(P_star, QR, 1, QRstar);
        u_starR=QRstar[0].u;
    }
    else if (P_star==PR) u_starR=uR;
    else  {
        rare_jump(P_star, QR, 0, QRstar);
        u_starR=QRstar[0].u;
    }

    if (P_star/PL<1)  {
        rare_jump(P_star, QL, 1, QLstar);
        u_starL=QLstar[0].u;
    }
    else if (P_star==PL) u_starL=uL;
    else  {
        shock_jump(P_star, QL, 0, QLstar);
        u_starL=QLstar[0].u;
    }

    zero_val=u_starL/(u_starR)-1.0;

    delete QLstar[0].gas;
    delete QRstar[0].gas;
    return zero_val;
}

/// \brief Arb. EOS Riemann calculator
///
/// Same as f_zero but is used for a correct P_star to find the internal
/// flow states.
int final_states ( double P_star,
                   struct L_flow_state * QL,
                   struct L_flow_state * QR,
                   struct L_flow_state * QLstar,
                   struct L_flow_state * QRstar,
                   double * Us_L,
                   double * Us_R )
{
    double PL, PR, uL, uR;
    double u_starL, u_starR;

    PL=QL[0].gas->p; uL=QL[0].u;
    PR=QR[0].gas->p; uR=QR[0].u;

    /*
     * Determine wave pattern for current internal pressure guess
     */

    if (P_star/PR>1)  {
        *Us_R=shock_jump(P_star, QR, 1, QRstar);
        u_starR=QRstar[0].u;
    }
    else if (P_star==PR) u_starR=uR;
    else  {
        rare_jump(P_star, QR, 0, QRstar);
        u_starR=QRstar[0].u;
    }

    if (P_star/PL<1)  {
        rare_jump(P_star, QL, 1, QLstar);
        u_starL=QLstar[0].u;
    }
    else if (P_star==PL)  {
        u_starL=uL;
    }
    else  {
        *Us_L=shock_jump(P_star, QL, 0, QLstar);
        u_starL=QLstar[0].u;
    }

    return 0;
}

/// \brief Steps AND records flow state through an unsteady expansion
///
/// Same as rare_jump but records flow state at equi-spaced pressure intervals
/// in order to describe complete wave pattern profile.
int step_through_rw( struct L_flow_state * QI,
		     struct L_flow_state * QF,
		     struct ersp_solution * base_sol,
		     int * count,
		     int rare_left,
		     double x_mid,
		     double time)
{
    double delta_p;
    double rho_scale=1.0e-3, rho_f=0.0, u;
    double T=0.0, rho=0.0, C_v, dT_dp, dT_drho_s;
    double f, T_ib, dT_drho_ib, dT_dp_ib,df_drho, f_new;
    double d_rho;
    int i=1, run=1, sign;
    double fn_abs, f_abs;
    int dens_too_high=1;
    double rho_limit=0.1;
    double drho_sign=-1.0;
    int status;
    L_flow_state Q1[1], Q[1];
    Gas_model *gmodel = get_gas_model_ptr();
    Q1[0].gas = new Gas_data(gmodel);
    Q[0].gas = new Gas_data(gmodel);

    /* First find the desired pressure decrement */
    delta_p=(QI[0].gas->p - QF[0].gas->p)/ERSP_NX;

    /* Now step through the expansion fan as usual
     * but append the solution to base_sol every time
     * p < p_last + delta_p
     */

    if (rare_left==0) sign=1;
    else sign=-1;

    Q1[0].gas=QI[0].gas;
    Q[0].gas = QI[0].gas ;

    Q[0].gas->T[0]=QI[0].gas->T[0]; Q[0].gas->rho=QI[0].gas->rho; C_v=gmodel->Cv(*(QI[0].gas));
    rho=QI[0].gas->rho; T=QI[0].gas->T[0]; Q[0].u = QI[0].u;
    gmodel->eval_thermo_state_rhoT(*(Q[0].gas));
    u = QI[0].u - QI[0].gas->a;
    fs2bs(&x_mid,&u,&time,QI,&(base_sol[*count]));
    while (run==1)  {
        d_rho=drho_sign*rho*rho_scale;
        C_v=gmodel->Cv(*(Q[0].gas));
        Q1[0].gas->copy_values_from(*(Q[0].gas));
        dT_dp = gmodel->dTdp_const_rho(*(Q1[0].gas), status);
        dT_drho_s = T/(C_v*pow(rho,2)*dT_dp);
        if (dens_too_high==0)  {
            Q[0].gas->rho=rho; Q[0].gas->p=QF[0].gas->p;
            gmodel->eval_thermo_state_rhop(*(Q[0].gas)); T_ib=Q[0].gas->T[0];
            Q1[0].gas->copy_values_from(*(Q[0].gas)); dT_drho_ib = gmodel->dTdrho_const_p(*(Q1[0].gas), status);
            Q1[0].gas->copy_values_from(*(Q[0].gas)); dT_dp_ib = gmodel->dTdp_const_rho(*(Q1[0].gas), status);
            f=T-T_ib; df_drho=dT_drho_s-dT_drho_ib; 
            T+=dT_drho_s*d_rho; T_ib+=dT_drho_ib*d_rho;
            rho+=d_rho;
            f_new=T-T_ib;
            if (sign_check(f,f_new)==0)  {
                rho_f=(rho-d_rho)-f/df_drho;
                T-=dT_drho_s*(rho_f-(rho-d_rho));
                d_rho = rho_f - (rho-d_rho);	
                run=0; 
            } else {
                fn_abs=fabs(f_new);
                f_abs=fabs(f);
                if (fn_abs>f_abs) drho_sign*=-1.0;
            }
        }
        else if (dens_too_high==1)  {
            T+=dT_drho_s*d_rho;
            rho+=d_rho;
            /* so that bad temperature warnings do not occur */
            rho_limit=QF[0].gas->p/(gmodel->R(*(QI[0].gas))*200.0);
            if (rho<rho_limit) dens_too_high=0;
        }
        /* Calculate expanded velocity */
        Q[0].gas->rho=rho; Q[0].gas->T[0]=T;
        gmodel->eval_thermo_state_rhoT(*(Q[0].gas));
        Q[0].u += sign*(Q[0].gas->a / Q[0].gas->rho) * d_rho;
        if (Q[0].gas->p<=(QI[0].gas->p-i*delta_p))   {
            i++; (*count)++; u = Q[0].u - Q[0].gas->a;
            fs2bs(&x_mid,&u,&time,Q,&(base_sol[*count]));
        }
    }

    delete Q[0].gas;
    delete Q1[0].gas;
    return 0;
}

/// \brief Main routine for computing a solution to the exact Riemann problem 
///
/// Given initial LHS & RHS states solves the exact Riemann problem for an 
/// arbitrary EOS.  Three fundamental velocities (RWT, CS and SW) are returned
/// along with the option of returning the whole discretised profile as a structure.
/// \param QL point to LHS flow state
/// \param QR point to RHS flow state
/// \param left_vel Calculated rarefaction wave tail velocity
/// \param cs_vel Calculated contact surface velocity
/// \param right_vel Calculated shock-wave velocity
/// \param base_sol Pointer to solution structure
/// \param option A switch (see below)
int exact_riemann( struct L_flow_state * QL,
                   struct L_flow_state * QR,
		   double * left_vel,
		   double * cs_vel,
		   double * right_vel,
		   struct ersp_solution * base_sol,
		   int option )
{
    double P_star_0, P_star_1, P_star_n=0.0, P_star_nm1, P_star_np1=0.0;
    double tol=1.0e-6, errors, abs_err;
    int count=0, i, extra_points;
    double Us_L[1],Us_R[1];
    double x_mid = 0.0, time =1.0, u=0.0, factor;
    L_flow_state QLstar[1], QRstar[1];
    struct pg_int_state S_pg[1]; 
    Gas_model *gmodel = get_gas_model_ptr();
    QLstar[0].gas = new Gas_data(gmodel);
    QRstar[0].gas = new Gas_data(gmodel);

    /*****************************************************************************
     * This is exact_riemann version 2.0 by Daniel F. Potter, November 2006   
     *
     * Option:
     * 0 = LUT EOS, don't write base_sol
     * 1 = LUT EOS, write base_sol
     * 2 = PG EOS, don't write base_sol
     *
     ******************************************************************************/

    gmodel->eval_thermo_state_pT(*(QL[0].gas));
    gmodel->eval_thermo_state_pT(*(QR[0].gas));

    /*****************************************************************************
     * Perfect gas solution
     *****************************************************************************/

    perfect_gas( QL, QR, S_pg );
    if (option==2) goto pg_only;

    /*****************************************************************************
     * Arbitrary EOS solution
     *****************************************************************************/

    P_star_0 = 0.95*S_pg[0].pstar;
    P_star_1 = 0.98*S_pg[0].pstar;
    P_star_np1=P_star_1;
    errors=f_zero(P_star_1,QL,QR);
    abs_err=fabs(errors);

    # if (DEBUG >= 1)
    printf("The initial guess for the real solver has an error of %g \n",errors);
    # endif

    /* Iterate for internal pressure using velocity as iterative variable */
    while (abs_err>tol)  {
        if (count==0)  {
            P_star_n=P_star_1; P_star_nm1=P_star_0;
        } else  {
            P_star_nm1=P_star_n; P_star_n=P_star_np1;
        }
        P_star_np1=P_star_n-f_zero(P_star_n,QL,QR)*(P_star_n-P_star_nm1)/
                   (f_zero(P_star_n,QL,QR)-f_zero(P_star_nm1,QL,QR));
        abs_err=fabs(f_zero(P_star_np1,QL,QR));
        count++;
        # if (DEBUG >= 1)
        printf("at iteration %d the abs_err is %g\n",count,abs_err);
        # endif
        if (count==10) break;
    }

    final_states(P_star_np1,QL,QR,QLstar,QRstar,Us_L,Us_R);
    printf("************* Exact Riemann solution found ***************\n\
* Iterations = %d, Relative error in velocity = %0.3e *\n",
           count,abs_err);
    printf("* Us = %0.1f m/s; P* = %0.1f kPa; u* = %0.1f m/s\n",
           *Us_R, 0.001*P_star_np1, QLstar[0].u);
    printf("**********************************************************\n");

    *cs_vel=(QRstar[0].u+QLstar[0].u)/2;  /* Average to be pedantic */
    *left_vel=(QL[0].u - QL[0].gas->a);
    *right_vel=*Us_R;

    pg_only:

    if (option==2)   {
        *cs_vel=S_pg[0].v_cs;
        *left_vel=S_pg[0].v_rwt;
        *right_vel=S_pg[0].v_sw;
    } else if (option==1)  {
        /* Now create property profile for user specified range and time */
        count=0;

        if (QL[0].gas->p > P_star_np1)  {
            /* Left slug state */
            step_through_rw(QL,QLstar,base_sol,&count,1,x_mid,time);
            /* Bottom of expansion fan */
            QLstar[0].gas->p=P_star_np1;
            final_states(QLstar[0].gas->p,QL,QR,QLstar,QRstar,Us_L,Us_R);
            ++count;
            u=QLstar[0].u - QLstar[0].gas->a;
            fs2bs(&x_mid,&u,&time,QLstar,&(base_sol[count]));

            /* Some more points to fatten the test slug region */
            extra_points = 5;
            for (i=0; i<extra_points; i++)  {
                ++count;
                factor = (5-i)/6.0;
                u=QLstar[0].u - factor * QLstar[0].gas->a;
                fs2bs(&x_mid,&u,&time,QLstar,&(base_sol[count]));
            }
            /* top of contact surface */
            ++count; u=QLstar[0].u;
            fs2bs(&x_mid,&u,&time,QLstar,&(base_sol[count]));
        }

        if (QL[0].gas->p < P_star_np1)  {
            /* Bottom of shock */
            ++count; u=*Us_L;
            fs2bs(&x_mid,&u,&time,QL,&(base_sol[count]));

            /* top of the shock */
            ++count; u=*Us_L;
            fs2bs(&x_mid,&u,&time,QLstar,&(base_sol[count]));

            /* bottom of contact surface */
            ++count; u=QLstar[0].u;
            fs2bs(&x_mid,&u,&time,QLstar,&(base_sol[count]));
        }

        if (QR[0].gas->p > P_star_np1)  {
            /* Bottom of contact surface */
            ++count; u = QRstar[0].u;
            fs2bs(&x_mid,&u,&time,QRstar,&(base_sol[count]));

            /* Bottom of expansion fan */
            ++count; u=QRstar[0].u + QRstar[0].gas->a;
            fs2bs(&x_mid,&u,&time,QRstar,&(base_sol[count]));
            base_sol[count].x_p = x_mid+(QRstar[0].u + QRstar[0].gas->a)*time;

            /* Expansion fan right */
            ++count;
            step_through_rw(QL,QLstar,base_sol,&count,0,x_mid,time);
        }

        if (QR[0].gas->p < P_star_np1)  {
            /* Top of contact surface */
            ++count; u = QRstar[0].u;
            fs2bs(&x_mid,&u,&time,QRstar,&(base_sol[count]));

            /* top of the shock */
            ++count; u=*Us_R;
            fs2bs(&x_mid,&u,&time,QRstar,&(base_sol[count]));

            /* bottom of the shock */
            ++count; u=*Us_R;
            fs2bs(&x_mid,&u,&time,QR,&(base_sol[count]));
        }
    }

    delete QRstar[0].gas;
    delete QLstar[0].gas;
    return 0;
}

/*=================================================================*/

/* Now functions used by L1d2 to solve the RSP at a diaphragm */

/// \brief Allocate memory for RSD structure
/// \param RSD : Pointer to Riemann simulation data structure
int RSD_alloc( struct riemann_simulation_data *RSD )
{
    Gas_model *gmodel = get_gas_model_ptr();
    int i, tb=0;

    RSD->rs_lcells = (struct L_cell *) calloc(RIEMANN_CELLS,sizeof(struct L_cell));
    tb += RIEMANN_CELLS * sizeof(struct L_cell);
    RSD->rs_rcells = (struct L_cell *) calloc(RIEMANN_CELLS,sizeof(struct L_cell));
    tb += RIEMANN_CELLS * sizeof(struct L_cell);
    for (i=0; i<RIEMANN_CELLS; i++) {
        RSD->rs_lcells[i].gas = new Gas_data(gmodel);
        RSD->rs_rcells[i].gas = new Gas_data(gmodel);
    }
    RSD->base_sol = (struct ersp_solution *) calloc(ERSP_NX+20,sizeof(struct ersp_solution));
    tb += RIEMANN_CELLS * sizeof(struct ersp_solution);
    RSD->QL[0].gas = new Gas_data(gmodel);
    RSD->QR[0].gas = new Gas_data(gmodel);
    // FIX-ME PJ
    // byte counting and memory leak.
    printf("RSD_alloc(): allocated %d bytes of memory.\n", tb);
    return SUCCESS;
} /* end function RSD_alloc */

/// \brief Interpolates the flow properties in a discretised wave-pattern
///
/// The Riemann wave pattern is described by a structure of discretised points.
/// This function interpolates between these points to solve for intermediate states.
/// \param x : Axial location of interest
/// \param x_mid : Wave pattern origin
/// \param delta_t : Time elapsed since rupture
/// \param rho : Interpolated density
/// \param e : Interpolated energy
/// \param u : Interpolated velocity
int riemann_interpolate( double x,
			 double x_mid,
			 double delta_t,
			 struct ersp_solution * base_sol,
			 double * rho,
			 double * e,
			 double * u)
{

    int i, i_current, i_last=0;
    double x_current, dx, drho_dx, de_dx, du_dx, delta_x, x_last;

    /* Note there are always ERSP_NX+10 entries in the base_sol structure */

    for (i=0; i<ERSP_NX+10; i++)  {
        i_current=i;
        x_current=base_sol[i].x_p*delta_t+x_mid;

        if (x_current>x) break;
        x_last=x_current;
        i_last=i;
    }

    if (i_current==0)  {
        /* If the point is before the tail of the expansion wave
           force the solution to the LHS state */
        *rho = base_sol[0].rho_p;
        *e = base_sol[0].e_p;
        *u = base_sol[0].u_p;
    }
    else if (i_current==ERSP_NX+9 && i_last==ERSP_NX+9)  {
        /* The point is ahead of the shockfront, so force the solution to 
           the RHS state */
        *rho = base_sol[ERSP_NX+9].rho_p;
        *e = base_sol[ERSP_NX+9].e_p;
        *u = base_sol[ERSP_NX+9].u_p;
    } else {
        dx=(base_sol[i_current].x_p - base_sol[i_last].x_p)*delta_t;
        delta_x = x - (base_sol[i_last].x_p*delta_t+x_mid);
        drho_dx=(base_sol[i_current].rho_p - base_sol[i_last].rho_p)/dx;
        de_dx=(base_sol[i_current].e_p - base_sol[i_last].e_p)/dx;
        du_dx=(base_sol[i_current].u_p - base_sol[i_last].u_p)/dx;

        *rho = base_sol[i_last].rho_p + drho_dx * delta_x;
        *e = base_sol[i_last].e_p + de_dx * delta_x;
        *u = base_sol[i_last].u_p + du_dx * delta_x;
    }

return 0;
}

/// \brief Find the location of a cells right face according to RSP profile
/// \param x_left : left face axial location
/// \param x_diaphragm : pointer to all slug data
/// \param delta_t : time lapse since simulated rupture
/// \param m_cell : the cell under investigation
/// \param base_sol : the rsp base solution structure
/// \param cell_sol : the interpolated properties for the new cell centroid
/// \returns : success or failure
int find_x_pos_right( double x_left, 
                      double x_diaphragm, 
                      double delta_t, 
                      struct L_cell * m_cell, 
                      struct ersp_solution * base_sol,
                      struct ersp_solution * cell_sol)
{
    double x_right_guess[3], x_mid_guess[3], fz[3];
    double rho,e,u;
    struct L_flow_state QD[1];
    Gas_model *gmodel = get_gas_model_ptr();
    QD[0].gas = new Gas_data(gmodel);

    riemann_interpolate(x_left,x_diaphragm,delta_t,base_sol,&rho,&e,&u);
    x_right_guess[1]=(m_cell[0].mass + rho*m_cell[0].area*x_left)/(rho*m_cell[0].area);
    x_mid_guess[1] = 0.5*(x_right_guess[1] + x_left);
    fz[1]=1.0-((x_right_guess[1]-x_left)*m_cell[0].area)*rho/m_cell[0].mass; 
    riemann_interpolate(x_mid_guess[1],x_diaphragm,delta_t,base_sol,&rho,&e,&u);
    x_right_guess[2]=(m_cell[0].mass + rho*m_cell[0].area*x_left)/(rho*m_cell[0].area);
    x_mid_guess[2] = 0.5*(x_right_guess[2] + x_left);
    fz[2]=1.0-((x_right_guess[2]-x_left)*m_cell[0].area)*rho/m_cell[0].mass; 

    while (fabs(fz[2])>1.0e-6)  {
        x_right_guess[0]=x_right_guess[1]; x_mid_guess[0]=x_mid_guess[1]; fz[0]=fz[1];
        x_right_guess[1]=x_right_guess[2]; x_mid_guess[1]=x_mid_guess[2]; fz[1]=fz[2];
        x_right_guess[2]=x_right_guess[1]-fz[1]*(x_right_guess[1]-x_right_guess[0])/(fz[1]-fz[0]);
        x_mid_guess[2] = 0.5*(x_right_guess[2] + x_left);
        riemann_interpolate(x_mid_guess[2],x_diaphragm,delta_t,base_sol,
                            &rho,&e,&u);
        fz[2]=1.0-((x_right_guess[2]-x_left)*m_cell[0].area)*rho/m_cell[0].mass; 
    }

    QD[0].gas->copy_values_from(*(m_cell[0].gas));	/* set mass fractions */
    QD[0].gas->rho=m_cell[0].mass / (m_cell[0].area * ( x_right_guess[2] - x_left));
    QD[0].gas->e[0]=e;
    gmodel->eval_thermo_state_rhoe(*(QD[0].gas));
    gmodel->eval_transport_coefficients(*(QD[0].gas));

    cell_sol[0].x_p=x_right_guess[2];
    cell_sol[0].rho_p=QD[0].gas->rho;
    cell_sol[0].e_p=QD[0].gas->e[0];
    cell_sol[0].u_p=u;
    cell_sol[0].P_p=QD[0].gas->p;
    cell_sol[0].T_p=QD[0].gas->T[0];

    delete QD[0].gas;
    return SUCCESS;
}

/// \brief Find the location of a cells left face according to RSP profile
/// \param x_right : right face axial location
/// \param x_diaphragm : pointer to all slug data
/// \param delta_t : time lapse since simulated rupture
/// \param m_cell : the cell under investigation
/// \param base_sol : the rsp base solution structure
/// \param cell_sol : the interpolated properties for the new cell centroid
/// \returns : success or failure
int find_x_pos_left( double x_right, 
		     double x_diaphragm, 
		     double delta_t, 
		     struct L_cell * m_cell, 
		     struct ersp_solution * base_sol, 
		     struct ersp_solution * cell_sol)
{
    double x_left_guess[3], x_mid_guess[3], fz[3];
    double rho,e,u;
    struct L_flow_state QD[1];
    Gas_model *gmodel = get_gas_model_ptr();
    QD[0].gas = new Gas_data(gmodel);

    riemann_interpolate(x_right,x_diaphragm,delta_t,base_sol,&rho,&e,&u);
    x_left_guess[1]=(rho*m_cell[0].area*x_right-m_cell[0].mass)/(rho*m_cell[0].area);
    x_mid_guess[1] = 0.5*(x_left_guess[1] + x_right);
    fz[1]=1.0-((x_right-x_left_guess[1])*m_cell[0].area)*rho/m_cell[0].mass; 
    riemann_interpolate(x_mid_guess[1],x_diaphragm,delta_t,base_sol,&rho,&e,&u);
    x_left_guess[2]=(rho*m_cell[0].area*x_right-m_cell[0].mass)/(rho*m_cell[0].area);
    x_mid_guess[2] = 0.5*(x_right + x_left_guess[2]);
    fz[2]=1.0-((x_right-x_left_guess[2])*m_cell[0].area)*rho/m_cell[0].mass; 

    while (fabs(fz[2])>1.0e-6)  {
        x_left_guess[0]=x_left_guess[1]; x_mid_guess[0]=x_mid_guess[1]; fz[0]=fz[1];
        x_left_guess[1]=x_left_guess[2]; x_mid_guess[1]=x_mid_guess[2]; fz[1]=fz[2];
        x_left_guess[2]=x_left_guess[1]-fz[1]*(x_left_guess[1]-x_left_guess[0])/(fz[1]-fz[0]);
        x_mid_guess[2] = 0.5*(x_left_guess[2] + x_right);
        riemann_interpolate(x_mid_guess[2],x_diaphragm,delta_t,base_sol,
                            &rho,&e,&u);
        fz[2]=1.0-((x_right-x_left_guess[2])*m_cell[0].area)*rho/m_cell[0].mass; 
    }

    QD[0].gas->copy_values_from(*(m_cell[0].gas));
    QD[0].gas->rho=m_cell[0].mass / (m_cell[0].area * ( x_right - x_left_guess[2]));
    QD[0].gas->e[0]=e;
    gmodel->eval_thermo_state_rhoe(*(QD[0].gas));
    gmodel->eval_transport_coefficients(*(QD[0].gas));

    cell_sol[0].x_p=x_left_guess[2];
    cell_sol[0].rho_p=QD[0].gas->rho;
    cell_sol[0].e_p=QD[0].gas->e[0];
    cell_sol[0].u_p=u;
    cell_sol[0].P_p=QD[0].gas->p;
    cell_sol[0].T_p=QD[0].gas->T[0];

    delete QD[0].gas;
    return SUCCESS;
}

/// \brief Initialise diaphragm sub-problem
/// \param RSD : pointer to the riemann slug data structure
/// \param A : pointer to the gas slug data structure
/// \param D : pointer to diaphragm data structure
/// \returns : success or failure
int L_riemann_initialise( riemann_simulation_data * RSD, 
                          std::vector<slug_data> &A,
                          diaphragm_data * D )
{
    double right_mass, left_mass, right_err, left_err;
    double left_end_dx_old, right_end_dx_old;
    int i,j, ncells_left[6], ncells_right[6];

    /* Prep the RSD structure */
    RSD->patch=true; RSD->solve=true;
    RSD->RSP_dt = D->RSP_dt;
    RSD->rsr = D->right_slug_id;
    RSD->rsr_end_id = D->right_slug_end_id;
    RSD->right_end_dx = D->right_slug_dx;
    RSD->rsl = D->left_slug_id;
    RSD->rsl_end_id = D->left_slug_end_id;
    RSD->left_end_dx = D->left_slug_dx;

    /* Assign diaphragm position to this structure */
    RSD->x_diaphragm = A[RSD->rsl].Cell[A[RSD->rsl].ixmax].x;
    /* Pull out LHS & RHS states for Riemann solver */
    L_slug_end_properties( &(A[RSD->rsr]), RSD->rsr_end_id, RSD->right_end_dx,
                           &right_mass, RSD->QR );
    L_slug_end_properties( &(A[RSD->rsl]), RSD->rsl_end_id, RSD->left_end_dx, 
                           &left_mass, RSD->QL );
    /* Obtain an initial guess */
    exact_riemann(RSD->QL,RSD->QR,&(RSD->left_vel),&(RSD->cs_vel),
                  &(RSD->right_vel),RSD->base_sol,2);


    /* Need to iterate here to ensure the Riemann solution is generated by the 
     * average of the properties embodied in the cells that will be included
     * in the solution.  This is very important as without it there will be 
     * large discontinuities at either extremity of the RSP domain.
     */

    for (j=0; j<=5; j++)  {
        left_err = 1.0;
        i=0; left_mass=0.0;
        RSD->left_end_dx = 0.0;
        /* Find what cells will be included in the LHS slug */
        while (left_err>0.0)  {
            L_slug_end_properties( &(A[RSD->rsl]), RSD->rsl_end_id, 
                                   RSD->left_end_dx, &left_mass, RSD->QL );
            left_err=(-1)*RSD->left_vel*RSD->RSP_dt-RSD->left_end_dx;
            left_end_dx_old = RSD->left_end_dx;
            RSD->left_end_dx = A[RSD->rsl].Cell[A[RSD->rsl].ixmax].x-
                               A[RSD->rsl].Cell[A[RSD->rsl].ixmax-i].xmid;
            i++;
            }
        ncells_left[j]=i;
        right_err = 1.0;
        i=0; right_mass=0.0;
        RSD->right_end_dx = 0.0;
        /* Find what cells will be included in the RHS slug */
        while (right_err>0.0)  {
            L_slug_end_properties( &(A[RSD->rsr]), RSD->rsr_end_id, RSD->right_end_dx,
                                   &right_mass, RSD->QR );
            right_err=RSD->right_vel*RSD->RSP_dt-RSD->right_end_dx;
            right_end_dx_old = RSD->right_end_dx;
            RSD->right_end_dx = A[RSD->rsr].Cell[A[RSD->rsr].ixmin+i].xmid-
                                A[RSD->rsl].Cell[A[RSD->rsl].ixmax].x;
            i++;
        }
        ncells_right[j]=i;
        printf("   Iteration %d LHS cells = %d; RHS cells = %d \n",
                j,ncells_left[j],ncells_right[j]);

        if (j>0)   {
            if (ncells_left[j-1]==ncells_left[j] && ncells_right[j-1]==ncells_right[j])  {
                printf("\n#### converged solution after %d iteration/s ####\n\n",j);
                RSD->ncells_left=ncells_left[j];
                RSD->ncells_right=ncells_right[j];
                if (RSD->ncells_left > RIEMANN_CELLS || RSD->ncells_right > RIEMANN_CELLS)  {
                    printf("ERROR: RSP region is too large, reduce RSP_dt\n");
                    return FAILURE;
                }
                if (RSD->ncells_left > ERSP_NX || RSD->ncells_right > ERSP_NX)  {
                    printf("ERROR: Insufficient sampling resolution, increase ERSP_NX\n");
                    return FAILURE;
                }
                exact_riemann(RSD->QL,RSD->QR,&(RSD->left_vel),&(RSD->cs_vel),
                              &(RSD->right_vel),RSD->base_sol,1);
                break;
            } else if (j==5)   {
                printf("ERROR: solution not converging, excessive flow state variation\n");
                return FAILURE;
            } else  {
            /* keep iterating with PG wave speeds*/
            exact_riemann(RSD->QL,RSD->QR,&(RSD->left_vel),&(RSD->cs_vel),
                          &(RSD->right_vel),RSD->base_sol,2);
            }
        } else  {
            /* Don't check for convergence on first run */
            exact_riemann(RSD->QL,RSD->QR,&(RSD->left_vel),&(RSD->cs_vel),
            &(RSD->right_vel),RSD->base_sol,2);
        }
    }

#   if DEBUG >= 1
    printf("\nmass left and right are %e and %e\n", left_mass, right_mass);
    printf("Properties at diaphragm %d rupture:\n\
            left_p = %0.1f, left_T = %0.2f, left_u = %0.1f \n\
            right_p = %0.1f, right_T = %0.2f, right_u = %0.1f\n",
            RSD->jd, RSD->QL[0].gas->p, RSD->QL[0].gas->T[0], RSD->QL[0].u, 
            RSD->QR[0].gas->p, RSD->QR[0].gas->T[0], RSD->QR[0].u);
    printf("left vel is %0.1f, cs_vel is %0.1f and right vel is %0.1f\n",
            RSD->left_vel, RSD->cs_vel, RSD->right_vel);
#   endif

    return SUCCESS;
}

/// \brief Compute a solution for RSP (doesn't patch the new data in)
/// \param RSD : pointer to the riemann slug data structure
/// \param A : pointer to gas slug data
/// \returns : success or failure
int L_riemann_solve(riemann_simulation_data *RSD, std::vector<slug_data> &A)
{
    struct ersp_solution *lhs_cells, *rhs_cells, *cell_sol;
    int m, i;
    double shock_pos,cell_pos,cs_pos,x_left=0.0,x_right=0.0;
    double mass_included;
    double rwh_pos, rwh_mass, sa_dens, ra_dens;
    double nra_dens, global_shift=0.0, sa_length=0.0;
    double nls_dens, rare_adj_length, right_shift;
    double drho, delta_rho, de, du;
    int rc, sc, m_extra=0;
    struct L_flow_state QD[1];
    Gas_model *gmodel = get_gas_model_ptr();
    QD[0].gas = new Gas_data(gmodel);

    /* need to allocate memory to cells_sols here */
    lhs_cells = (struct ersp_solution *) calloc(RIEMANN_CELLS, 
                      sizeof(struct ersp_solution));
    rhs_cells = (struct ersp_solution *) calloc(RIEMANN_CELLS, 
                       sizeof(struct ersp_solution));
    cell_sol = (struct ersp_solution *) calloc(1, 
                sizeof(struct ersp_solution));
    /* Make a structure of these LHS cells
     * (note the order + 5 to allow for corrections) */
    for (m=0; m<RSD->ncells_left+5; m++)  {
        RSD->rs_lcells[m]=A[RSD->rsl].Cell[A[RSD->rsl].ixmax-m];
    }
    /* Make a structure of these RHS cells (note the order and 
     * + 5 to allow for corrections) */
    for (m=0; m<RSD->ncells_right+5; m++)  {
        RSD->rs_rcells[m]=A[RSD->rsr].Cell[A[RSD->rsr].ixmin+m];
    }
    /* Locate the cell the shock is currently in */
    shock_pos = RSD->RSP_dt*RSD->right_vel+RSD->x_diaphragm;
    mass_included=0.0;
    for (m=0; m<RSD->ncells_right; m++)  {
        /* Note that cell_pos is the right hand edge of the cell */
        cell_pos = RSD->rs_rcells[m].x;
        mass_included += RSD->rs_rcells[m].mass;
        if (cell_pos>shock_pos) break;
        if (m==RSD->ncells_right-1) {
            printf("ERROR: shock has moved outside predicted range\n");
            return FAILURE;
        }
    }
    sc=m-1;
    cs_pos = RSD->RSP_dt*RSD->cs_vel+RSD->x_diaphragm;
    if (sc<=2)  {
        /* this occurs when shock hasn't travelled past few cells yet, 
         * don't do anything */
        printf("ERROR: Patching window is too small [sc=%d], increase RSP_dt...\n",sc);
        return FAILURE;
    } else  {
        for (m=0;m<=sc;m++)  {
            if (m==0) {
                x_left=cs_pos;
            } else {
                x_left=x_right;
            }
            find_x_pos_right(x_left,RSD->x_diaphragm,RSD->RSP_dt,
                             &(RSD->rs_rcells[m]),RSD->base_sol,cell_sol);
            x_right=cell_sol[0].x_p;
            rhs_cells[m]=cell_sol[0];
            /* make x the cell mid-point */
            rhs_cells[m].x_p=0.5*(x_right+x_left);
            rhs_cells[m].dx_p=x_right-x_left;
            if (m==sc)   {
                /* check next cell to see if its density makes sense */
                sa_dens = RSD->rs_rcells[m+1].mass / 
	        (RSD->rs_rcells[m+1].area * (RSD->rs_rcells[m+1].x - x_right));
                sa_length = (RSD->rs_rcells[m+1].x - x_right);
            }
        }
    }
    /* Now do same for left hand side of contact surface */
    rwh_pos = RSD->RSP_dt*RSD->left_vel+RSD->x_diaphragm;
    /* assume constant area */
    rwh_mass = (RSD->x_diaphragm - rwh_pos)*RSD->rs_lcells[0].area*RSD->QL[0].gas->rho;
    mass_included=0.0;
    /* need to use mass for this one as property profile is not constant */
    for (m=0; m<RSD->ncells_left; m++)  {
        mass_included += RSD->rs_lcells[m].mass;
        if (mass_included>rwh_mass) break;
        if (m==RSD->ncells_left-1)  {
            printf("ERROR: LHS cell searching has failed...\n");
            return FAILURE;
        }
    }
    rc=m+1;
    if (rc<=5)  {
        /* this occurs when expansion wave hasn't travelled past first few cells yet */
        printf("ERROR: Patching window is too small [rc=%d], increase RSP_dt...\n",sc);
        return FAILURE;
    } else  {
        for (m=0; m<=rc; m++)  {
            if (m==0) {
                x_right=cs_pos;
            } else {
                x_right=x_left;
            }
            find_x_pos_left(x_right,RSD->x_diaphragm,RSD->RSP_dt,
                            &(RSD->rs_lcells[m]),RSD->base_sol,cell_sol);
            x_left=cell_sol[0].x_p;
            lhs_cells[m]=cell_sol[0];
            lhs_cells[m].x_p=0.5*(x_left+x_right);
            lhs_cells[m].dx_p=x_right-x_left;
            if (m==rc)   {
                /* check next cell to see if its density makes sense */
                ra_dens = RSD->rs_lcells[m+1].mass / (RSD->rs_lcells[m+1].area
				     * (x_left-(2.0*RSD->rs_lcells[m+1].xmid 
				     - RSD->rs_lcells[m+1].x)));
                rare_adj_length = x_left-(2.0*RSD->rs_lcells[m+1].xmid - RSD->rs_lcells[m+1].x);
                /* Now need to see if a global left or right shift is required to fix 
                 * rarefaction adjacent cell */
                global_shift = rare_adj_length-RSD->rs_lcells[m+1].mass / 
			       (0.5*(RSD->QL[0].gas->rho + lhs_cells[m].rho_p) 
			       * RSD->rs_lcells[m+1].area );
            }
        }
    }
    nra_dens =  RSD->rs_rcells[sc+1].mass / (RSD->rs_rcells[sc+1].area 
                          * (sa_length+global_shift));
    /* propagate global_shift to all cells involved */
    for (m=0; m<=rc; m++)   {
        lhs_cells[m].x_p-=global_shift;
        }
    for (m=0; m<=sc; m++)   {
        rhs_cells[m].x_p-=global_shift;
        }
    /* fill out rarefaction adjacent cell with new (forced) properties */
    lhs_cells[rc+1].dx_p = lhs_cells[rc].x_p-0.5*
		lhs_cells[rc].dx_p-
		(2.0*RSD->rs_lcells[rc+1].xmid-RSD->rs_lcells[rc+1].x);
    lhs_cells[rc+1].x_p=2.0*RSD->rs_lcells[rc+1].xmid-
		RSD->rs_lcells[rc+1].x + 0.5 * lhs_cells[rc+1].dx_p;
    lhs_cells[rc+1].rho_p=RSD->rs_lcells[rc+1].mass/
		(RSD->rs_lcells[rc+1].area*lhs_cells[rc+1].dx_p);
    lhs_cells[rc+1].u_p=0.5*(lhs_cells[rc].u_p+RSD->QL[0].u);
    lhs_cells[rc+1].e_p=0.5*(lhs_cells[rc].e_p+RSD->QL[0].gas->e[0]);
    QD[0].gas->copy_values_from(*(RSD->QL[0].gas));
    QD[0].gas->rho = lhs_cells[rc+1].rho_p;
    QD[0].gas->e[0] = lhs_cells[rc+1].e_p;
    gmodel->eval_thermo_state_rhoe(*(QD[0].gas));
    gmodel->eval_transport_coefficients(*(QD[0].gas));
    lhs_cells[rc+1].P_p=QD[0].gas->p;
    lhs_cells[rc+1].T_p=QD[0].gas->T[0];

    /* Alter shock adjacent cell accordingly */
    rhs_cells[sc+1].dx_p = RSD->rs_rcells[sc+1].x
		    - (rhs_cells[sc].x_p + rhs_cells[sc].dx_p);

    /* Now fit cells into existing slugs without creating bad adj. cells */
    if (nra_dens<RSD->QR[0].gas->rho)   {
        /* retract shock cell by cell until all densities are ok */
        for (m=0; m<5; m++)   {
            /* right_shift relates to LHS face of the cell adjacent to the shock */
            right_shift = rhs_cells[sc+1-m].dx_p-RSD->rs_rcells[sc+1-m].mass/
		          (RSD->QR[0].gas->rho*RSD->rs_rcells[sc+1-m].area);

            rhs_cells[sc-m].dx_p+=right_shift;
            rhs_cells[sc-m].x_p+=0.5*right_shift;
            rhs_cells[sc+1-m].dx_p-=right_shift;
            rhs_cells[sc+1-m].x_p=rhs_cells[sc-m].x_p+
			    0.5*rhs_cells[sc-m].dx_p+
			    0.5*rhs_cells[sc+1-m].dx_p;
            nls_dens = RSD->rs_rcells[sc-m].mass/
			    (RSD->rs_rcells[sc-m].area *
			    (rhs_cells[sc-m].dx_p));
            rhs_cells[sc+1-m].u_p = RSD->QR[0].gas->p;
            rhs_cells[sc+1-m].e_p = RSD->QR[0].gas->e[0];
            rhs_cells[sc+1-m].T_p = RSD->QR[0].gas->T[0];
            rhs_cells[sc+1-m].P_p = RSD->QR[0].gas->p;
            rhs_cells[sc+1-m].rho_p = RSD->QR[0].gas->rho;

            if (nls_dens<rhs_cells[sc-m].rho_p && 
				nls_dens>RSD->QR[0].gas->rho)  {
                # if (DEBUG>=1)
                printf("RHS: ad_cell density = %e by retracting %d cell\n",
                       nls_dens,m+1);
                # endif
                /* Assign properties to sc - m cell */
                rhs_cells[sc-m].rho_p=nls_dens;
                drho = (RSD->QR[0].gas->rho-rhs_cells[sc-1-m].rho_p);
                delta_rho = (rhs_cells[sc-m].rho_p - 
			     rhs_cells[sc-1-m].rho_p);
                de = (RSD->QR[0].gas->e[0]-rhs_cells[sc-1-m].e_p);
                du = (RSD->QR[0].u-rhs_cells[sc-1-m].u_p);
                QD[0].gas->copy_values_from(*(RSD->QR[0].gas));
                QD[0].gas->rho = rhs_cells[sc-m].rho_p;
                QD[0].gas->e[0] = de / drho * delta_rho +  rhs_cells[sc-1-m].e_p;
                gmodel->eval_thermo_state_rhoe(*(QD[0].gas));
		gmodel->eval_transport_coefficients(*(QD[0].gas));
                rhs_cells[sc-m].u_p = du / drho * delta_rho + 
			    rhs_cells[sc-1-m].u_p;
                rhs_cells[sc-m].e_p = QD[0].gas->e[0];
                rhs_cells[sc-m].T_p = QD[0].gas->T[0];
                rhs_cells[sc-m].P_p = QD[0].gas->p;
                break;
            }
            else if (nls_dens>rhs_cells[sc-m].rho_p)   {
                printf("ERROR: Shock retraction failed\n");
                return FAILURE;
            }
            else if (nls_dens<RSD->QR[0].gas->rho)  {
	        /* Retract by another cell... 
                 * Properties already assigned to +1 cell that is now fixed */
                if (m==4) {
                    printf("ERROR: Shock retraction failed\n");
                    return FAILURE;
                }
            }
        }
    }
    else if (nra_dens>rhs_cells[sc].rho_p)  {
        /* extend shock by one cell */
        # if (DEBUG>=1)
        printf("Shock side requires shifting as adjacent cell has a too high a density\n");
        # endif
        /* extend shock cell by cell until all densities are ok */
        for (m=0; m<3; m++)   {
            /* right_shift relates to RHS face of the cell adjacent to the shock */
            right_shift = RSD->rs_rcells[sc+1+m].mass/
		    (rhs_cells[sc+m].rho_p*RSD->rs_rcells[sc+1+m].area)
		     - rhs_cells[sc+1+m].dx_p;
            rhs_cells[sc+1+m].dx_p+=right_shift;
            rhs_cells[sc+1+m].x_p+=0.5*right_shift;
            rhs_cells[sc+2+m].dx_p=RSD->rs_rcells[sc+2+m].x - 
			    RSD->rs_rcells[sc+1+m].x -right_shift;
            rhs_cells[sc+2+m].x_p=rhs_cells[sc+1+m].x_p+
			    0.5*rhs_cells[sc+1+m].dx_p+
			    0.5*rhs_cells[sc+2+m].dx_p;

            nls_dens = RSD->rs_rcells[sc+2+m].mass/
			    (RSD->rs_rcells[sc+2+m].area * 
			    (rhs_cells[sc+2+m].dx_p));

            rhs_cells[sc+1+m].u_p = rhs_cells[sc+m].u_p;
            rhs_cells[sc+1+m].e_p = rhs_cells[sc+m].e_p;
            rhs_cells[sc+1+m].T_p = rhs_cells[sc+m].T_p;
            rhs_cells[sc+1+m].P_p = rhs_cells[sc+m].P_p;
            rhs_cells[sc+1+m].rho_p = rhs_cells[sc+m].rho_p;

            if (nls_dens<rhs_cells[sc+1+m].rho_p 
		    && nls_dens>RSD->QR[0].gas->rho)   {
                m_extra=m+1;
                # if (DEBUG>=1)
                printf("LHS: ad_cell density = %e by retracting %d cell\n", 
                       nls_dens, m_extra);
                # endif
                /* Assign properties to sc + 2 + m cell */
                rhs_cells[sc+2+m].rho_p=nls_dens;
                drho = (RSD->QR[0].gas->rho-rhs_cells[sc+1+m].rho_p);
                delta_rho = (rhs_cells[sc+2+m].rho_p - 
			    rhs_cells[sc+1+m].rho_p);
                de = (RSD->QR[0].gas->e[0]-rhs_cells[sc+1+m].e_p);
                du = (RSD->QR[0].u-rhs_cells[sc+1+m].u_p);
                QD[0].gas->copy_values_from(*(RSD->QR[0].gas));
                QD[0].gas->rho = rhs_cells[sc+2+m].rho_p;
                QD[0].gas->e[0] = de / drho * delta_rho +  rhs_cells[sc+1+m].e_p;
                gmodel->eval_thermo_state_rhoe(*(QD[0].gas));
		gmodel->eval_transport_coefficients(*(QD[0].gas));
                rhs_cells[sc+2+m].u_p = du / drho * delta_rho + 
					    rhs_cells[sc+1+m].u_p;
                rhs_cells[sc+2+m].e_p = QD[0].gas->e[0];
                rhs_cells[sc+2+m].T_p = QD[0].gas->T[0];
                rhs_cells[sc+2+m].P_p = QD[0].gas->p;
                break;
            }
            else if (nls_dens>rhs_cells[sc+1+m].rho_p)   {
                /* Need to extend shock by another cell */
                /* Properties already assigned to +1 cell that is now fixed */
                if (m==2)  {
                    printf("ERROR: Shock extension failed\n");
                    return FAILURE;
                }
            }
            else if (nls_dens<RSD->QR[0].gas->rho)   {
            printf("ERROR: adcell density is too high\n");
            return FAILURE;
            }
        }
    }
    else   {
        /* The new adjacent cell density is okay!
         * Give the adjacent cell some flow properties via 
         * linear interpolation between bounding cells */

        rhs_cells[sc+1].rho_p = nra_dens;
        rhs_cells[sc+1].dx_p = sa_length+global_shift;
        rhs_cells[sc+1].x_p = RSD->rs_rcells[sc+1].x -
	    	 0.5*rhs_cells[sc+1].dx_p;

        drho = (RSD->QR[0].gas->rho-rhs_cells[sc].rho_p);
        delta_rho = rhs_cells[sc+1].rho_p - rhs_cells[sc].rho_p;
        de = (RSD->QR[0].gas->e[0]-rhs_cells[sc].e_p);
        du = (RSD->QR[0].u-rhs_cells[sc].u_p);
        QD[0].gas->copy_values_from(*(RSD->QR[0].gas));
        QD[0].gas->rho = rhs_cells[sc+1].rho_p;
        QD[0].gas->e[0] = de / drho * delta_rho +  rhs_cells[sc].e_p;
        gmodel->eval_thermo_state_rhoe(*(QD[0].gas));
	gmodel->eval_transport_coefficients(*(QD[0].gas));
        rhs_cells[sc+1].u_p = du / drho * delta_rho + rhs_cells[sc].u_p;
        rhs_cells[sc+1].e_p = QD[0].gas->e[0];
        rhs_cells[sc+1].T_p = QD[0].gas->T[0];
        rhs_cells[sc+1].P_p = QD[0].gas->p;
    }

    /*  Append the quantities of mass and volume to the solution structure */
    for (m=0;m<=sc+1+m_extra;m++)   {
        rhs_cells[m].mass = RSD->rs_rcells[m].mass;
        rhs_cells[m].vol = RSD->rs_rcells[m].area * rhs_cells[m].dx_p;
    }

    for (m=0;m<=rc+1;m++)   {
        lhs_cells[m].mass = RSD->rs_lcells[m].mass;
        lhs_cells[m].vol = RSD->rs_lcells[m].area * lhs_cells[m].dx_p;
    }

    /* Now wite over appropriate values in RSD structure */
    for (i=0; i<=rc+1; i++)   {
        /* Over-write old cell values */
        RSD->rs_lcells[i].xmid = lhs_cells[i].x_p;
        RSD->rs_lcells[i].x = lhs_cells[i].x_p+0.5*lhs_cells[i].dx_p;
        RSD->rs_lcells[i].gas->p = lhs_cells[i].P_p;
        RSD->rs_lcells[i].u = lhs_cells[i].u_p;
        RSD->rs_lcells[i].gas->T[0] = lhs_cells[i].T_p;
        RSD->rs_lcells[i].gas->rho = lhs_cells[i].rho_p;
        RSD->rs_lcells[i].gas->e[0] = lhs_cells[i].e_p;
        RSD->rs_lcells[i].volume = lhs_cells[i].vol;
    }
    for (i=0; i<=sc+1+m_extra; i++)  {
        /* Over-write old cell values */
        RSD->rs_rcells[i].xmid = rhs_cells[i].x_p;
        RSD->rs_rcells[i].x = rhs_cells[i].x_p+0.5*rhs_cells[i].dx_p;
        RSD->rs_rcells[i].gas->p = rhs_cells[i].P_p;
        RSD->rs_rcells[i].u = rhs_cells[i].u_p;
        RSD->rs_rcells[i].gas->T[0] = rhs_cells[i].T_p;
        RSD->rs_rcells[i].gas->rho = rhs_cells[i].rho_p;
        RSD->rs_rcells[i].gas->e[0] = rhs_cells[i].e_p;
        RSD->rs_rcells[i].volume = rhs_cells[i].vol;	
    }

    /* Append shock and rarefaction wave cells to RSD structure */

    RSD->sc = sc;
    RSD->rc = rc;
    RSD->m_extra = m_extra;

    /* Free allocated memory */

    free(lhs_cells);
    free(rhs_cells);
    free(cell_sol);

    delete QD[0].gas;
    return SUCCESS;
}	/* end L_riemann_solve */

/// \brief Patch in the computed RSP solution
/// \param RSD : pointer to the riemann slug data structur
/// \param A[] : pointer to gas slug structures
/// \returns : success or failure
int L_riemann_patch(riemann_simulation_data *RSD,
		    std::vector<slug_data> &A)
{
    int i;
    /* Patch in Riemann solution */
    for (i=0; i<=RSD->rc+1; i++)   {
        L_copy_cell_data(&(RSD->rs_lcells[i]),
                         &(A[RSD->rsl].Cell[A[RSD->rsl].ixmax-i]),1);
    }
    for (i=0; i<=RSD->sc+1+RSD->m_extra; i++)   {
        L_copy_cell_data(&(RSD->rs_rcells[i]),
                         &(A[RSD->rsr].Cell[A[RSD->rsr].ixmin+i]),1);
    }

    L_encode_conserved(&(A[RSD->rsr]));	/* Mass, momentum, energy */
    L_encode_conserved(&(A[RSD->rsl]));
    L_decode_conserved(&(A[RSD->rsr]));	/* Entropy, gas properties */
    L_decode_conserved(&(A[RSD->rsl]));

    L_exchange_bc_data(&(A[RSD->rsl]), &(A[RSD->rsr]));
    /* Fix the 2nd ghost cell */
    A[RSD->rsr].Cell[1].x = A[RSD->rsl].Cell[A[RSD->rsl].ixmax].x;

    return SUCCESS;
}
