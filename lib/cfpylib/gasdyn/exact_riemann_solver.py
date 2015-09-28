from gaspy import *
from math import *

tol=1.0e-10

# The managed gas model lives here.
gmodel = None

def set_gas_model_ptr(gmptr):
    global gmodel
    gmodel = gmptr

class pg_int_state:
    def __init__(self):
        self.star = 0.0
        self.pstar_d = 0.0
        self.astar = 0.0
        self.v_rwt = 0.0
        self.v_cs = 0.0
        self.v_sw = 0.0

class L_flow_state:
    def __init__(self,_gas=None,_u=0.0):
        if _gas==None:
            if gmodel==None:
                raise RuntimeError("cannot create L_flow_state without a gas model")
            else:
                self.gas = Gas_data(gmodel)
        else:
            self.gas = Gas_data(_gas)
        self.u = _u


def pow(a,b):
    return a**b

def sign_check( a, b ):
    c=a/b
    if (c<0): return 0;
    else: return 1;

def shock_wave( QI, ustar, left ):
    uI = QI.u;
    gI = gmodel.gamma(QI.gas);
    pI = QI.gas.p;
    aI = QI.gas.a;
    cI = gI*pI/aI;

    if (left==0): wI = (gI+1)/4*(ustar-uI)/aI+sqrt(1+pow((gI+1)/4*(ustar-uI)/aI,2));
    else: wI = (gI+1)/4*(ustar-uI)/aI-sqrt(1+pow((gI+1)/4*(ustar-uI)/aI,2));

    pIstar = pI+cI*(ustar-uI)*wI;
    pIstar_d = 2*cI*pow(wI,3)/(1+pow(wI,2));
    aIstar = aI*sqrt(((gI+1)+(gI-1)*pIstar/pI)/((gI+1)+(gI-1)*pI/pIstar));

    QIstar_pg = pg_int_state()

    QIstar_pg.pstar = pIstar;
    QIstar_pg.pstar_d = pIstar_d;
    QIstar_pg.astar = aIstar;

    return wI, QIstar_pg

def rarefaction(  QI, ustar, left ):
    sign=1.0;
    
    uI = QI.u;
    gI = gmodel.gamma(QI.gas);
    pI = QI.gas.p;
    aI = QI.gas.a;

    if (left==0):  sign=1.0;
    elif (left==1):  sign=-1.0;
    aIstar = aI+sign*(gI-1)/2*(ustar-uI);
    pIstar = pI*pow(aIstar/aI,2*gI/(gI-1));
    pIstar_d = sign*gI*pIstar/aIstar;

    QIstar_pg = pg_int_state()

    QIstar_pg.pstar = pIstar;
    QIstar_pg.pstar_d = pIstar_d;
    QIstar_pg.astar = aIstar;

    return QIstar_pg

def perfect_gas ( QL, QR ):

    uL = QL.u; pL = QL.gas.p
    aL = QL.gas.a; gL = gmodel.gamma(QL.gas)

    uR = QR.u; pR = QR.gas.p;
    aR = QR.gas.a; gR = gmodel.gamma(QR.gas)

    # Determine wave pattern a-priori

    u_rcvr = uL+2*aL/(gL-1)+2*aR/(gR-1)

    if (pR > pL):
        sigma = gR;
        u_scn = uL-(aL/gL)*(pR/pL-1)/sqrt((gL+1)/(2*gL)*pR/pL+(gL-1)/(2*gL))
        u_ncr = uL+2*aR/(gR-1)*(1-pow(pL/pR,(gR-1)/(2*gR)))
        if (uR < u_scn):  designation = 0             # #  SCS */
        elif (uR == u_scn):  designation = 99      # #  SCN */
        elif (uR < u_ncr):  designation = 1        # #  SCR */
        elif (uR == u_ncr):  designation = 98      # #  NCR */
        elif (uR < u_rcvr):  designation = 2       # #  RCR */
        elif (uR == u_rcvr): designation = 2       # #  RCVR (NA) */
        else: designation = 2;                         # #  RCVCR (NA) */
    else:  ##  pR <= pL */
        sigma = gL
        u_ncs = uL-(aR/gR)*(pL/pR-1)/sqrt((gR+1)/(2*gR)*pL/pR+(gR-1)/(2*gR))
        u_rcn = uL+2*aL/(gL-1)*(1-pow(pR/pL,(gL-1)/(2*gL)))
        if (uR < u_ncs):  designation = 0             # #  SCS */
        elif (uR == u_ncs):  designation = 96      # #  NCS */
        elif (uR < u_rcn):   designation = 3       # #  RCS */
        elif (uR == u_rcn):  designation = 95      # #  RCN */
        elif (uR < u_rcvr):  designation = 2       # #  RCR */
        elif (uR == u_rcvr): designation = 2        # #  RCVR (NA) */
        else: designation = 2                          # #  RCVCR (NA) */ 

    # Initial guesses from Gottlieb and Groth (1994)

    uLbar = uL+2.0/(gL-1.0)*aL
    uRbar = uR-2.0/(gR-1.0)*aR

    z=(gL-1)/(gR-1)*aR/aL*pow(pL/pR,(sigma-1)/(2*sigma))
    u0star = (uLbar*z+uRbar)/(1+z)

    # Unique iterative loop required for each wave pattern

    count=1
    epsilon=1
    ustar=u0star
    while (epsilon > tol):
        if ((designation==0) or (designation==1) or (designation==99)):
            WSL, QLstar_pg = shock_wave( QL, ustar, 1 )
            vL=uL+aL*WSL
        if ((designation==96) or (designation==3) or (designation==0)):
            WSR, QRstar_pg = shock_wave( QR, ustar, 0 )
            vR = uR+aR*WSR
        if ((designation==2) or (designation==3) or (designation==95)):
            QLstar_pg = rarefaction( QL, ustar, 1 )
        if ((designation==1) or (designation==2) or (designation==98)):
            QRstar_pg = rarefaction( QR, ustar, 0 )
        if ((designation==96) or (designation==98)) :
            QLstar_pg.pstar=QL.gas.p
            QLstar_pg.pstar_d=0
        if ((designation==95) or (designation==99)):
            QRstar_pg.pstar=QR.gas.p
            QRstar_pg.pstar_d=0

        ustar -= (QLstar_pg.pstar-QRstar_pg.pstar)/(QLstar_pg.pstar_d-QRstar_pg.pstar_d)
        epsilon=abs(1.0-QLstar_pg.pstar/QRstar_pg.pstar)
        count+=1

    pstar=(QLstar_pg.pstar+QRstar_pg.pstar)/2.0

    gL=gmodel.gamma(QL.gas); RL=gmodel.R(QL.gas)
    gR=gmodel.gamma(QR.gas); RR=gmodel.R(QR.gas)

    T_starL=pow(QLstar_pg.astar,2)/(gL*RL);
    rho_starL=(QLstar_pg.pstar)/(RL*T_starL);
    T_starR=pow(QRstar_pg.astar,2)/(gR*RR);
    rho_starR=(QRstar_pg.pstar)/(RR*T_starR);

    S_pg = pg_int_state()

    S_pg.pstar=pstar;
    S_pg.v_rwt = QL.u - QL.gas.a;
    S_pg.v_cs = ustar;
    S_pg.v_sw = vR;

    return S_pg

#  Now an exact Riemann solver for an arbitrary EOS */

def shock_energy ( Us, v_i, P_i, rho_i, P_star, Qi ):
    Q = L_flow_state()
    Q.gas = Gas_data(gmodel);

    Q.gas.copy_values_from(Qi.gas);           #  makes mass fractions correct */
    V_i=v_i-Us;
    V_star=V_i/abs(V_i)*abs(((P_i-P_star)+rho_i*pow(V_i,2))/(rho_i*V_i));
    rho_star=rho_i*V_i/V_star;
    e=Qi.gas.e[0]; rho=Qi.gas.rho; p=Qi.gas.p;

    Q.gas.rho = rho_star;
    Q.gas.p = P_star;
    gmodel.eval_thermo_state_rhop(Q.gas);
    gmodel.eval_transport_coefficients(Q.gas);
    e_star=Q.gas.e[0];

    delta_h = - (e + p/rho) + ( e_star + P_star/rho_star ) ;
    zero_val=(0.5*(pow(V_star,2)-pow(V_i,2))+delta_h)/(delta_h);

    return zero_val

def shock_jump ( P_star,QI,shock_right ):
    Q1 = L_flow_state()

    #  usually will be a right facing shock */
    if (shock_right==0): sign=-1;
    else: sign=1;

    rhoI=QI.gas.rho; TI=QI.gas.T[0]; uI=QI.u;
    PI=QI.gas.p; gI=gmodel.gamma(QI.gas);
    C_vI=gmodel.Cv(QI.gas);
    Q1.gas=Gas_data(QI.gas);
    dT_drhoI = gmodel.dTdrho_const_p_py(Q1.gas)
    Q1.gas=Gas_data(QI.gas)
    dT_dpI = gmodel.dTdp_const_rho_py(Q1.gas)

    Us_PG=sign*sqrt((PI/(2.0*rhoI))*(P_star/PI*(gI+1)+gI-1))+uI;
    Us_APRX=Us_PG*1.02;
    count=0;
    errors=shock_energy(Us_PG,uI,PI,rhoI,P_star,QI);
    abs_err=abs(errors);
    first_err=abs_err;

    Us_n=Us_PG;
    Us_nm1=Us_APRX;
    Us_np1=Us_PG;

    shortcut = False

    while (abs_err>tol):
        #  a check for divergence */
        if (count==4 and abs_err>first_err):
            ig_shift_flag+=1;
            if (ig_shift_flag==500)  :
                printf("WARNING: bad initial guess in shock_jump\n");
                Us_np1 = Us_PG;
                shortcut=True
            else:
                count=0;
                Us_nm1=Us_PG*pow(1.001,pow(-1.0,ig_shift_flag)*ig_shift_flag);
                Us_n=Us_PG;
        elif (count!=0):
            Us_nm1=Us_n;
            Us_n=Us_np1;
        if not shortcut:
            se_n=shock_energy(Us_n,uI,PI,rhoI,P_star,QI);
            se_nm1=shock_energy(Us_nm1,uI,PI,rhoI,P_star,QI);
            Us_np1=Us_n-se_n*(Us_n-Us_nm1)/(se_n-se_nm1);
            errors=shock_energy(Us_np1,uI,PI,rhoI,P_star,QI);
            abs_err=abs(errors);
            count+=1;

    Us_I=Us_np1;
    VI=uI-Us_np1;
    V_starI=VI/abs(VI)*abs(((PI-P_star)+rhoI*pow(VI,2))/(rhoI*VI)); 
    u_starI=V_starI+Us_I;
    rho_starI=rhoI*VI/V_starI;

    QIstar = L_flow_state()

    QIstar.gas=Gas_data(QI.gas);
    QIstar.u=u_starI;
    QIstar.gas.rho=rho_starI;
    QIstar.gas.p=P_star;
    gmodel.eval_thermo_state_rhop(QIstar.gas);
    gmodel.eval_transport_coefficients(QIstar.gas);

    return Us_I,QIstar

def rare_jump ( P_star, QI, rare_left):
    u_starI=QI.u; rho_scale=1.0e-3; rho_f=0.0;
    T=0.0; rho=0.0
    run=1; dens_too_high=1
    rho_limit=0.1; drho_sign=-1.0;
    Q1 = L_flow_state()
    Q = L_flow_state()

    if (rare_left==0): sign=1;
    else: sign=-1;


    Q1.gas.copy_values_from(QI.gas);
    Q.gas.copy_values_from(QI.gas);               #  Ensures correct mass fraction */
    #  Step along the isentrope until P_star is reached */
    Q.gas.T[0]=QI.gas.T[0]; Q.gas.rho=QI.gas.rho; C_v=gmodel.Cv(QI.gas);
    rho=QI.gas.rho; T=QI.gas.T[0]; u_starI = QI.u;
    gmodel.eval_thermo_state_rhoT(Q.gas);
    gmodel.eval_transport_coefficients(Q.gas);
    while (run==1):
        d_rho=drho_sign*rho*rho_scale;
        #  Isentrope diffs */
        C_v=gmodel.Cv(Q.gas);
        Q1.gas=Gas_data(QI.gas);
        dT_dp = gmodel.dTdp_const_rho_py(Q1.gas);
        dT_drho_s = T/(C_v*pow(rho,2)*dT_dp);
        #  Calculate step for isobar only if a reasonable density exists */
        if (dens_too_high==0):
            Q.gas.rho=rho; Q.gas.p=P_star;
            gmodel.eval_thermo_state_rhop(Q.gas);
            T_ib=Q.gas.T[0];
            Q1.gas.copy_values_from(Q.gas);
            dT_drho_ib = gmodel.dTdrho_const_p_py(Q1.gas);
            Q1.gas.copy_values_from(Q.gas);
            dT_dp_ib = gmodel.dTdp_const_rho_py(Q1.gas);
            f=T-T_ib;
            df_drho=dT_drho_s-dT_drho_ib;               #  combine diffs */
            T+=dT_drho_s*d_rho;                         #  update both temperatures */ 
            T_ib+=dT_drho_ib*d_rho;
            rho+=d_rho;                                 #  step density forward */
            f_new=T-T_ib;
            if (sign_check(f,f_new)==0):
                rho_f=(rho-d_rho)-f/df_drho;            #  linear interpolation */
                T-=dT_drho_s*(rho_f-(rho-d_rho));
                d_rho = rho_f - (rho-d_rho);	        #  Fix d_rho for integration */ 
                run=0;                                  #  terminate run */
            else:
                #  check for right direction if stepping is to continue */
                fn_abs=abs(f_new);
                f_abs=abs(f);
                if (fn_abs>f_abs): drho_sign*=-1.0;      #  reverse density step */
        elif (dens_too_high==1):
            T+=dT_drho_s*d_rho; rho+=d_rho;
            #  so that bad temperature warnings do not occur */
            rho_limit=P_star/(gmodel.R(QI.gas)*200.0);	
            if (rho<rho_limit):  dens_too_high=0
        #  Calculate expanded velocity */
        Q.gas.rho=rho; Q.gas.T[0]=T;
        gmodel.eval_thermo_state_rhoT(Q.gas);
        u_starI += sign*(Q.gas.a / Q.gas.rho) * d_rho

    QIstar = L_flow_state(Q.gas)

    QIstar.gas.rho=rho_f;
    QIstar.gas.T[0]=T;
    QIstar.u=u_starI;
    QIstar.gas.p=P_star;
    gmodel.eval_thermo_state_rhop(QIstar.gas);
    gmodel.eval_transport_coefficients(QIstar.gas);

    return QIstar

def f_zero ( P_star, QL, QR ):
    PL=QL.gas.p; uL=QL.u;
    PR=QR.gas.p; uR=QR.u;

# Determine wave pattern for current internal pressure guess

    if (P_star/PR>1):
        Us_R, QRstar = shock_jump(P_star, QR, 1);
        u_starR=QRstar.u;
    elif (P_star==PR): u_starR=uR;
    else:
        QRstar = rare_jump(P_star, QR, 0);
        u_starR=QRstar.u;

    if (P_star/PL<1):
        QLstar = rare_jump(P_star, QL, 1);
        u_starL=QLstar.u
    elif (P_star==PL): u_starL=uL;
    else:
        Us_L = QLstar = shock_jump(P_star, QL, 0)
        u_starL=QLstar.u

    zero_val=u_starL/(u_starR)-1.0

    return zero_val;

def final_states ( P_star, QL, QR ):
    PL=QL.gas.p; uL=QL.u;
    PR=QR.gas.p; uR=QR.u;

    # Determine wave pattern for current internal pressure guess
    Us_L = 0.0; Us_R = 0.0

    if (P_star/PR>1):
        Us_R,QRstar=shock_jump(P_star, QR, 1)
        u_starR=QRstar.u;
    elif (P_star==PR): u_starR=uR;
    else:
        QRstar=rare_jump(P_star, QR, 0)
        u_starR=QRstar.u

    if (P_star/PL<1):
        QLstar=rare_jump(P_star, QL, 1);
        u_starL=QLstar.u;
    elif (P_star==PL):
        u_starL=uL;
    else:
        Us_L,QLstar=shock_jump(P_star, QL, 0);
        u_starL=QLstar.u;

    return QLstar, QRstar, Us_L, Us_R

def exact_riemann( QL, QR ):
    P_star_n=0.0; P_star_np1=0.0;
    tol=1.0e-6
    count=0

    gmodel.eval_thermo_state_pT(QL.gas)
    gmodel.eval_thermo_state_pT(QR.gas)

    #*****************************************************************************
    # * Perfect gas solution
    #*****************************************************************************/

    S_pg = perfect_gas( QL, QR )

    #*****************************************************************************
    # * Arbitrary EOS solution
    #*****************************************************************************/

    P_star_0 = 0.95*S_pg.pstar
    P_star_1 = 0.98*S_pg.pstar
    P_star_np1=P_star_1
    errors=f_zero(P_star_1,QL,QR)
    abs_err=abs(errors)

    print "The initial guess for the real solver has an error of %g \n" % errors

    # Iterate for internal pressure using velocity as iterative variable */
    while (abs_err>tol):
        if (count==0):
            P_star_n=P_star_1; P_star_nm1=P_star_0;
        else:
            P_star_nm1=P_star_n; P_star_n=P_star_np1;
        P_star_np1=P_star_n-f_zero(P_star_n,QL,QR)*(P_star_n-P_star_nm1)/(f_zero(P_star_n,QL,QR)-f_zero(P_star_nm1,QL,QR))
        abs_err=abs(f_zero(P_star_np1,QL,QR))
        count+=1
        print "at iteration %d the abs_err is %g" % ( count,abs_err )
        if (count==10): break

    QLstar,QRstar,Us_L,Us_R = final_states(P_star_np1,QL,QR)
    print "************* Exact Riemann solution found ***************"
    print "* Iterations = %d, Relative error in velocity = %0.3e *" %(count,abs_err)
    print "* Us = %0.1f m/s; P* = %0.1f kPa; u* = %0.1f m/s"%(Us_R, 0.001*P_star_np1, QLstar.u)
    print "**********************************************************\n"

    return Us_R, QRstar, QLstar

