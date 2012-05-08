/** \file MBWR-EOS.cxx
 *  \ingroup gas
 *
 *  \author Peter Blyton
 *  \version 18-Apr-2012
 */

#include <iostream>
#include <sstream>
#include <cmath>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "physical_constants.hh"
#include "MBWR-EOS.hh"

using namespace std;

const double MBWR_PROPS_CONVERGE_TOL = 0.001;
const int MBWR_MAX_ITERATIONS = 20;

MBWR_EOS::
MBWR_EOS(lua_State *L)
    : Equation_of_state()
{
    lua_getglobal(L, "species"); // Put the list of species table to TOS

    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "MBWR_EOS::MBWR_EOS():\n";
        ost << "Error in the declaration of species: a table is expected.\n";
        input_error(ost);
    }

    if ( lua_objlen(L, -1) != 1 ) {
        ostringstream ost;
        ost << "MBWR_EOS::MBWR_EOS():\n";
        ost << "Error in the declaration of species: single species expected.\n";
        input_error(ost);
    }

    lua_rawgeti(L, -1, 1); // Put the specific species name to TOS
    const char* sp = luaL_checkstring(L, -1);
    lua_pop(L, 1); // pop specific species name off stack
    lua_pop(L, 1); // pop list of species table off stack

    lua_getglobal(L, sp); // Now bring the specific species table to TOS
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "MBWR_EOS::MBWR_EOS()\n";
        ost << "Error locating information table for species: " << sp << endl;
        input_error(ost);
    }

    double M;
    M = get_positive_value(L, -1, "M");
    M_.push_back(M);
    
    rho_c_ = get_positive_value(L, -1, "rho_c");

    lua_getfield(L, -1, "MBWR_EOS_coeffs"); // Bring EOS coeffs to TOS
    b_ = get_vector(L, -1, "b");
    lua_pop(L, 1); // pop EOS coeffs off stack
    lua_pop(L, 1); // pop specific species table off stack
}

MBWR_EOS::
~MBWR_EOS() {}

int
MBWR_EOS::
s_eval_pressure(Gas_data &Q)
{
    Q.p = MBWR_pressure(Q.rho, Q.T[0], rho_c_, M_[0], b_);
    return SUCCESS;
}

int
MBWR_EOS::
s_eval_temperature(Gas_data &Q)
{
    int status;
    Q.T[0] = MBWR_temperature(Q.rho, Q.p, rho_c_, M_[0], b_, status);
    return status;
}

int
MBWR_EOS::
s_eval_density(Gas_data &Q)
{
    int status;
    Q.rho = MBWR_density(Q.T[0], Q.p, rho_c_, M_[0], b_, status);
    return status;
}

double
MBWR_EOS::
s_gas_constant(const Gas_data &Q, int &status)
{   // Return the "equivalent" R for the ideal gas equation, needed by CFD code.
    status = SUCCESS;
    return Q.p / (Q.rho * Q.T[0]);
}

double 
MBWR_EOS::
s_prho_ratio(const Gas_data &Q, int isp)
{
    return MBWR_pressure(Q.rho, Q.T[0], rho_c_, M_[0], b_)/Q.rho;
}

double
MBWR_EOS::
s_dTdp_const_rho(const Gas_data &Q, int &status)
{
    return 1.0/MBWR_dpdT(Q.rho, Q.T[0], rho_c_, M_[0], b_);
}

double
MBWR_EOS::
s_dTdrho_const_p(const Gas_data &Q, int &status)
{   // Use the partial derivative cyclic relation.
    return -s_dTdp_const_rho(Q, status)*s_dpdrho_const_T(Q, status);
}

double
MBWR_EOS::
s_dpdrho_const_T(const Gas_data &Q, int &status)
{
    status = SUCCESS;
    return MBWR_dpdrho(Q.rho, Q.T[0], rho_c_, M_[0], b_);
}

double
MBWR_EOS::
s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
{
    status = SUCCESS;
    return MBWR_dpdrho(Q.rho, Q.T[0], rho_c_, M_[0], b_);
}

double
MBWR_EOS::
s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
{
    status = SUCCESS;
    return MBWR_dpdT(Q.rho, Q.T[0], rho_c_, M_[0], b_);
}

double
MBWR_EOS::
s_integral_const_T_energy(const Gas_data &Q)
{
    return MBWR_integral_const_T_energy(Q.rho, Q.T[0], rho_c_, M_[0], b_);
}

double
MBWR_EOS::
s_dintegral_const_T_energy_dT(const Gas_data &Q)
{
    return MBWR_dintegral_const_T_energy_dT(Q.rho, Q.T[0], rho_c_, M_[0], b_);
}

double
MBWR_EOS::
s_integral_const_T_entropy(const Gas_data &Q)
{
    return MBWR_integral_const_T_entropy(Q.rho, Q.T[0], rho_c_, M_[0], b_);
}

double MBWR_pressure(double rho, double T, double rho_c, double M, const std::vector<double> &b)
{
    // MBWR p-rho-T equation of state.
    double R = PC_R_u/100; // Universal gas constant in units L.bar/(mol.K)
    rho = rho*1e-3/M; // Density from kg/m^3 to mol/L
    rho_c = rho_c*1e-3/M; // Critical density from kg/m^3 to mol/L
    double Tinv = 1/T;
    double gamma = 1/pow(rho_c, 2);
    double calc_exp = exp(-gamma*pow(rho, 2));
    double p = rho*R*T
        + pow(rho, 2)*(b[1]*T + b[2]*sqrt(T) + b[3] + Tinv*(b[4] + Tinv*b[5]))
        + pow(rho, 3)*(b[6]*T + b[7] + Tinv*(b[8] + Tinv*b[9]))
        + pow(rho, 4)*(b[10]*T + b[11] + Tinv*b[12])
        + pow(rho, 5)*b[13]
        + pow(rho, 6)*Tinv*(b[14] + Tinv*b[15])
        + pow(rho, 7)*b[16]*Tinv
        + pow(rho, 8)*Tinv*(b[17] + Tinv*b[18])
        + pow(rho, 9)*pow(Tinv, 2)*b[19]
        + pow(rho, 3)*pow(Tinv, 2)*(b[20] + Tinv*b[21])*calc_exp
        + pow(rho, 5)*pow(Tinv, 2)*(b[22] + pow(Tinv, 2)*b[23])*calc_exp
        + pow(rho, 7)*pow(Tinv, 2)*(b[24] + Tinv*b[25])*calc_exp
        + pow(rho, 9)*pow(Tinv, 2)*(b[26] + pow(Tinv, 2)*b[27])*calc_exp
        + pow(rho, 11)*pow(Tinv, 2)*(b[28] + Tinv*b[29])*calc_exp
        + pow(rho, 13)*pow(Tinv, 2)*(b[30] + Tinv*(b[31] + Tinv*b[32]))*calc_exp;
    return p*1e5; // Convert bar to Pa
}

double MBWR_dpdT(double rho, double T, double rho_c, double M, const std::vector<double> &b)
{
    // Derivative of pressure p-rho-T equation w.r.t temperature.
    // This is needed for the integrals for internal energy and entropy.
    // Also used to iteratively solve the EOS for temperature using Newtons method.
    double R = PC_R_u/100; // Universal gas constant, L.bar/(mol.K) = J/(100 mol.K)
    rho = rho*1e-3/M; // Density from kg/m^3 to mol/L
    rho_c = rho_c*1e-3/M; // Critical density from kg/m^3 to mol/L
    double Tinv = 1/T;
    double gamma = 1/pow(rho_c, 2);
    double calc_exp = exp(gamma*pow(rho, 2));
    double dpdT = rho*R
        + pow(rho, 2)*(b[1] + 0.5*b[2]/sqrt(T) - pow(Tinv, 2)*(b[4] + 2*Tinv*b[5]))
        + pow(rho, 3)*(b[6] - pow(Tinv, 2)*(b[8] + 2*Tinv*b[9]))
        + pow(rho, 4)*(b[10] - pow(Tinv, 2)*b[12])
        - pow(rho, 6)*pow(Tinv, 2)*(b[14] + 2*Tinv*b[15])
        - pow(rho, 7)*b[16]*pow(Tinv, 2)
        - pow(rho, 8)*pow(Tinv, 2)*(b[17] + 2*Tinv*b[18])
        - pow(rho, 9)*2*pow(Tinv, 3)*b[19]
        - pow(rho, 3)*pow(Tinv, 3)*(2*b[20] + 3*Tinv*b[21])*calc_exp
        - pow(rho, 5)*pow(Tinv, 3)*(2*b[22] + 4*pow(Tinv, 2)*b[23])*calc_exp
        - pow(rho, 7)*pow(Tinv, 3)*(2*b[24] + 3*Tinv*b[25])*calc_exp
        - pow(rho, 9)*pow(Tinv, 3)*(2*b[26] + 4*pow(Tinv, 2)*b[27])*calc_exp
        - pow(rho, 11)*pow(Tinv, 3)*(2*b[28] + 3*Tinv*b[29])*calc_exp
        - pow(rho, 13)*pow(Tinv, 3)*(2*b[30] + Tinv*(3*b[31] + 4*Tinv*b[32]))*calc_exp;
    return dpdT*1e5; // Convert bar/K to Pa/K
}

double MBWR_dpdrho(double rho, double T, double rho_c, double M, const std::vector<double> &b)
{
    // Derivative of pressure p-rho-T equation w.r.t density.
    // Used to iteratively solve the EOS for density using Newtons method.
    double R = PC_R_u/100; // Universal gas constant in units L.bar/(mol.K)
    rho = rho*1e-3/M; // Density from kg/m^3 to mol/L
    rho_c = rho_c*1e-3/M; // Critical density from kg/m^3 to mol/L
    double Tinv = 1/T;
    double gamma = 1/pow(rho_c, 2);
    double calc_exp = exp(-gamma*pow(rho, 2));
    double dpdrho = R*T
        + 2*rho*(b[1]*T + b[2]*sqrt(T) + b[3] + Tinv*(b[4] + Tinv*b[5]))
        + 3*pow(rho, 2)*(b[6]*T + b[7] + Tinv*(b[8] + Tinv*b[9]))
        + 4*pow(rho, 3)*(b[10]*T + b[11] + Tinv*b[12])
        + 5*pow(rho, 4)*b[13]
        + 6*pow(rho, 5)*Tinv*(b[14] + Tinv*b[15])
        + 7*pow(rho, 6)*b[16]*Tinv
        + 8*pow(rho, 7)*Tinv*(b[17] + Tinv*b[18])
        + 9*pow(rho, 8)*pow(Tinv, 2)*b[19]
        + pow(rho, 2)*pow(Tinv, 2)*(b[20] + Tinv*b[21])*(3 - 2*gamma*pow(rho, 2))*calc_exp
        + pow(rho, 4)*pow(Tinv, 2)*(b[22] + pow(Tinv, 2)*b[23])*(5 - 2*gamma*pow(rho, 2))*calc_exp
        + pow(rho, 6)*pow(Tinv, 2)*(b[24] + Tinv*b[25])*(7 - 2*gamma*pow(rho, 2))*calc_exp
        + pow(rho, 8)*pow(Tinv, 2)*(b[26] + pow(Tinv, 2)*b[27])*(9 - 2*gamma*pow(rho, 2))*calc_exp
        + pow(rho, 10)*pow(Tinv, 2)*(b[28] + Tinv*b[29])*(11 - 2*gamma*pow(rho, 2))*calc_exp
        + pow(rho, 12)*pow(Tinv, 2)*(b[30] + Tinv*(b[31] + Tinv*b[32]))*(13 - 2*gamma*pow(rho, 2))*calc_exp;
    return dpdrho*1e5*1e-3/M; // Convert bar.L/mol to Pa.m^3/kg (= J/kg)
}

double MBWR_temperature(double rho, double p, double rho_c, double M, const std::vector<double> &b, int &status)
{
    // Evaluate temperature from density and pressure. Must solve EOS iteratively.
    // Temperature from ideal gas EOS is a good first guess.
    double T_guess = p / (rho * PC_R_u / M);
    int i = 0;
    double T_temp;
    status = SUCCESS;
    do { // Newtons method to solve for temperature.
        if (i >= MBWR_MAX_ITERATIONS) {
            status = ITERATION_ERROR;
            cout << "MBWR_temperature():\n";
            cout << "    Iterations did not converge.\n";
            return T_guess;
        }
        T_temp = T_guess;
        T_guess -= (MBWR_pressure(rho, T_guess, rho_c, M, b) - p)/MBWR_dpdT(rho, T_guess, rho_c, M, b);
        i++;
    } while (fabs(T_temp - T_guess) > MBWR_PROPS_CONVERGE_TOL);
    return T_guess;
}

double MBWR_density(double T, double p, double rho_c, double M, const std::vector<double> &b, int &status)
{
    // Evaluate density from temperature and pressure. Must solve EOS iteratively.
    // Density from ideal gas EOS is a good first guess.
    double rho_guess = p * M / (PC_R_u * T);
    int i = 0;
    double rho_temp;
    status = SUCCESS;
    do { // Newtons method to solve for density.
        if (i >= MBWR_MAX_ITERATIONS) {
            status = ITERATION_ERROR;
            cout << "MBWR_density():\n";
            cout << "    Iterations did not converge.\n";
            return rho_guess;
        }
        rho_temp = rho_guess;
        rho_guess -= (MBWR_pressure(rho_guess, T, rho_c, M, b) - p)/MBWR_dpdrho(rho_guess, T, rho_c, M, b);
        i++;
    } while (fabs(rho_temp - rho_guess) > MBWR_PROPS_CONVERGE_TOL);
    return rho_guess;
}

double MBWR_integral_const_T_energy(double rho, double T, double rho_c, double M, const std::vector<double> &b)
{
    // Integral of (rho^-2)*(P - T(dP/dT)) from zero to current density, w.r.t density.
    // For use in the evaluation of internal energy.
    rho = rho*1e-3/M; // Density from kg/m^3 to mol/L
    rho_c = rho_c*1e-3/M; // Critical density from kg/m^3 to mol/L
    double Tinv = 1/T;
    double Tinvs = pow(Tinv, 2);
    double gamma = 1/pow(rho_c, 2);
    double calc_exp = exp(-gamma*pow(rho, 2));
    // Simple polynomial terms
    double integral = rho*(0.5*b[2]*sqrt(T) + b[3] + Tinv*(2*b[4] + 3*Tinv*b[5]))
        + pow(rho, 2)*(b[7] + Tinv*(2*b[8] + 3*Tinv*b[9]))/2
        + pow(rho, 3)*(b[11] + 2*Tinv*b[12])/3
        + pow(rho, 4)*b[13]/4
        + pow(rho, 5)*Tinv*(2*b[14] + 3*Tinv*b[15])/5
        + pow(rho, 6)*2*b[16]*Tinv/6
        + pow(rho, 7)*Tinv*(2*b[17] + 3*Tinv*b[18])/7
        + pow(rho, 8)*Tinvs*3*b[19]/8;
    // More complex integrated terms
    double A08 = Tinvs*(3*b[20] + 4*Tinv*b[21]);
    double A09 = Tinvs*(3*b[22] + 5*Tinvs*b[23]);
    double A10 = Tinvs*(3*b[24] + 4*Tinv*b[25]);
    double A11 = Tinvs*(3*b[26] + 5*Tinvs*b[27]);
    double A12 = Tinvs*(3*b[28] + 4*Tinv*b[29]);
    double A13 = Tinvs*(3*b[30] + Tinv*(4*b[31] + 5*Tinv*b[32]));
    double rhos = pow(rho, 2);
    double gamin = 1/gamma;
    double gamins = pow(gamin, 2);
    double gaminc = gamin*gamins;
    integral += A08*gamin*(1 - calc_exp)/2
        + A09*gamin*(gamin - calc_exp*(rhos + gamin))/2
        + A10*gamin*(gamins - calc_exp*(rhos*(gamin + rhos/2) + gamins))
        + A11*gamin*(3*gaminc - calc_exp*(rhos*(3*gamins + rhos*(3*gamin/2 + rhos/2)) + 3*gaminc))
        + A12*gamin*(12*gamins*gamins - calc_exp*(rhos*(12*gaminc + rhos*(6*gamins + rhos*(2*gamin + rhos/2))) + 12*gamins*gamins))
        + A13*gamin*(60*gaminc*gamins - calc_exp*(rhos*(60*gamins*gamins + rhos*(30*gaminc + rhos*(10*gamins + rhos*(5*gamin/2 + rhos/2)))) + 60*gaminc*gamins));
    return integral*1e5*1e-3/M; // Convert bar.L/mol to Pa.m^3/kg (= J/kg)
}

double MBWR_dintegral_const_T_energy_dT(double rho, double T, double rho_c, double M, const std::vector<double> &b)
{
    // Derivative of MBWR_integral_const_T_energy w.r.t temperature.
    // For use in the evaluation of Cv.
    rho = rho*1e-3/M; // Density from kg/m^3 to mol/L
    rho_c = rho_c*1e-3/M; // Critical density from kg/m^3 to mol/L
    double Tinv = 1/T;
    double Tinvs = pow(Tinv, 2);
    double Tinvc = Tinv*Tinvs;
    double gamma = 1/pow(rho_c, 2);
    double calc_exp = exp(-gamma*pow(rho, 2));
    // Simple polynomial terms
    double dintegral_dT = rho*(0.25*b[2]/sqrt(T) - Tinvs*(2*b[4] + 6*Tinv*b[5]))
        - pow(rho, 2)*(Tinvs*(2*b[8] + 6*Tinv*b[9]))/2
        - pow(rho, 3)*2*Tinvs*b[12]/3
        - pow(rho, 5)*Tinvs*(2*b[14] + 6*Tinv*b[15])/5
        - pow(rho, 6)*2*b[16]*Tinvs/6
        - pow(rho, 7)*Tinvs*(2*b[17] + 6*Tinv*b[18])/7
        - pow(rho, 8)*Tinvc*6*b[19]/8;
    // More complex integrated terms
    double A08 = -Tinvc*(6*b[20] + 12*Tinv*b[21]);
    double A09 = -Tinvc*(6*b[22] + 20*Tinvs*b[23]);
    double A10 = -Tinvc*(6*b[24] + 12*Tinv*b[25]);
    double A11 = -Tinvc*(6*b[26] + 20*Tinvs*b[27]);
    double A12 = -Tinvc*(6*b[28] + 12*Tinv*b[29]);
    double A13 = -Tinvc*(6*b[30] + Tinv*(12*b[31] + 20*Tinv*b[32]));
    double rhos = pow(rho, 2);
    double gamin = 1/gamma;
    double gamins = pow(gamin, 2);
    double gaminc = gamin*gamins;
    dintegral_dT += A08*gamin*(1 - calc_exp)/2
        + A09*gamin*(gamin - calc_exp*(rhos + gamin))/2
        + A10*gamin*(gamins - calc_exp*(rhos*(gamin + rhos/2) + gamins))
        + A11*gamin*(3*gaminc - calc_exp*(rhos*(3*gamins + rhos*(3*gamin/2 + rhos/2)) + 3*gaminc))
        + A12*gamin*(12*gamins*gamins - calc_exp*(rhos*(12*gaminc + rhos*(6*gamins + rhos*(2*gamin + rhos/2))) + 12*gamins*gamins))
        + A13*gamin*(60*gaminc*gamins - calc_exp*(rhos*(60*gamins*gamins + rhos*(30*gaminc + rhos*(10*gamins + rhos*(5*gamin/2 + rhos/2)))) + 60*gaminc*gamins));
    return dintegral_dT*1e5*1e-3/M; // Convert bar.L/mol to Pa.m^3/kg (= J/kg)
}

double MBWR_integral_const_T_entropy(double rho, double T, double rho_c, double M, const std::vector<double> &b)
{
    // Integral of (rho^-2)*(rho*R - (dP/dT)) from zero to current density, w.r.t density.
    // For use in the evaluation of entropy.
    rho = rho*1e-3/M; // Density from kg/m^3 to mol/L
    rho_c = rho_c*1e-3/M; // Critical density from kg/m^3 to mol/L
    double Tinv = 1/T;
    double Tinvs = pow(Tinv, 2);
    double Tinvc = Tinv*Tinvs;
    double gamma = 1/pow(rho_c, 2);
    double calc_exp = exp(-gamma*pow(rho, 2));
    // Simple polynomial terms
    double integral = - rho*(b[1] + 0.5*b[2]/sqrt(T) - Tinvs*(b[4] + 2*Tinv*b[5]))
        - pow(rho, 2)*(b[6] - Tinvs*(b[8] + 2*Tinv*b[9]))/2
        - pow(rho, 3)*(b[10] - Tinvs*b[12])/3
        + pow(rho, 5)*Tinvs*(b[14] + 2*Tinv*b[15])/5
        + pow(rho, 6)*Tinvs*b[16]/6
        + pow(rho, 7)*Tinvs*(b[17] + 2*Tinv*b[18])/7
        + pow(rho, 8)*2*Tinvc*b[19]/8;
    // More complex integrated terms
    double B07 = Tinvc*(2*b[20] + 3*Tinv*b[21]);
    double B08 = Tinvc*(2*b[22] + 4*Tinvs*b[23]);
    double B09 = Tinvc*(2*b[24] + 3*Tinv*b[25]);
    double B10 = Tinvc*(2*b[26] + 4*Tinvs*b[27]);
    double B11 = Tinvc*(2*b[28] + 3*Tinv*b[29]);
    double B12 = Tinvc*(2*b[30] + Tinv*(3*b[31] + 4*Tinv*b[32]));
    double rhos = pow(rho, 2);
    double gamin = 1/gamma;
    double gamins = pow(gamin, 2);
    double gaminc = gamin*gamins;
    integral += B07*gamin*(1 - calc_exp)/2
        + B08*gamin*(gamin - calc_exp*(rhos + gamin))/2
        + B09*gamin*(gamins - calc_exp*(rhos*(gamin + rhos/2) + gamins))
        + B10*gamin*(3*gaminc - calc_exp*(rhos*(3*gamins + rhos*(3*gamin/2 + rhos/2)) + 3*gaminc))
        + B11*gamin*(12*gamins*gamins - calc_exp*(rhos*(12*gaminc + rhos*(6*gamins + rhos*(2*gamin + rhos/2))) + 12*gamins*gamins))
        + B12*gamin*(60*gaminc*gamins - calc_exp*(rhos*(60*gamins*gamins + rhos*(30*gaminc + rhos*(10*gamins + rhos*(5*gamin/2 + rhos/2)))) + 60*gaminc*gamins));
    return integral*1e5*1e-3/M; // Convert bar.L/(mol.K) to Pa.m^3/(kg.K) (= J/(kg.K))
}
