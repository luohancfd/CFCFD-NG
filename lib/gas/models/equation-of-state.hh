// Author: Rowan J. Gollan
// Version: 24-May-2008
//            Initial coding.
//

#ifndef EQUATION_OF_STATE_HH
#define EQUATION_OF_STATE_HH

#include "../../util/source/useful.h"
#include "gas_data.hh"
#include "gas-model.hh"

class Equation_of_state {
public:
    Equation_of_state() {}
    virtual ~Equation_of_state() {}

    int eval_pressure(Gas_data &Q)
    { return s_eval_pressure(Q); }

    int eval_temperature(Gas_data &Q)
    { return s_eval_temperature(Q); }

    int eval_density(Gas_data &Q) 
    { return s_eval_density(Q); }

    // For evaluating internal energy from enthalpy
    double prho_ratio(const Gas_data &Q, int isp) 
    { return s_prho_ratio(Q, isp); }

    double dTdp_const_rho(const Gas_data &Q, int &status)
    { return s_dTdp_const_rho(Q, status); }

    double dTdrho_const_p(const Gas_data &Q, int &status)
    { return s_dTdrho_const_p(Q, status); }

    double dpdrho_const_T(const Gas_data &Q, int &status)
    { return s_dpdrho_const_T(Q, status); }
    
    double dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
    { return s_dpdrho_i_const_T(Q, isp, status); }

    double dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
    { return s_dpdT_i_const_rho(Q, itm, status); }
    
    double gas_constant(const Gas_data &Q, int &status)
    { status = SUCCESS; return calculate_gas_constant(Q.massf, M_); }

    double molecular_weight(const Gas_data &Q)
    { return calculate_molecular_weight(Q.massf, M_); }

protected:
    std::vector<double> M_;

private:
    virtual int s_eval_pressure(Gas_data &Q) = 0;
    virtual int s_eval_temperature(Gas_data &Q) = 0;
    virtual int s_eval_density(Gas_data &Q) = 0;
    virtual double s_prho_ratio(const Gas_data &Q, int isp) = 0;
    virtual double s_dTdp_const_rho(const Gas_data &Q, int &status) = 0;
    virtual double s_dTdrho_const_p(const Gas_data &Q, int &status) = 0;
    virtual double s_dpdrho_const_T(const Gas_data &Q, int &status) = 0;
    virtual double s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status) = 0;
    virtual double s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status) = 0;
    virtual double s_gas_constant(const Gas_data &Q, int &status) = 0;
};

#endif
