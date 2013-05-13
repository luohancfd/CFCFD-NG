// Author: Rowan J. Gollan
// Date: 12-Sep-2008

#ifndef REACTION_HH
#define REACTION_HH

#include <string>
#include <valarray>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "../models/gas-model.hh"
#include "reaction-rate-coeff.hh"
#include "equilibrium-constant.hh"

class Reaction {
public:
    Reaction(lua_State *L, Gas_model &g);
    virtual ~Reaction();
    
    void compute_k_f(const Gas_data &Q)
    { k_f_ = s_compute_k_f(Q); }

    double k_f()
    { return k_f_; }

    void compute_k_b(const Gas_data &Q)
    { k_b_ = s_compute_k_b(Q); }

    double k_b()
    { return k_b_; }

    void compute_forward_rate(const std::valarray<double> &y)
    { w_f_ = s_compute_forward_rate(y); }
    double w_f()
    { return w_f_; }

    void compute_backward_rate(const std::valarray<double> &y)
    { w_b_ = s_compute_backward_rate(y); }
    double w_b()
    { return w_b_; }
    
    double production(int isp);

    double loss(int isp);

    bool compute_kf_first(void)
    { return compute_kf_first_; }

    int get_nu(int isp);

protected:
    virtual double s_compute_k_f(const Gas_data &Q);
    virtual double s_compute_k_b(const Gas_data &Q);
    virtual double s_compute_forward_rate(const std::valarray<double> &y) = 0;
    virtual double s_compute_backward_rate(const std::valarray<double> &y) = 0;
    
private:
    double w_f_;
    double w_b_;
    double k_f_;
    double k_b_;
    std::map<int, int> nu_;
    bool compute_kf_first_;
    bool different_rcts_;
    Reaction_rate_coefficient *frc_;
    Reaction_rate_coefficient *brc_;
    Equilibrium_constant *ec_;
    Gas_data *Q_;
};

Reaction* create_Reaction(lua_State *L, Gas_model &g);
Reaction* get_reaction_from_file(int ir, std::string cfile, Gas_model &g);
#endif
