// Author: Rowan J. Gollan
// Date: 06-Nov-2008
// Place: Poquoson, Virginia, USA
// Note:
//   This is a port of PJs look-up table
//   implementation with some cosmetic changes:
//     - gzipped text files, instead of binary format
//     - vector storage instead of arrays
//     - reworked to fit in new class framework
//

#ifndef LOOK_UP_TABLE_HH
#define LOOK_UP_TABLE_HH

#include <string>

#include "gas_data.hh"
#include "gas-model.hh"

class Look_up_table : public Gas_model {
public:
    Look_up_table(std::string cfile);
    ~Look_up_table();
private:
    int determine_interpolants(const Gas_data &Q, int &ir, int &ie,
			       double &lrfrac, double &efrac);
    int s_eval_thermo_state_rhoe(Gas_data &Q);
    int s_eval_transport_coefficients(Gas_data &Q);
    int s_eval_diffusion_coefficients(Gas_data &Q);
    double s_molecular_weight(int isp); 
    double s_internal_energy(const Gas_data &Q, int isp);
    double s_enthalpy(const Gas_data &Q, int isp);
    double s_entropy(const Gas_data &Q, int isp);
    double s_dedT_const_v(const Gas_data &Q, int &status);
    double s_dhdT_const_p(const Gas_data &Q, int &status);
    double s_gas_constant(const Gas_data &Q, int &status);

    int iesteps_, irsteps_;
    double emin_, emax_, de_;
    double lrmin_, lrmax_, dlr_;

    matrix Cv_hat_;
    matrix Cv_;
    matrix R_hat_;
    matrix g_hat_;
    matrix mu_hat_;
    matrix k_hat_;
};

Gas_model* create_look_up_table_gas_model(std::string cfile);

#endif
