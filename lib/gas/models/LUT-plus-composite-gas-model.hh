// Author: PJ, DFP, RJG (see below)
// Date: 05-Dec-2008
// Place: Hampton, Virginia, USA
// Note:
//   This is a port of PJs implementation of the 
//   LUT + ideal gas mix model.  DFP updated this
//   previously for libgas2.  This implementation
//   fits with the new Gas_model class.
//

#ifndef LUT_PLUS_COMPOSITE_GAS_MODEL_HH
#define LUT_PLUS_COMPOSITE_GAS_MODEL_HH

#include <string>

#include "gas_data.hh"
#include "gas-model.hh"
#include "look-up-table.hh"
#include "composite-gas-model.hh"


class LUT_plus_composite : public Gas_model {
public:
    LUT_plus_composite(const std::string cfile);
    ~LUT_plus_composite();
private:
    Look_up_table *LUT_;
    Composite_gas_model *CGM_;
    Gas_data *Q_LUT_;
    Gas_data *Q_CGM_;

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
};

Gas_model* create_LUT_plus_composite_gas_model(const std::string cfile);

#endif
