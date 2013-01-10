// Author: Rowan J. Gollan
// Version: 24-May-2008
//            Initial coding.
//

#ifndef THERMAL_BEHAVIOUR_MODEL_HH
#define THERMAL_BEHAVIOUR_MODEL_HH

#include "physical_constants.hh"
#include "gas_data.hh"
#include "equation-of-state.hh"

#include "../../nm/source/segmented-functor.hh"

class Thermal_behaviour_model {
public:
    Thermal_behaviour_model() {}
    virtual ~Thermal_behaviour_model() {}

    virtual int get_number_of_modes()
    { return 1; }

    int decode_conserved_energy(Gas_data &Q, const std::vector<double> &rhoe)
    { return s_decode_conserved_energy(Q, rhoe); }

    int encode_conserved_energy(const Gas_data &Q, std::vector<double> &rhoe)
    { return s_encode_conserved_energy(Q, rhoe); }

    double dhdT_const_p(const Gas_data &Q, Equation_of_state *EOS_, int &status)
    { return s_dhdT_const_p(Q, EOS_, status); }

    double dedT_const_v(const Gas_data &Q, Equation_of_state *EOS_, int &status)
    { return s_dedT_const_v(Q, EOS_, status); }

    int eval_energy(Gas_data &Q, Equation_of_state *EOS_)
    { return s_eval_energy(Q, EOS_); }

    int eval_temperature(Gas_data &Q, Equation_of_state *EOS_)
    { return s_eval_temperature(Q, EOS_); }

    double eval_energy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
    { return s_eval_energy_isp(Q, EOS_, isp); }

    double eval_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
    { return s_eval_enthalpy_isp(Q, EOS_, isp); }

    double eval_entropy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
    { return s_eval_entropy_isp(Q, EOS_, isp); }
    
    double eval_modal_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp, int itm)
    { return s_eval_modal_enthalpy_isp(Q, EOS_, isp, itm); }
    
    double eval_modal_Cv(Gas_data &Q, Equation_of_state *EOS_, int itm)
    { return s_eval_modal_Cv(Q, EOS_, itm); }

    double eval_modal_massf(const Gas_data &Q, int itm)
    {return s_eval_modal_massf(Q, itm); }

private:
    virtual int s_decode_conserved_energy(Gas_data &Q, const std::vector<double> &rhoe) = 0;
    virtual int s_encode_conserved_energy(const Gas_data &Q, std::vector<double> &rhoe) = 0;
    virtual double s_dhdT_const_p(const Gas_data &Q, Equation_of_state *EOS_, int &status) = 0;
    virtual double s_dedT_const_v(const Gas_data &Q, Equation_of_state *EOS_, int &status) = 0;
    virtual int s_eval_energy(Gas_data &Q, Equation_of_state *EOS_) = 0;
    virtual int s_eval_temperature(Gas_data &Q, Equation_of_state *EOS_) = 0;
    virtual double s_eval_energy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp) = 0;
    virtual double s_eval_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp) = 0;
    virtual double s_eval_entropy_isp(const Gas_data &Q, Equation_of_state *, int isp) = 0;
    virtual double s_eval_modal_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp, int itm );
    virtual double s_eval_modal_Cv(Gas_data &Q, Equation_of_state *EOS_, int itm);
    virtual double s_eval_modal_massf(const Gas_data &Q, int itm);
};

int tbm_decode_conserved_energy(std::vector<double> &e, 
				const std::vector<double> &rhoe, 
				const double &rho);
int tbm_encode_conserved_energy(std::vector<double> &rhoe, 
				const std::vector<double> &e, 
				const double &rho);
double tbm_dhdT_const_p(const std::vector<Segmented_functor *> &Cp_, 
			const std::vector<double> &massf, 
			const std::vector<double> &T);

#endif
