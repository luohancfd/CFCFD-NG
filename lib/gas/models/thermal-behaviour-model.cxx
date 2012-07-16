// brief: Contains common functions for thermal behaviour models
// version: 22 Sep - initial copying
//

#include "../../util/source/useful.h"

#include "thermal-behaviour-model.hh"

using namespace std;

int
tbm_decode_conserved_energy(vector<double> &e, 
			    const vector<double> &rhoe, 
			    const double &rho)
{    
    // The assumption is that use of this model
    // implies only one thermal mode.
    e[0] = rhoe[0]/rho;
    return SUCCESS;
}

int
tbm_encode_conserved_energy(vector<double> &rhoe, 
			    const vector<double> &e, 
			    const double &rho)
{
    // The assumption is that use of this model
    // implies only one thermal mode.
    rhoe[0] = rho*e[0];
    return SUCCESS;
}

double
tbm_dhdT_const_p(const vector<Segmented_functor *> &Cp_, 
		 const vector<double> &massf,
		 const vector<double> &T)
{
    // Reference:
    // Cengel and Boles (2002)
    // Thermodynamics: an Engineering Approach, 3rd edition
    // 4th Ed.
    // McGraw Hill
    // Equation 11-35 on p. 615

    // The assumption is that use of this model
    // implies only one thermal mode.
    
    double Cp = 0.0;
    for ( size_t isp = 0; isp < Cp_.size(); ++isp ) {
	Cp += massf[isp]*((*Cp_[isp])(T[0]));
    }
    return Cp;
}

double
Thermal_behaviour_model::
s_eval_modal_enthalpy_isp( const Gas_data &Q, Equation_of_state *EOS_, int isp, int itm )
{
    // NOTE: this function should never be used as Noneq_thermal_behaviour
    //       implements its own version
    cout << "Thermal_behaviour_model::s_eval_modal_enthalpy_isp()" << endl
         << "Something has gone wrong as this function is not intended for use." << endl;
    exit( FAILURE );
}


double
Thermal_behaviour_model::
s_eval_modal_Cv(Gas_data &Q, Equation_of_state *EOS_, int itm)
{
    // NOTE: all thermal behaviour models except Noneq_thermal_behaviour should implement this function
    int status;
    return s_dedT_const_v(Q,EOS_,status);
}
