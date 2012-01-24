// Author: PJ
// Date: Dec 2007 -- Jan 2008

#include <iostream>
#include <sstream>
#include <cmath>
#include "kernel.hh"
#include "c-flow-condition.hh"

using namespace std;

//---------------------------------------------------------------------- 
// Class for defining flow conditions, mainly for use in
// the new boundary condition classes within the CFD codes.
// Because of this intended use, it uses the managed gas model
// rather than a GasModel object directly. 
//
// PJ, Dec 2007 -- Jan 2008, first code

CFlowCondition::CFlowCondition( Gas_model *gmodel,
				double p, 
				double u, double v, double w,
				const std::vector<double> T,
				const std::vector<double> massf,
				const std::string label,
				double tke, double omega,
				double mu_t, double k_t,
				int S )
    : u(u), v(v), w(w), label(label),
      tke(tke), omega(omega), mu_t(mu_t), k_t(k_t), S(S)
{
    gas = new Gas_data(gmodel);
    gas->p = p;
    int nmodes = gmodel->get_number_of_modes();
    if ( nmodes != (int)T.size() ) {
	cerr << "CFlowCondition(): nmodes=" << nmodes 
	     << ", size of T=" << T.size() << endl;
    }
    gas->T = T;
    int nsp = gmodel->get_number_of_species();
    if ( nsp != (int)massf.size() ) {
	cerr << "CFlowCondition(): nsp=" << nsp 
	     << ", size of massf=" << massf.size() << endl;
    }
    double sum_massf = 0.0;
    for ( int isp = 0; isp < nsp; ++isp ) sum_massf += massf[isp];
    if ( fabs(sum_massf - 1.0) > 1.0e-5 ) {
	cerr << "CFlowCondition() warning: mass fractions sum to " 
	     << sum_massf << endl;
    }
    gas->massf = massf;
    gmodel->eval_thermo_state_pT(*gas);
    gmodel->eval_transport_coefficients(*gas);
}

CFlowCondition::CFlowCondition( const CFlowCondition &cfc ) 
    : u(cfc.u), v(cfc.v), w(cfc.w),
      label(cfc.label), tke(cfc.tke), omega(cfc.omega),
      mu_t(cfc.mu_t), k_t(cfc.k_t), S(cfc.S)
{
    Gas_model *gmodel = get_gas_model_ptr();
    gas = new Gas_data(gmodel);
    gas->copy_values_from(*(cfc.gas));
    // The following two calls may be redundant since we have just
    // made a complete copy of the gas_data.
    gmodel->eval_thermo_state_pT(*gas);
    gmodel->eval_transport_coefficients(*gas);
}

CFlowCondition::~CFlowCondition() 
{
    delete gas;
}

/// \brief CFlowCondition compact string representation
///        as would be usable in a Python script.
string CFlowCondition::str() const
{
    ostringstream ost;
    ost << "CFlowCondition( p=" << gas->p
	<< ", u=" << u << ", v=" << v << ", w=" << w;
    ost << ", T=[";
    for ( size_t i = 0; i < gas->T.size(); ++i ) {
	ost << gas->T[i] << ",";
    }
    ost << "]"; // end of temperatures
    ost << ", massf=[";
    for ( size_t i = 0; i < gas->massf.size(); ++i ) {
	ost << gas->massf[i] << ",";
    }
    ost << "]"; // end of mass fractions
    ost << ", label=\"" << label << "\"";
    ost << ", mu_t=" << mu_t << ", k_t=" << k_t
	<< ", tke=" << tke << ", omega=" << omega;
    ost << ", S=" << S;
    ost << " )";
    return ost.str();
}

string CFlowCondition::write_to_ini_str( int indx ) const
// Write the data to a string using INI-file format.
{
    ostringstream ost;
    ost.setf(ios_base::scientific, ios_base::floatfield);
    ost << "[flow/" << indx << "]" << endl;
    ost << "label = " << label << endl;
    ost << "p = " << gas->p << endl;
    ost << "nmodes = " << gas->T.size() << endl;
    ost << "T = ";
    for ( size_t i = 0; i < gas->T.size(); ++i ) {
	ost << gas->T[i] << " ";
    }
    ost << endl;
    ost << "nsp = " << gas->massf.size() << endl;
    ost << "massf = ";
    for ( size_t i = 0; i < gas->massf.size(); ++i ) {
	ost << gas->massf[i] << " ";
    }
    ost << endl;
    ost << "u = " << u << endl;
    ost << "v = " << v << endl;
    ost << "w = " << w << endl;
    ost << "tke = " << tke << endl;
    ost << "omega = " << omega << endl;
    ost << "mu_t = " << mu_t << endl;
    ost << "k_t = " << k_t << endl;
    ost << "S = " << S << endl;
    return ost.str();
}


/// \brief Overload stream output for CFlowCondition objects
ostream& operator<<( ostream &os, const CFlowCondition &cfc )
{
    os << cfc.str();
    return os;
}

