// Author: PJ
// Date: Dec 2007 -- Jan 2008

#ifndef C_FLOW_CONDITION_HH
#define C_FLOW_CONDITION_HH

#include <sstream>
#include <vector>
#include <string>

#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"

/// \brief CFlowCondition class for use in boundary conditions.
///
/// PJ, Dec 2007 - Jan 2008, sufficient code to support new BC objects

class CFlowCondition {
public:
    Gas_data *gas;
    double u, v, w;
    double Bx, By, Bz, psi, divB;
    std::string label;
    double tke, omega;
    double mu_t, k_t;
    int S; // shock indicator

public:
    CFlowCondition( Gas_model *gmodel,
		    double p, 
		    double u, double v, double w, 
		    const std::vector<double> T,
		    const std::vector<double> massf,
		    const std::string label="",
		    double tke=0.0, double omega=1.0,
		    double mu_t=0.0, double k_t=0.0,
		    int S=0,
		    double Bx=0.0, double By=0.0, double Bz=0.0,
		    double psi=0.0, double divB=0.0);
    CFlowCondition( const CFlowCondition &cfc );
    CFlowCondition();
    CFlowCondition & operator=(const CFlowCondition &cfc);
    ~CFlowCondition();
    /// \brief Returns a string representation of the CFlowCondition.
    std::string str() const;
    // void print_data() const;
    std::string write_to_ini_str( int indx ) const;
    std::string write_to_json_str( int indx ) const;
    // std::string write_to_mbcns_p_str( int indx ) const;
};

// Helper functions.

/// \brief Output a string representation of the generic gas model.
std::ostream& operator<<( std::ostream &os, const CFlowCondition &cfc );

#endif
