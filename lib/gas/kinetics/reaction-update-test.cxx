// Author: Rowan J. Gollan
// Date: 30-Oct-2008
// Place: Hampton, Viriginia, USA
//

#include <string>
#include <fstream>

#include "../models/physical_constants.hh"
#include "../models/gas-model.hh"
#include "reaction-update.hh"

using namespace std;

int main()
{
    string gasfile("reaction-update-test-gas-model.lua");
    string reacfile("reaction-update-test-reactions.lua");

    Gas_model *g = create_gas_model(gasfile);
    Reaction_update *r = create_Reaction_update(string(reacfile), *g);
    
    Gas_data Q(g);
    vector<double> c, M;
    c.resize(3, 0.0);
    M.resize(3, 0.0);
    for ( size_t isp = 0; isp < 3; ++isp ) {
	M[isp] = g->molecular_weight(isp);
    }

    // Problem setup
    // -- constants
    double M_H2 = g->molecular_weight(0);
    double M_I2 = g->molecular_weight(1);
    
    // -- initial conditions
    double c_0 = 4.54;
    double T_0 = 700.0;
    double t_final = 100000.0;
    double t_inc = 100.0;
    double dt_suggest = 1.0;

    // -- derived values
    double p_0 = 2.0*c_0*PC_R_u*T_0;
    double M_mix = 0.5*M_H2 + 0.5*M_I2;
    double R_mix = PC_R_u_kmol / M_mix;
    double rho_0 = p_0 / (R_mix*T_0);
    double f_H2 = 0.5 * (M_H2/M_mix);
    double f_I2 = 0.5 * (M_I2/M_mix);

    Q.T[0] = 700.0;
    Q.p = p_0;
    Q.massf[0] = f_H2;
    Q.massf[1] = f_I2;

    g->eval_thermo_state_pT(Q);

    ofstream ofp;
    ofp.open("reaction-update-test.data");
    
    double t = 0.0;
    while ( t <= t_final ) {
    	r->update_state(Q, t_inc, dt_suggest);
	g->eval_thermo_state_rhoT(Q);
	t += t_inc;
	convert_massf2conc(Q.rho, Q.massf, M, c);
	ofp << t << " " << c[2] << endl;
    }
    ofp.close();

    delete g;
    delete r;

    return 0;
    
}
