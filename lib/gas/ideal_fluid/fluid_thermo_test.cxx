// fluid_thermo_test.cxx

#include <iostream>
#include <string>
#include "fluid_thermo.hh"
using std::cout;
using std::endl;

int main ()
{
    cout << "fluid_thermo_test:" << endl;

    PureFluid *co2 = new CarbonDioxide();
    cout << co2->name << endl;

    cout << "Pressure from density and temperature:" << endl;
    double rho = 1.0/0.00779;
    double T = 300.0;
    double dT = 0.01;
    cout << "rho=" << rho << " T=" << T << " p=" << co2->p_rhoT(rho,T);
    cout << " expected p=5.0e6" << endl;
    cout << "    dpdT=" << co2->dpdT_rhoT(rho,T) 
	 << " finite-diff=" << 0.5*(co2->p_rhoT(rho,T+dT) - co2->p_rhoT(rho,T-dT))/dT
	 << " rhoR=" << rho*co2->R << endl;
    double u, h, s;
    co2->us_rhoT(rho, T, u, h, s);
    cout << "    h=" << h << " u=" << u << " s=" << s << " expected h=366980. s=1335.1" << endl;

    rho = 1.0/0.01824;
    T = 500.0;
    cout << "rho=" << rho << " T=" << T << " p=" << co2->p_rhoT(rho,T);
    cout << " expected p=5.0e6" << endl;
    cout << "    dpdT=" << co2->dpdT_rhoT(rho,T) 
	 << " finite-diff=" << 0.5*(co2->p_rhoT(rho,T+dT) - co2->p_rhoT(rho,T-dT))/dT
	 << " rhoR=" << rho*co2->R << endl;
    co2->us_rhoT(rho, T, u, h, s);
    cout << "    h=" << h << " u=" << u << " s=" << s << " expected h=600430. s=1940.7" << endl;

    rho = 1.0/0.03039;
    T = 800.0;
    cout << "rho=" << rho << " T=" << T << " p=" << co2->p_rhoT(rho,T);
    cout << " expected p=5.0e6" << endl;
    cout << "    dpdT=" << co2->dpdT_rhoT(rho,T) 
	 << " finite-diff=" << 0.5*(co2->p_rhoT(rho,T+dT) - co2->p_rhoT(rho,T-dT))/dT
	 << " rhoR=" << rho*co2->R << endl;
    co2->us_rhoT(rho, T, u, h, s);
    cout << "    h=" << h << " u=" << u << " s=" << s << " expected h=941100. s=2472.4" << endl;

    rho = 1.0/0.1893;
    T = 1000.0;
    cout << "rho=" << rho << " T=" << T << " p=" << co2->p_rhoT(rho,T);
    cout << " expected p=1.0e6" << endl;
    cout << "    dpdT=" << co2->dpdT_rhoT(rho,T) 
	 << " finite-diff=" << 0.5*(co2->p_rhoT(rho,T+dT) - co2->p_rhoT(rho,T-dT))/dT
	 << " rhoR=" << rho*co2->R << endl;
    co2->us_rhoT(rho, T, u, h, s);
    cout << "    h=" << h << " u=" << u << " s=" << s << " expected h=1186610. s=3051.5" << endl;

    cout << "Try utility function..." << endl;
    State st;
    st.p = 4.9993e6;
    st.T = 800.0;
    co2->eval_state(&st, "pt");
    cout << "PT: rho=" << st.rho << " T=" << st.T << " st.p=" << st.p << endl;
    st.s = 2472.44;
    st.T = 800.0;
    co2->eval_state(&st, "st");
    cout << "ST: rho=" << st.rho << " T=" << st.T << " s=" << st.s << endl;
    st.rho = 32.9056;
    st.u = 789170.0;
    co2->eval_state(&st, "du");
    cout << "DU: rho=" << st.rho << " T=" << st.T << " u=" << st.u << endl;
    st.p = 4.9993e6;
    st.h = 941098.0;
    co2->eval_state(&st, "ph");
    cout << "PH: rho=" << st.rho << " T=" << st.T << " p=" << st.p << " h=" << st.h << endl;
    st.p = 4.9993e6;
    st.s = 2472.44;
    co2->eval_state(&st, "ps");
    cout << "PS: rho=" << st.rho << " T=" << st.T << " p=" << st.p << " s=" << st.s << endl;

    cout << "Saturated liquid density:" << endl;
    T = 280.0;
    cout << "T=" << T << " rho_f=" << co2->sat_liquid_density(T);
    cout << " expected rho_f=895.48" << endl;

    cout << "Saturated pressure:" << endl;
    T = 287.5;
    rho = 1.0/0.00641;
    cout << "T=" << T << " p_sat=" << co2->p_sat(T);
    cout << " expected p_sat=5.0e6" << endl;

    cout << "Specific heat, constant volume, ideal gas:" << endl;
    T = 287.5;
    cout << "T=" << T << " cv0=" << co2->cv0(T) << endl;

    cout << "Done." << endl;
    return 0;
}
