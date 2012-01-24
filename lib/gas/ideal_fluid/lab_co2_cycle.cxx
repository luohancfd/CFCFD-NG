// lab_co2_cycle.cxx
// Supercritical CO2 Brayton cycle for QEGC lab.
//
// Modelled on Hal's MATLAB code (co2_co2.m)
// PJ Apr,May 2008
//
// Output is to stdout but is commented so that it is compatible with GnuPlot.

#include <vector>
#include <iostream>
#include <string>
#include "fluid_thermo.hh"
using std::cout;
using std::endl;

int main()
{
    double Treservoir = 273.0 + 235.0; // hot rock temperature
    double Preservoir = 25.01e6;       // pressure at the bottom of the bore
                                       // as computed by Hal's matlab code
    double Tcooler = 320.0;            // ambient air temperature
    double Pcooler = 8.0e6;            // pressure in heat exchanger on surface

    std::vector<State> st; // set of CO2 states
    int N = 16;
    st.resize(N);
    double eta_c = 0.70; // compressor efficiency
    double eta_t = 0.70; // turbine efficiency
    

    cout << "# Supercritical CO2 Brayton Cycle." << endl;
    cout << "# (using Hal's pressures and temperatures to anchor cycle)" << endl;
    cout << "# Pressure in MPa, T in degrees K, enthalpy in kJ/kg, entropy in kJ/kg.K" << endl;
    cout << "" << endl;
    PureFluid *fluid = new CarbonDioxide();

    cout << "# [0] Inlet to compressor as used by Hal." << endl;
    st[0].p = Pcooler; st[0].T = Tcooler;
    fluid->eval_state(&st[0], "pt");
    cout << "# [0] " << st[0].str() << endl;

    cout << "# [1] Outlet from compressor." << endl;
    st[1].p = 13.39e6; st[1].s = st[0].s; // specified
    st[1].rho = st[0].rho * 1.5; st[1].T = st[0].T * 1.5; // guess
    fluid->eval_state(&st[1], "ps" ,1);
    cout << "# [1] " << st[1].str() << endl;

    cout << "# [11] imperfect compressor: efficiency=" << eta_c << endl;
    st[11].h = 1.0/eta_c * (st[1].h - st[0].h) + st[0].h;
    st[11].p = st[1].p;
    st[11].rho = st[1].rho; st[11].T = st[1].T; // guess
    fluid->eval_state(&st[11], "ph" ,1);
    cout << "# [11] " << st[11].str() << endl;

    cout << "# [2] Bottom of bore (no change in elevation in lab)." << endl;
    st[2].p = st[1].p; st[2].s = st[0].s; // specified
    st[2].rho = st[1].rho; st[2].T = st[1].T; // must be...
    fluid->eval_state(&st[2], "ps" ,1);
    cout << "# [2] " << st[2].str() << endl;

    cout << "# [3] End of hot reservoir, entry to upward bore." << endl;
    st[3].p = 13.39e6; st[3].T = 444.0; // specified from thermosiphon cycle
    st[3].rho = st[2].rho*st[2].T/st[3].T; // guess
    fluid->eval_state(&st[3], "pt" ,1);
    cout << "# [3] " << st[3].str() << endl;

    cout << "# [4] Exit of bore at surface, inlet to turbine." << endl;
    st[4].p = 13.39e6; st[4].s = st[3].s; // specified, isentropic
    st[4].rho = st[3].rho; st[4].T = st[3].T; // must be...
    fluid->eval_state(&st[4], "ps" ,1);
    cout << "# [4] " << st[4].str() << endl;

    cout << "# [5] Exit of turbine." << endl;
    st[5].p = Pcooler; st[5].s = st[4].s; // specified, isentropic
    st[5].rho = st[4].rho/2.0; st[5].T = st[4].T/1.5; // guess
    fluid->eval_state(&st[5], "ps" ,1);
    cout << "# [5] " << st[5].str() << endl;

    cout << "# [15] imperfect turbine: efficiency=" << eta_t << endl;
    st[15].h = st[4].h - eta_t * (st[4].h - st[5].h);
    st[15].p = Pcooler;
    st[15].rho = st[5].rho; st[15].T = st[5].T; // guess
    fluid->eval_state(&st[15], "ph" ,1);
    cout << "# [15] " << st[15].str() << endl;

    cout << "------------------------------" << endl;
    cout << "With ideal turbomachinery:" << endl;
    double work_from_turbine = st[4].h - st[5].h;
    double heat_added = st[3].h - st[2].h;
    // double heat_added = (st[3].s - st[2].s) * 0.5 * (st[3].T + st[2].T); // Hal's approach.
    double thermal_eff = work_from_turbine / heat_added;
    double heat_rejected = st[5].h - st[0].h;
    double work_by_pump = st[1].h - st[0].h;
    double cycle_eff = (work_from_turbine - work_by_pump) / heat_added;
    cout << "# work-from-turbine=" << work_from_turbine/1000.0 << " kJ/kg" << endl;
    cout << "# heat-added=" << heat_added/1000.0 << " kJ/kg" << endl;
    cout << "# geothermal-efficiency=" << thermal_eff << endl;
    cout << "# cycle-thermal-efficiency=" << cycle_eff << endl;
    cout << "# work-done-by-pump=" << work_by_pump/1000.0 << " kJ/kg" << endl;
    cout << "# work-done-down-well=" << (st[2].h-st[1].h)/1000.0 << " kJ/kg" << endl;
    cout << "# work-done-up-well=" << (st[4].h-st[3].h)/1000.0 << " kJ/kg" << endl;
    cout << "# heat-rejected=" << heat_rejected/1000.0 << " kJ/kg" << endl;

    cout << "------------------------------" << endl;
    cout << "With imperfect turbomachinery:" << endl;
    work_from_turbine = st[4].h - st[15].h;
    heat_added = st[3].h - st[11].h;
    thermal_eff = work_from_turbine / heat_added;
    heat_rejected = st[15].h - st[0].h;
    work_by_pump = st[11].h - st[0].h;
    cycle_eff = (work_from_turbine - work_by_pump) / heat_added;
    cout << "# work-from-turbine=" << work_from_turbine/1000.0 << " kJ/kg" << endl;
    cout << "# heat-added=" << heat_added/1000.0 << " kJ/kg" << endl;
    cout << "# geothermal-efficiency=" << thermal_eff << endl;
    cout << "# cycle-thermal-efficiency=" << cycle_eff << endl;
    cout << "# work-done-by-pump=" << work_by_pump/1000.0 << " kJ/kg" << endl;
    cout << "# work-done-down-well=" << (st[2].h-st[1].h)/1000.0 << " kJ/kg" << endl;
    cout << "# work-done-up-well=" << (st[4].h-st[3].h)/1000.0 << " kJ/kg" << endl;
    cout << "# heat-rejected=" << heat_rejected/1000.0 << " kJ/kg" << endl;

    cout << "# Done." << endl;
    return 0;
}
