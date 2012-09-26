// l_kernel.hh

#ifndef L_KERNEL_HH
#define L_KERNEL_HH

#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/kinetics/reaction-update.hh"
#include "../../../lib/gas/kinetics/energy-exchange-update.hh"

// Simulation-control parameters
class SimulationData {
public:
    int test_case;
    int nslug;          /* number of gas slugs        */
    int npiston;        /* number of pistons          */
    int ndiaphragm;     /* number of diaphragms       */

    int max_step;       /* global iteration limit     */
    double sim_time;    /* present simulation time    */
    double max_time;    /* final solution time, s     */
    double dt_init;     /* initial time step          */
    double dt_global;   /* simulation time step       */
    double dt_allow;    /* allowable global time step */
    double CFL;         /* target CFL (worst case)    */
    int Torder;         /* order of time-stepping     */
    int Xorder;         /* order of reconstruction    */
    int fr_chem;        /* flag to activate finite-rate chemistry */
    double k;           /* "thermal conductivity" for */
                        /* damping temperature glitch */

    int n_dt_plot;      /* number changes in plot intervals */
    std::vector<double> t_change; /* times at which dt_plot changes */
    std::vector<double> dt_plot; /* interval for writing soln  */
    std::vector<double> dt_his;  /* interval for writing sample */
    int hnloc;          /* number of history locations */
    std::vector<double> hxloc;   /* history locations          */
    int hncell;         /* history cell count (all slugs) */

    SimulationData(std::string config_file_name, int echo_input);
    ~SimulationData();
}; // end SimulationData


// Managed gas models
Gas_model *set_gas_model_ptr(Gas_model *gmptr);
Gas_model *get_gas_model_ptr();
int set_reaction_update(std::string file_name);
Reaction_update *get_reaction_update_ptr();
int set_energy_exchange_update( std::string file_name );
Energy_exchange_update *get_energy_exchange_update_ptr();
void L_set_case_id(int id);
int L_get_case_id(void);

#endif
