// l_cell.hh
// Refactored to contain mostly cell-related code, 26-Sep-2012

#ifndef L_CELL_HH
#define L_CELL_HH

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"

// Gas states at the interfaces.
// This is used for interfacing to the Riemann solver.
struct L_flow_state
{
    Gas_data *gas;
    double u;       /* normal velocity, m/s           */
};


// Data stored in each Lagrangian cell.
class LCell {
public:
    /* GEOMETRY -- INTERFACES */
    double x;           /* Interface position, m          */
                        /* This is the interface to the   */
                        /* right of the cell midpoint.    */
    double area;        /* Interface area, m**2           */
    double pface;       /* Interface pressure (Riemann)   */
    double uface;       /* Interface velocity (Riemann)   */
    double qstar;       /* axial heat flux (fudge)        */
    double volume;      /* Cell volume, m**3              */
    double xmid;        /* location of midpoint           */
    double T_Wall;      /* specified wall temperature     */
    double K_over_L;    /* pipe-fitting loss coeff.       */

    /* CELL-AVERAGE VARIABLES */
    Gas_data *gas;   /* core values of gas properties             */
    Gas_data *ref;   /* reference values for some viscous effects */
    double u;              /* bulk velocity, m/s                      */
    double shear_stress;   /* Wall shear stress, N/m**2                 */
    double heat_flux;      /* Heat flux at wall, W/m**2                 */
    double entropy;        /* Entropy referenced to 1 atm and 300K      */

    /* CONSERVED VARIABLES */
    double mass;        /* cell mass                      */
    double moment;      /* X-momentum/unit volume         */
    double Energy;      /* Total Energy/unit volume       */

    /* OTHER INTEGRATED VARIABLES */
    double L_bar;       /* This length scale is the       */
                        /* Integral of local velocity     */
                        /* for use in Mirel's model of    */
                        /* the tube-wall Boundary Layer   */

    /* A record of the cell state */
    /* (for adaptive stepping) */
    double x_old, mass_old, moment_old, Energy_old, L_bar_old;

    /* TIME DERIVATIVES */
    double DxDt[NL];    /* updates for interface posn.    */
    double DmDt[NL];    /* Mass                           */
    double DmomDt[NL];  /* X-momentum                     */
    double DEDt[NL];    /* Total Energy                   */

    double DLDt[NL];    /* Length scale, L_bar            */

    /* PRODUCTION VECTOR */
    double Q_m;         /* Mass from sources or sinks     */
    double Q_mom;       /* X-Momentum from body forces or */
                        /* wall friction                  */
    double Q_E;         /* Total Energy production or     */
                        /* transfer to the wall           */

    double dt_chem;     // Need to remember the suggested timestep.
    double dt_therm;    // and for nonequilibrium thermodynamics.

    LCell(Gas_model* gmodel);
    LCell(const LCell& c);
    ~LCell();
    int encode_conserved();
    int decode_conserved();
};


#define BLEND_PUT 0
#define BLEND_RHOUE 1

int L_copy_cell_data(LCell *source, LCell *target, int copy_extras);
int L_blend_cells(LCell *cA, LCell *cB, LCell *c, double alpha, int blend_type );

#endif
