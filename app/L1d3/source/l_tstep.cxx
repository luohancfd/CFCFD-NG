/** \file l_tstep.cxx
 * \ingroup l1d3
 * \brief Time-stepping routines for l1d.cxx.
 *
 * \author PA Jacobs
 */

/*-----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/kinetics/reaction-update.hh"
#include "l_kernel.hh"
#include "l1d.hh"
#include "l_tstep.hh"
#include "../../../lib/nm/source/qd_power.h"
#include "../../../lib/nm/source/qd_log10.h"


LCell::LCell(Gas_model* gmodel)
{
    gas = new Gas_data(gmodel);
    ref = new Gas_data(gmodel);
}


LCell::LCell(const LCell& c)
{
    x = c.x;
    area = c.area;
    pface = c.pface;
    uface = c.uface;
    qstar = c.qstar;
    volume = c.volume;
    xmid = c.xmid;
    T_Wall = c.T_Wall;
    K_over_L = c.K_over_L;
    gas = new Gas_data(*(c.gas));
    ref = new Gas_data(*(c.ref));
    u = c.u;
    shear_stress = c.shear_stress;
    heat_flux = c.heat_flux;
    entropy = c.entropy;
    mass =c.mass;
    moment = c.moment;
    Energy = c.Energy;
    L_bar = c.L_bar;
    x_old = c.x_old;
    mass_old = c.mass_old;
    moment_old = c.moment_old;
    Energy_old = c.Energy_old;
    L_bar_old = c.L_bar_old;
    for ( int i = 0; i < NL; ++i ) {
	DxDt[i] = c.DxDt[i];
	DmDt[i] = c.DmDt[i];
	DmomDt[i] = c.DmomDt[i];
	DEDt[i] = c.DEDt[i];
	DLDt[i] = c.DLDt[i];
    }
    Q_m = c.Q_m;
    Q_mom = c.Q_mom;
    Q_E = c.Q_E;
    dt_chem = c.dt_chem;
    dt_therm = c.dt_therm;
}


LCell::~LCell()
{
    delete ref;
    delete gas;
}


int LCell::encode_conserved()
// Compute the conserved quantities from the primary variables.
{
    // Mass.  Assuming that the cell volume is correct.
    mass = gas->rho * volume;
    // X-momentum.
    moment = mass * u;
    // Total Energy = mass * (specific internal energy + kinetic energy/unit mass).
    Energy = mass * (gas->e[0] + 0.5*u*u);
    return SUCCESS;
} // end encode_conserved()


int LCell::decode_conserved()
// Compute the primary variables from the conserved quantities.
{
    Gas_model *gmodel = get_gas_model_ptr();
    gas->rho = mass / volume;
    u = moment / mass;
    double ke = 0.5*u*u;
    gas->e[0] = Energy/mass - ke;
    if ( gas->rho <= 0.0 || mass <= 0.0 || fabs(u) > 1.0e5 ) {
	printf("LCell_decode_conserved: Bad value for density, mass or velocity\n");
	printf("    rho=%g, mass=%g, u=%g in cell\n", gas->rho, mass, u);
	printf("    xmid=%g, x.right=%g\n", xmid, x );
	printf("    mass=%g, velocity=%g, volume=%g\n", mass, u, volume );
	printf("    momemtum=%g, total Energy=%g e[0]=%g\n", moment, Energy, gas->e[0] );
	gas->print_values();
	exit(BAD_CELLS_ERROR);
    }
    // Fill out the other thermo variables
    gmodel->eval_thermo_state_rhoe(*(gas));
    gmodel->eval_transport_coefficients(*(gas));
    // c->entropy = gmodel->mixture_entropy(*(c->gas)); 
    // FIX-ME -- would like the generic call to entropy to work
    // Entropy referenced to 1 atm and 300K */
    double gam = gmodel->gamma(*(gas));
    entropy = gmodel->Cv(*(gas)) * log(pow(gas->T[0]/300.0, gam) *
				       pow(gas->p/101.3e3, (1.0-gam)));
    return SUCCESS;
} // end decode_conserved()


/// \brief Copy all of the gas-dynamic data from the source cell to the destination cell.
int L_copy_cell_data(LCell *source, LCell *destination, int copy_extras)
{
    // Basic copy: just enough for the application of boundary
    // conditions to the end of the gas slugs.
    destination->gas->copy_values_from(*(source->gas));
    destination->xmid = source->xmid;
    destination->L_bar = source->L_bar;
    destination->u = source->u;
    if (copy_extras == 1) {
	// Now, for the extra items that will be needed for cell refinement.
        destination->x = source->x;
        destination->area = source->area;
        destination->pface = source->pface;
        destination->uface = source->uface;
        destination->volume = source->volume;

        destination->T_Wall = source->T_Wall;
        destination->K_over_L = source->K_over_L;
        destination->shear_stress = source->shear_stress;
        destination->heat_flux = source->heat_flux;
        destination->entropy = source->entropy;

        destination->mass = source->mass;
        destination->moment = source->moment;
        destination->Energy = source->Energy;
    }
    return 0;
} // end function L_copy_cell_data()


int L_blend_cells(LCell* cA, LCell* cB, LCell* c, double alpha, int blend_type)
{
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();

    if ( blend_type == BLEND_PUT ) {
	c->u     = (1.0 - alpha) * cA->u     + alpha * cB->u;
	c->gas->p = (1.0 - alpha) * cA->gas->p + alpha * cB->gas->p;
	c->gas->T[0] = (1.0 - alpha) * cA->gas->T[0] + alpha * cB->gas->T[0];
	for ( int isp = 0; isp <= nsp; ++isp ) {
	    c->gas->massf[isp] = (1.0 - alpha) * cA->gas->massf[isp] + 
		            alpha * cB->gas->massf[isp];
	}
	gmodel->eval_thermo_state_pT(*(c->gas));
	gmodel->eval_transport_coefficients(*(c->gas));
    } else {
	c->u       = (1.0 - alpha) * cA->u       + alpha * cB->u;
	c->gas->rho = (1.0 - alpha) * cA->gas->rho + alpha * cB->gas->rho;
	c->gas->e[0]= (1.0 - alpha) * cA->gas->e[0]+ alpha * cB->gas->e[0];
	c->gas->T[0]= (1.0 - alpha) * cA->gas->T[0]+ alpha * cB->gas->T[0];
	for ( int isp = 0; isp <= nsp; ++isp ) {
	    c->gas->massf[isp] = (1.0 - alpha) * cA->gas->massf[isp] + 
		            alpha * cB->gas->massf[isp];
	}
	gmodel->eval_thermo_state_rhoe(*(c->gas));
	gmodel->eval_transport_coefficients(*(c->gas));
    }
    return 0;
} // end L_blend_cells()
