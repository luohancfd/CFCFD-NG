// l_cell.cxx
// Contains Lagrangian-cell related code.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <numeric>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/kinetics/reaction-update.hh"
#include "l_kernel.hh"
#include "l1d.hh"
#include "l_cell.hh"
#include "../../../lib/nm/source/qd_power.h"
#include "../../../lib/nm/source/qd_log10.h"

LFlowState::LFlowState(Gas_model* gmodel)
{
    gas = new Gas_data(gmodel);
    u = 0.0;
}


LFlowState::LFlowState(const LFlowState& fs)
{
    gas = new Gas_data(*(fs.gas));
    u = fs.u;
}


LFlowState::~LFlowState()
{
    delete gas;
}


LCell::LCell(Gas_model* gmodel)
{
    // Let most elements default to zero values.
    x = 0.0;
    area = 0.0;
    pface = 0.0;
    uface = 0.0;
    qstar = 0.0;
    volume = 0.0;
    xmid = 0.0;
    T_Wall = 0.0;
    K_over_L = 0.0;
    gas = new Gas_data(gmodel);
    ref = new Gas_data(gmodel);
    u = 0.0;
    shear_stress = 0.0;
    heat_flux = 0.0;
    entropy = 0.0;
    mass = 0.0;
    moment = 0.0;
    Energy = 0.0;
    L_bar = 0.0;
    x_old = 0.0;
    mass_old = 0.0;
    moment_old = 0.0;
    Energy_old = 0.0;
    L_bar_old = 0.0;
    for ( int i = 0; i < NL; ++i ) {
	DxDt[i] = 0.0;
	DmDt[i] = 0.0;
	DmomDt[i] = 0.0;
	DEDt[i] = 0.0;
	DLDt[i] = 0.0;
    }
    Q_m = 0.0;
    Q_mom = 0.0;
    Q_E = 0.0;
    dt_chem = 0.0;
    dt_therm = 0.0;
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
    mass = c.mass;
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
    double e = accumulate(gas->e.begin(), gas->e.end(), 0.0);
    Energy = mass * (e + 0.5*u*u);
    return SUCCESS;
} // end encode_conserved()


int LCell::decode_conserved()
// Compute the primary variables from the conserved quantities.
{
    Gas_model *gmodel = get_gas_model_ptr();
    gas->rho = mass / volume;
    u = moment / mass;
    double ke = 0.5*u*u;
    double e_non_translation = accumulate(gas->e.begin()+1, gas->e.end(), 0.0);
    // Translational energy is them just what remains after we remove the others.
    gas->e[0] = Energy/mass - ke - e_non_translation;
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
int LCell::copy_data_from(LCell& source, int copy_extras)
{
    // Basic copy: just enough for the application of boundary
    // conditions to the end of the gas slugs.
    gas->copy_values_from(*(source.gas));
    xmid = source.xmid;
    L_bar = source.L_bar;
    u = source.u;
    if (copy_extras == 1) {
	// Now, for the extra items that will be needed for cell refinement.
        x = source.x;
        area = source.area;
        pface = source.pface;
        uface = source.uface;
        volume = source.volume;

        T_Wall = source.T_Wall;
        K_over_L = source.K_over_L;
        shear_stress = source.shear_stress;
        heat_flux = source.heat_flux;
        entropy = source.entropy;

        mass = source.mass;
        moment = source.moment;
        Energy = source.Energy;
    }
    return SUCCESS;
} // end copy_data()


std::string LCell::write_iface_values_to_string()
// Write the flow solution for a cell to a string.
{
    // The new format for L1d3 puts everything onto one line.
    std::ostringstream ost;
    ost.setf(std::ios_base::scientific);
    ost.precision(12);
    ost << x << " " << area; 
    // Don't put the newline char on the end.
    return ost.str();
} // end of write_iface_values_to_string()


int LCell::scan_iface_values_from_string(char *bufptr)
// Scan a string, extracting the data for an interface between cells.
// There isn't any checking of the file content.
// If anything gets out of place, the result is wrong data.
{
    // Look for a new-line character and truncate the string there.
    char *cptr = strchr(bufptr, '\n');
    if ( cptr != NULL ) cptr = '\0';
    // Now, we should have a string with only numbers separated by spaces.
    x = atof(strtok( bufptr, " " )); // tokenize on space characters
    area = atof(strtok( NULL, " " ));
    return SUCCESS;
} // end scan_iface_values_from_string()


std::string LCell::write_cell_values_to_string()
// Write the flow solution for a cell to a string.
{
    // The new format for L1d3 puts everything onto one line.
    std::ostringstream ost;
    ost.setf(std::ios_base::scientific);
    ost.precision(12);
    ost << xmid << " " 
	<< volume << " " 
	<< u << " " 
	<< L_bar << " " 
	<< gas->rho << " " 
	<< gas->p << " " 
	<< gas->a << " " 
	<< shear_stress << " " 
	<< heat_flux << " " 
	<< entropy;
    // Species mass fractions.
    size_t nsp = gas->massf.size();
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	ost << " " << gas->massf[isp];
    }
    if ( nsp > 1 ) ost << " " << dt_chem;
    // Individual energies (in e, T pairs)
    size_t nmodes = gas->T.size();
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	ost << " " << gas->e[imode] << " " << gas->T[imode];
    }
    if ( nmodes > 1 ) ost << " " << dt_therm;
    // Don't put the newline char on the end.
    return ost.str();
} // end of write_cell_values_to_string()


int LCell::scan_cell_values_from_string(char *bufptr)
// Scan a string, extracting the data for a 
// There isn't any checking of the file content.
// If anything gets out of place, the result is wrong data.
{
    // Look for a new-line character and truncate the string there.
    char *cptr = strchr(bufptr, '\n');
    if ( cptr != NULL ) cptr = '\0';
    // Now, we should have a string with only numbers separated by spaces.
    xmid = atof(strtok( bufptr, " " )); // tokenize on space characters
    volume = atof(strtok( NULL, " " ));
    u = atof(strtok( NULL, " " ));
    L_bar = atof(strtok( NULL, " " ));
    gas->rho = atof(strtok( NULL, " " ));
    gas->p = atof(strtok( NULL, " " ));
    gas->a = atof(strtok( NULL, " " ));
    shear_stress = atof(strtok( NULL, " " ));
    heat_flux = atof(strtok( NULL, " " ));
    entropy = atof(strtok( NULL, " " ));
    size_t nsp = gas->massf.size();
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	gas->massf[isp] = atof(strtok( NULL, " " ));
    }
    if ( nsp > 1 ) dt_chem = atof(strtok( NULL, " " ));
    size_t nmodes = gas->T.size();
    for ( size_t imode = 0; imode < nmodes; ++imode ) {
	gas->e[imode] = atof(strtok( NULL, " " ));
	gas->T[imode] = atof(strtok( NULL, " " ));
    }
    if ( nmodes > 1 ) dt_therm = atof(strtok( NULL, " " ));
    return SUCCESS;
} // end scan_cell_values_from_string()

//----------------------------------------------------------------------

int L_blend_cells(LCell& cA, LCell& cB, LCell& c, double alpha, int blend_type)
{
    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t nmodes = gmodel->get_number_of_modes();

    if ( blend_type == BLEND_PUT ) {
	c.u     = (1.0 - alpha) * cA.u     + alpha * cB.u;
	c.gas->p = (1.0 - alpha) * cA.gas->p + alpha * cB.gas->p;
	c.gas->T[0] = (1.0 - alpha) * cA.gas->T[0] + alpha * cB.gas->T[0];
	for ( size_t isp = 0; isp <= nsp; ++isp ) {
	    c.gas->massf[isp] = (1.0 - alpha) * cA.gas->massf[isp] + 
		            alpha * cB.gas->massf[isp];
	}
	gmodel->eval_thermo_state_pT(*(c.gas));
	gmodel->eval_transport_coefficients(*(c.gas));
    } else {
	c.u       = (1.0 - alpha) * cA.u       + alpha * cB.u;
	c.gas->rho = (1.0 - alpha) * cA.gas->rho + alpha * cB.gas->rho;
	for ( size_t isp = 0; isp <= nsp; ++isp ) {
	    c.gas->massf[isp] = (1.0 - alpha) * cA.gas->massf[isp] + alpha * cB.gas->massf[isp];
	}
	for ( size_t imode = 0; imode < nmodes; ++imode ) {
	    c.gas->e[imode]= (1.0 - alpha) * cA.gas->e[imode]+ alpha * cB.gas->e[imode];
	    c.gas->T[imode]= (1.0 - alpha) * cA.gas->T[imode]+ alpha * cB.gas->T[imode];
	}
	gmodel->eval_thermo_state_rhoe(*(c.gas));
	gmodel->eval_transport_coefficients(*(c.gas));
    }
    return SUCCESS;
} // end L_blend_cells()
