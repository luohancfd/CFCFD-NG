/// \file block_filter.cxx
/// \brief Filter functions that work across block data.
/// \ingroup eilmer3
///
/// \version 23-Mar-2013 moved out of block.cxx
///

#include <iostream>
#include <unistd.h>
#include <math.h>
#include "cell.hh"
#include "kernel.hh"
#include "block.hh"


/// \brief Part of the spatial filter function.
/// 
/// \param PROP : name of the variable to be filtered.
#define SIMPLE_FILTER(PROP)						\
    for (i = imin; i <= imax; ++i) {					\
        for (j = jmin; j <= jmax; ++j) {				\
            cell = get_cell(i,j);					\
	    cN = get_cell(i,j+1);					\
	    cE = get_cell(i+1,j);					\
	    cS = get_cell(i,j-1);					\
	    cW = get_cell(i-1,j);					\
            (cell)->fs->PROP = (1.0 - alpha) * (cell)->fs->PROP +	\
		alpha * 0.25 * ((cN)->fs->PROP + (cE)->fs->PROP + (cS)->fs->PROP + (cW)->fs->PROP); \
        }								\
    }


/// \brief Flux Corrected Transport filter. 
/// Applies diffusion and anti-diffusion step
/// as a filter.
/// Part of the spatial filter function. Based on 
/// technical report by M. N. Macrossan, in turn based
/// on presentation by Kaplan and Oran.
/// 
/// \param PROP : name of the variable to be filtered.
#define FCT_DIFFUSION(PROP)				\
    CALC_DIFFUSION(PROP)				\
    APPLY_DIFFUSION(PROP)				\

#define FCT_ANTI_DIFFUSION(PROP)			\
    CALC_ANTI_DIFFUSION(PROP)				\
    APPLY_ANTI_DIFFUSION(PROP)				\

#define CALC_DIFFUSION(PROP)				\
    for (i = imin; i <= imax; ++i) {			\
        for (j = jmin; j <= jmax; ++j) {		\
            cell = get_cell(i,j);			\
	    cN = get_cell(i,j+1);			\
	    cE = get_cell(i+1,j);			\
	    cS = get_cell(i,j-1);			\
	    cW = get_cell(i-1,j);			\
            diffuse[i][j] = ( 1 - 4 * mu ) * (cell)->fs->PROP + \
		mu * ( (cN)->fs->PROP + (cE)->fs->PROP + (cS)->fs->PROP + (cW)->fs->PROP ); \
	}						\
    }

#define APPLY_DIFFUSION(PROP)				\
    for (i = imin; i <= imax; ++i) {			\
        for (j = jmin; j <= jmax; ++j) {		\
            cell = get_cell(i,j);			\
            (cell)->fs->PROP = diffuse[i][j];		\
	}						\
    }

#define CALC_ANTI_DIFFUSION(PROP)			\
    for (i = imin; i <= imax+1; ++i) {			\
        for (j = jmin; j <= jmax; ++j) {		\
	    p2 = get_cell(i+1,j)->fs->PROP;		\
	    p1 = get_cell(i,j)->fs->PROP;		\
	    m1 = get_cell(i-1,j)->fs->PROP;		\
	    m2 = get_cell(i-2,j)->fs->PROP;		\
            iadflux[i][j] = calc_anti_diffusive_flux(m2, m1,		\
						     p1, p2, mu);	\
	}						\
    }							\
    for (i = imin; i <= imax; ++i) {			\
        for (j = jmin; j <= jmax+1; ++j) {		\
	    p2 = get_cell(i,j+1)->fs->PROP;		\
	    p1 = get_cell(i,j)->fs->PROP;		\
	    m1 = get_cell(i,j-1)->fs->PROP;		\
	    m2 = get_cell(i,j-2)->fs->PROP;		\
            jadflux[i][j] = calc_anti_diffusive_flux(m2, m1,		\
						     p1, p2, mu);	\
	}								\
    }

#define APPLY_ANTI_DIFFUSION(PROP)			\
    for (i = imin; i <= imax; ++i) {			\
        for (j = jmin; j <= jmax; ++j) {		\
            cell = get_cell(i,j);			\
	    (cell)->fs->PROP += iadflux[i][j] + jadflux[i][j] - \
		iadflux[i+1][j] - jadflux[i][j+1];		\
	}							\
    }

double Block::calc_anti_diffusive_flux(double m2, double m1, double p1, double p2, double mu)
{
    double delm = m1 - m2;
    double del =  p1 - m1;
    double delp = p2 - p1;
    double S = copysign(1.0, del);
    double delmin = min(S * delm, S * delp);

    return S * max(0.0, min(mu * fabs(del), delmin));
}

int Block::apply_spatial_filter_diffusion(double mu, size_t npass, size_t dimensions)
/// \brief Filter the cell-centred primary variables.
///
/// This filtering is done on a block-by-block basis.
/// Valid flow data are needed in the ghost cells at the edges.
/// \param alpha : filter coefficient (closer to 1.0, more fudging)
/// \param npass : this many passes of the simple averager
//
// To do: We should fix for 3D or remove. I think that there are no cases that need it. 
{
    global_data &G = *get_global_data_ptr();
    FV_Cell *cell, *cN, *cE, *cS, *cW;
    size_t isp;
    size_t i, j;
    Gas_model *gm = get_gas_model_ptr();
    size_t nsp = gm->get_number_of_species();
    vector<vector<double> > diffuse;
    diffuse.resize(nni+4);
    for (size_t i = 0; i < diffuse.size(); ++i ) {
	diffuse[i].resize(nnj+4);
    }

    // Apply the "standard filter".
    //for ( ipass = 0; ipass < npass; ++ipass ) {
	//SIMPLE_FILTER(gas->rho)
	//SIMPLE_FILTER(gas->e[0])
	//SIMPLE_FILTER(vel.x)
        //SIMPLE_FILTER(vel.y)
        //for (isp = 0; isp < nsp; ++isp) {
        //    SIMPLE_FILTER(gas->massf[isp])
	FCT_DIFFUSION(gas->rho)
	FCT_DIFFUSION(gas->e[0])
	FCT_DIFFUSION(vel.x)
        FCT_DIFFUSION(vel.y)
        for (isp = 0; isp < nsp; ++isp) {
            FCT_DIFFUSION(gas->massf[isp])
        }
	//}

    // We should make the thermodynamic state consistent, at least.
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    cell = get_cell(i,j);
	    Gas_data *gas= cell->fs->gas;
	    gm->eval_thermo_state_rhoe(*gas);
	    if ( G.viscous ) gm->eval_transport_coefficients(*gas, gm);
	    if ( G.diffusion ) gm->eval_diffusion_coefficients(*gas);
	} // j loop
    } // i loop
    return SUCCESS;
} // end of apply_spatial_filter_diffusion()

int Block::apply_spatial_filter_anti_diffusion(double mu, size_t npass, size_t dimensions)
/// \brief Filter the cell-centred primary variables.
///
/// This filtering is done on a block-by-block basis.
/// Valid flow data are needed in the ghost cells at the edges.
/// \param alpha : filter coefficient
/// \param npass : this many passes
//
// To do: We should fix for 3D or remove. I think that there are no cases that need it. 
{
    global_data &G = *get_global_data_ptr();
    FV_Cell *cell;
    double m2, m1, p1, p2;
    size_t isp;
    size_t i, j;
    Gas_model *gm = get_gas_model_ptr();
    size_t nsp = gm->get_number_of_species();
    vector<vector<double> > iadflux;
    vector<vector<double> > jadflux;
    iadflux.resize(nni+5);
    for (size_t i = 0; i < iadflux.size(); ++i ) {
	iadflux[i].resize(nnj+4);
    }
    jadflux.resize(nni+4);
    for (size_t i = 0; i < jadflux.size(); ++i ) {
	jadflux[i].resize(nnj+5);
    }
    // Apply the anti-diffusion.
    FCT_ANTI_DIFFUSION(gas->rho)
    FCT_ANTI_DIFFUSION(gas->e[0])
    FCT_ANTI_DIFFUSION(vel.x)
    FCT_ANTI_DIFFUSION(vel.y)
    for (isp = 0; isp < nsp; ++isp) {
	FCT_ANTI_DIFFUSION(gas->massf[isp])
    }
    
    // We should make the thermodynamic state consistent, at least.
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    cell = get_cell(i,j);
	    Gas_data *gas= cell->fs->gas;
	    gm->eval_thermo_state_rhoe(*gas);
	    if ( G.viscous ) gm->eval_transport_coefficients(*gas, gm);
	    if ( G.diffusion ) gm->eval_diffusion_coefficients(*gas);
	} // j loop
    } // i loop
    return SUCCESS;
} // end of apply_spatial_filter_anti_diffusion()
