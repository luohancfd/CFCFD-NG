/** \file visc.hh
 * \ingroup eilmer3
 * \brief Prototypes for the viscous-flux functions.
 */

#ifndef VISC_HH
#define VISC_HH

int estimate_turbulence_viscosity(struct global_data *gdp, Block *bdp);
int viscous_flux_2D(Block *bdp);
int viscous_derivatives_2D(Block *bdp, size_t gtl);
int viscous_derivatives_edges(Block *bdp, size_t gtl);
int viscous_derivatives_corners(Block *bdp, size_t gtl);

#endif
