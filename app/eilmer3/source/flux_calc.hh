/** \file flux_calc.hh
 * \ingroup eilmer3
 * \brief Flux calculators for Euler terms.
 *
 * \author PJ
 * \version 29-Jun-04, initially extracted from mb_cns.h
 * \version July-2008, Elmer3 puts cell code, flux calculators and flow code back together.
 *
 */

#ifndef FLUX_CALC_HH
#define FLUX_CALC_HH
/* Go ahead and make some definitions... */

#include "../../../lib/geometry2/source/geom.hh"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "c-flow-condition.hh"
#include "cell.hh"

const int DEBUG_FLUX = 0;

/** \brief Type of flux calculation...
 *
 * \verbatim
 * FLUX_RIEMANN  : 3-stage approximate Riemann solver.
 * FLUX_AUSM     : AUSM flux-splitting
 * FLUX_EFM      : Mike Macrossan's EFM flux calculation
 * FLUX_AUSMDV   : Wada and Liou's flux calculator AIAA Paper 94-0083
 * FLUX_ADAPTIVE : EFM near shocks, AUSMDV otherwise
 * FLUX_AUSM_PLUS_UP : Liou's 2006 all-speed flux calculator
 * FLUX_HLLE     : MHD HLLE approximate Riemann solver
 * \endverbatim
 */
enum flux_calc_t {FLUX_RIEMANN, FLUX_AUSM, FLUX_EFM, FLUX_AUSMDV,
		  FLUX_ADAPTIVE, FLUX_AUSM_PLUS_UP, FLUX_HLLE};
flux_calc_t set_flux_calculator(std::string name);
flux_calc_t get_flux_calculator();
std::string get_flux_calculator_name(flux_calc_t calc); 

int compute_interface_flux(FlowState &Lft, FlowState &Rght, FV_Interface &IFace, double omegaz=0.0);
int set_flux_vector_in_local_frame(ConservedQuantities &F, const FlowState &fs);
int set_flux_vector_in_global_frame(FV_Interface &IFace, FlowState &fs, double omegaz);

/* rivp.c */
int rivp(FlowState &QL, FlowState &QR, FlowState &QIF, double &WSL, double &WSR);
int rivp_stage_3(FlowState &QL, FlowState &QR, FlowState &QLstar, FlowState &QRstar,
                 double WSL, double WSR, double geff, FlowState &QIF);
/* ausm.c */
int ausm(FlowState &QL, FlowState &QR, FlowState &QIF, double &WSL, double &WSR);
/* efm.c */
int efmflx(FlowState &Lft, FlowState &Rght, FV_Interface &IFace);
int exxef(double sn, double &exx, double &ef);
/* ausmdv.c */
int ausmdv(FlowState &Lft, FlowState &Rght, FV_Interface &IFace);
/* adaptive_flux.c */
int adaptive_flux(FlowState &Lft, FlowState &Rght, FV_Interface &IFace);
/* ausm_plus_up.c */
int ausm_plus_up(FlowState &Lft, FlowState &Rght, FV_Interface &IFace);
/* hlle.cxx */
int hlle(FlowState &Lft, FlowState &Rght, FV_Interface &IFace);

#endif
