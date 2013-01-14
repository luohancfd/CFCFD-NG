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

#define DEBUG_FLUX 0

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
#define FLUX_RIEMANN       0
#define FLUX_AUSM          1
#define FLUX_EFM           2
#define FLUX_AUSMDV        3
#define FLUX_ADAPTIVE      4
#define FLUX_AUSM_PLUS_UP  5
#define FLUX_HLLE          6

/*
 * Function prototypes...
 */
int set_flux_calculator(int iflux);
int get_flux_calculator(void);
int compute_interface_flux(FlowState &Lft, FlowState &Rght, FV_Interface &IFace, double omegaz=0.0);
int set_interface_flux(FV_Interface &IFace, FlowState *IFace_flow_state);
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
