/** \file equilibrium_air.hh
 *  \ingroup radiation
 *
 *  \author Rowan J. Gollan
 *  \version 13-Jan-07 : ported from work of Sep-Oct 2003
 *           06-Jul-09 : ported from old lib/radiation
 *
 **/

#ifndef EQUILIBRIUM_AIR_HH
#define EQUILIBRIUM_AIR_HH

#include <vector>

#include "../../util/source/useful.h"
#include "spectral_model.hh"

class EquilibriumAir : public RadiationSpectralModel {
public:
    EquilibriumAir();

    ~EquilibriumAir();
    
    std::string str() const;

private:

    double integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved );
    
    double variably_integrated_emission_for_gas_state( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u, bool spectrally_resolved )
    { return 0.0; }
    
    double Planck_absorption(Gas_data &Q);
    
    void spectra_for_gas_state( Gas_data &Q, CoeffSpectra &X );
    
    void spectral_distribution_for_gas_state(Gas_data &Q, std::vector<double> &nus);
    
    void write_line_widths( Gas_data &Q );
};

#endif
