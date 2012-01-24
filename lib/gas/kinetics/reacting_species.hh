/*  \file reacting_species.hh
 *  \brief Declarations for a reacting species.
 *
 *  \author Rowan J Gollan
 *  \version 24-Feb-2006
 */

#ifndef R_SPECIES_HH
#define R_SPECIES_HH

#include <vector>
#include <valarray>

#include "../../nm/source/no_fuss_linear_algebra.hh"
#include "reaction_pieces.hh"


class ReactingSpecies {
public:
    ReactingSpecies( const std::vector<ReactionPieces*> r_pieces );
    ReactingSpecies( const ReactingSpecies &r );
    virtual ~ReactingSpecies();
    ReactingSpecies* clone();
    
    int nreac() {return int(_r_pieces.size()); }
    double rate_for_reac( int ir, const std::valarray<double> &w_f, const std::valarray<double> &w_b );

    virtual double production( const std::valarray<double> &w_f, const std::valarray<double> &w_b );
    virtual double loss( const std::valarray<double> &w_f, const std::valarray<double> &w_b );

protected:
    std::vector<ReactionPieces*> _r_pieces;
};

#endif
