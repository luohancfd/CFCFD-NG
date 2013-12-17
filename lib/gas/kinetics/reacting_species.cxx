/*  \file reacting_species.cxx
 *  \brief Deginitions for a reacting species.
 *
 *  \author Rowan J Gollan
 *  \version 24-Feb-2006
 **/

#include <iostream>
#include <vector>
#include "../../nm/source/no_fuss_linear_algebra.hh"
#include "reaction_pieces.hh"
#include "reacting_species.hh"

using namespace std;

ReactingSpecies::ReactingSpecies( const vector<ReactionPieces*> r_pieces )
{
     for( size_t i = 0; i < r_pieces.size(); ++i ) {
 	_r_pieces.push_back( r_pieces[i]->clone() );
     }
}

ReactingSpecies::ReactingSpecies( const ReactingSpecies &r )
{
    for( size_t i = 0; i < r._r_pieces.size(); ++i ) {
	_r_pieces.push_back( r._r_pieces[i]->clone() );
    }
}

ReactingSpecies::~ReactingSpecies()
{
    for( size_t i = 0; i < _r_pieces.size(); ++i ) {
	delete _r_pieces[i];
    }
}

ReactingSpecies* ReactingSpecies::clone()
{
    return new ReactingSpecies(*this);
}

double
ReactingSpecies::
rate_for_reac( int ir, const vector<double> &w_f, const vector<double> &w_b )
{
    return (_r_pieces[ir]->production(w_f, w_b) - _r_pieces[ir]->loss(w_f, w_b));

}

double
ReactingSpecies::
production( const vector<double> &w_f, const vector<double> &w_b )
{
    double val = 0.0;
    for( size_t i = 0; i < _r_pieces.size(); ++i ) {
	val += _r_pieces[i]->production( w_f, w_b );
    }
    return val;
}

double 
ReactingSpecies::
loss( const vector<double> &w_f, const vector<double> &w_b )
{
    double val = 0.0;
    for( size_t i = 0; i < _r_pieces.size(); ++i ) {
	val += _r_pieces[i]->loss( w_f, w_b );
    }
    return val;
}

