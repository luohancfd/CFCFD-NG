#include <iostream>
#include <vector>

#include "../../nm/source/no_fuss_linear_algebra.hh"
#include "reaction_pieces.hh"

using namespace std;

ReactionPieces::ReactionPieces( int ir, int nu )
    : ir_( ir ), nu_( nu ) {}

ReactionPieces::ReactionPieces( const ReactionPieces &r )
    : ir_( r.ir_ ), nu_( r.nu_ ) {}

ReactionPieces::~ReactionPieces() {}

ReactionPiecesA::ReactionPiecesA( int ir, int nu )
    : ReactionPieces( ir, nu ) {}

ReactionPiecesA::ReactionPiecesA( const ReactionPiecesA &r )
    : ReactionPieces( r.ir_, r.nu_ ) {}

ReactionPiecesA::~ReactionPiecesA() {}

ReactionPiecesA*
ReactionPiecesA::
clone()
{
    return new ReactionPiecesA(*this);
}

double
ReactionPiecesA::
production( const vector<double> &w_f, const vector<double> &w_b )
{
    return nu_*w_f[ir_];
}

double
ReactionPiecesA::
loss( const vector<double> &w_f, const vector<double> &w_b )
{
    return nu_*w_b[ir_];
}

ReactionPiecesB::ReactionPiecesB( int ir, int nu )
    : ReactionPieces( ir, nu ) {}

ReactionPiecesB::ReactionPiecesB( const ReactionPiecesB &r )
    : ReactionPieces( r.ir_, r.nu_ ) {}

ReactionPiecesB::~ReactionPiecesB() {}

ReactionPiecesB*
ReactionPiecesB::
clone()
{
    return new ReactionPiecesB(*this);
}


double
ReactionPiecesB::
production( const vector<double> &w_f, const vector<double> &w_b )
{
     return -nu_*w_b[ir_];
}

double
ReactionPiecesB::
loss( const vector<double> &w_f, const vector<double> &w_b )
{
    return -nu_*w_f[ir_];
}


