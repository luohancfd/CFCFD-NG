/*  \file reaction_pieces.cxx
 *  \brief Declarations for the ReactionPieces class
 *
 *  \author Rowan J Gollan
 *  \version 24-Feb-2006
 **/

#ifndef R_PIECES_HH
#define R_PIECES_HH

#include <vector>

#include "../../nm/source/no_fuss_linear_algebra.hh"


class ReactionPieces {
public:
    ReactionPieces( int ir, int nu );
    ReactionPieces( const ReactionPieces &r );
    virtual ~ReactionPieces();
    virtual ReactionPieces* clone() = 0;

    int ir() { return ir_; };

    virtual double production( const std::vector<double> &w_f,
			       const std::vector<double> &w_b ) = 0;
    virtual double loss( const std::vector<double> &w_f,
			 const std::vector<double> &w_b ) = 0;
protected:
    int ir_;
    int nu_;
};

////////// POSITIVE NU /////////////////

class ReactionPiecesA : public ReactionPieces {
public:
    ReactionPiecesA( int ir, int nu );
    ReactionPiecesA( const ReactionPiecesA &r );
    virtual ~ReactionPiecesA();
    virtual ReactionPiecesA* clone();

    double production( const std::vector<double> &w_f,
		       const std::vector<double> &w_b );
    double loss( const std::vector<double> &w_f,
		 const std::vector<double> &w_b );
};

////////// NEGATIVE NU /////////////////

class ReactionPiecesB : public ReactionPieces {
public:
    ReactionPiecesB( int ir, int nu );
    ReactionPiecesB( const ReactionPiecesB &r );
    virtual ~ReactionPiecesB();
    virtual ReactionPiecesB* clone();

    double production( const std::vector<double> &w_f,
		       const std::vector<double> &w_b );
    double loss( const std::vector<double> &w_f,
		 const std::vector<double> &w_b );
};


#endif



