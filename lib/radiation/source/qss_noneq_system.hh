/** \file coupled_noneq_system.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 24-Mar-09: Second version of QSS nonequilibrium model
 *  \brief Declarations for nonequilibirum population model
 *
 **/

#ifndef QSS_NONEQ_SYSTEM_HH
#define QSS_NONEQ_SYSTEM_HH

#include <vector>
#include <string>
#include <valarray>

#include "../../../lib/gas_models2/source/gas.hh"
#include "../../util/source/config_parser.hh"
#include "cr_reactions.hh"
#include "radiator.hh"
#include "noneq_radiators.hh"

class QSSNoneqSystem {
public:
    QSSNoneqSystem( const string rad_name, const string input_file );
    ~QSSNoneqSystem();
    
public:
    int complete_initialisation( ConfigParser *cfg );
    int solve_system( Gas_data &Q );
    
protected:
    NoneqRadiator * ne_rad_;
    std::vector<CR_Reaction*> reactions_;
    Valmatrix dGdy_;
    std::valarray<double> C_;
    std::valarray<double> y_out_;
};

#endif

