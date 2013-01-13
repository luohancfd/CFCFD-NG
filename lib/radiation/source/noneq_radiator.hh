/** \file noneq_radiators.hh
 *  \ingroup radiation2
 *
 *  \author Daniel F. Potter
 *  \version 03-Apr-10
 *  \brief Definitions for nonequilibrium radiator class (port from lib/radiation)
 *
 **/
 
#ifndef NONEQ_RADIATOR_HH
#define NONEQ_RADIATOR_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "radiator.hh"
 
class NoneqElecState {
public:
     NoneqElecState( int index, std::string label ) : index_(index), label_(label) {};
    ~NoneqElecState() {};
    
public:
    int ne_index_;	// index in the whole noneq system
    int index_;		// electron level index for radiator
    std::string label_;
};

class NoneqRadiator {
public:
    NoneqRadiator( const std::string name, const std::string section, ConfigParser * cfg, bool radiators_present );
    ~NoneqRadiator();
public:
    size_t get_nelev_from_elev( int ie_index );
public:
    std::vector<NoneqElecState*> elec_states_;
    Radiator * rad_pointer_;
    std::string name_;
    int irad_;
    bool all_levels_noneq_;
    bool boltz_eqs_;
    vector<double> boltz_fractions_;
    double c_;		// moles / cm**3 of this radiator
};

NoneqRadiator * get_noneq_radiator_pointer( int ne_irad );

NoneqRadiator * get_noneq_radiator_pointer_from_name( std::string name );

NoneqRadiator * create_new_noneq_radiator( const string name, const string section, ConfigParser * cfg, bool radiators_present=false );

void clear_noneq_radiators( void );

void decode_noneq_label( std::string label, int &irad, int &ielev, int &ivlev, int &ne_index );

int get_number_of_noneq_radiators( void );

#endif
