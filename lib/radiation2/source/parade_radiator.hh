/** \file parade_radiator.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 11-Jan-12: Initial implementation
 *
 *  \brief Declarations for class describing a parade radiator
 *
 **/

#ifndef PARADE_RADIATOR_HH
#define PARADE_RADIATOR_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

class ParadeRadiator {
public:
    /// \brief Constructor
    ParadeRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    virtual ~ParadeRadiator();

public:
    std::string name;
    std::string type;
    int isp;
    int iTr;
    int iTv;
    int nsys;
    std::vector<std::string> sys_names;

    double m_w;
};

ParadeRadiator *create_new_parade_radiator( lua_State * L, const std::string name );

#endif
