/** \file spradian_radiator.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 11-Jan-12: Initial implementation
 *
 *  \brief Declarations for class describing a spradian radiator
 *
 **/

#ifndef SPRADIAN_RADIATOR_HH
#define SPRADIAN_RADIATOR_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

class SpradianRadiator {
public:
    /// \brief Constructor
    SpradianRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    virtual ~SpradianRadiator();

public:
    void set_conc_pointer( double * concp )
    { conc = concp; }

    void set_concentration( double rho_i );

    double get_concentration()
    { return *conc; }

public:
    std::string name;
    std::string type;
    int isp;
    int nsys;
    std::vector<std::string> sys_names;

    double m_w;

private:
    double * conc;
};

SpradianRadiator *create_new_spradian_radiator( lua_State * L, const std::string name );

#endif
