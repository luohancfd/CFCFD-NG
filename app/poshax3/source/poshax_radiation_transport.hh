/** \file poshax_radiation_transport.hh
 *  \ingroup poshax2
 *  \brief Header file for the poshax radiation transport model class and functions.
 **/

#ifndef POSHAX_RADIATION_TRANSPORT_HH
#define POSHAX_RADIATION_TRANSPORT_HH

#include <string>
#include <vector>

#include "../../../lib/radiation2/source/spectral_model.hh"
#include "../../../lib/util/source/lua_service.hh"
#include "../../../lib/gas/models/gas_data.hh"

class PoshaxRadiationTransportModel {
public:
    PoshaxRadiationTransportModel( void ) {}
    
    PoshaxRadiationTransportModel( lua_State *L );
    
    virtual ~PoshaxRadiationTransportModel();
    
    int set_radiation_spectral_model( std::string file_name );
    
    RadiationSpectralModel * get_rsm_pointer()
    { return rsm_; }

    virtual std::string str() const = 0;
    
    virtual double eval_Q_rad( Gas_data &Q ) = 0;
    
protected:
    int spectrally_resolved_;
    
    RadiationSpectralModel* rsm_;
};

class OpticallyThin : public PoshaxRadiationTransportModel {
public:
    OpticallyThin( lua_State *L );
    
    ~OpticallyThin();

    std::string str() const;
    
    int initialise();
    
    double eval_Q_rad( Gas_data &Q );
};

class OpticallyThick : public PoshaxRadiationTransportModel {
public:
    OpticallyThick( lua_State *L );
    
    ~OpticallyThick();

    std::string str() const;

    double eval_Q_rad( Gas_data &Q );
};


class OpticallyVariable : public PoshaxRadiationTransportModel {
public:
    OpticallyVariable( lua_State *L );
    
    ~OpticallyVariable();

    std::string str() const;

    double eval_Q_rad( Gas_data &Q );
    
private:
    double wavel_switch_;
    double Lambda_lower_;
    double Lambda_upper_;
};

PoshaxRadiationTransportModel * create_poshax_radiation_transport_model( const std::string file_name );

#endif
