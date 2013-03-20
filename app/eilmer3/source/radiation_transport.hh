/** \file radiation_transport.hh
 *  \ingroup eilmer3
 *  \brief Header file for the radiation transport model class and functions.
 **/

#ifndef RADIATION_TRANSPORT_HH
#define RADIATION_TRANSPORT_HH

#include <string>
#include <vector>

#include "../../../lib/radiation/source/spectral_model.hh"
#include "../../../lib/util/source/lua_service.hh"
#include "../../../lib/util/source/randomc.h"

#include "block.hh"
#include "cell_finder.hh"
#include "ray_tracing_pieces.hh"

#define EXIT_ON_RT_FAILURE 1

class RadiationTransportModel {
public:
    RadiationTransportModel( void ) {}
    
    RadiationTransportModel( lua_State *L );
    
    virtual ~RadiationTransportModel();

    virtual std::string str() const = 0;
    
    virtual void compute_Q_rad_for_flowfield() = 0;
    
    virtual int initialise() = 0;
    
    int set_radiation_spectral_model( std::string file_name, size_t nthreads );
    
protected:
    std::vector<RadiationSpectralModel*> rsm_;
};

class NoRadiationTransportModel : public RadiationTransportModel {
public:
    NoRadiationTransportModel( void );
    
    ~NoRadiationTransportModel();

    std::string str() const;
    
    int initialise();

    void compute_Q_rad_for_flowfield();
};

class OpticallyThin : public RadiationTransportModel {
public:
    OpticallyThin( lua_State *L );
    
    ~OpticallyThin();

    std::string str() const;
    
    int initialise();
    
    void compute_Q_rad_for_flowfield();

    void compute_Q_rad_for_block( Block * A );

private:
    int spectrally_resolved_;
};

class TangentSlab : public RadiationTransportModel {
public:
    TangentSlab( lua_State *L );
    
    ~TangentSlab();
    
    std::string str() const;
    
    int initialise();
    
    void compute_Q_rad_for_flowfield();
    
    void compute_Q_rad_for_block( Block * A );

private:
    bool exact_;
};

class DiscreteTransfer : public RadiationTransportModel {
public:
    DiscreteTransfer( lua_State *L );
    
    ~DiscreteTransfer();
    
    std::string str() const;
    
    int initialise();
    
    void compute_Q_rad_for_flowfield();
    
private:
    void initialise_rays_for_cell( RayTracingCell * cell, size_t nrays, size_t ndim, bool planar );
    
    void initialise_rays_for_interface( RayTracingInterface * interface, size_t nrays, size_t ndim, bool planar );
    
    int trace_ray( RayTracingRay * ray, size_t ib, size_t ic, size_t jc, size_t kc );
    
    size_t get_cell_index( Block * A, size_t i, size_t j, size_t k )
    { return (k-A->kmin)*(A->nnj*A->nni)+(j-A->jmin)*A->nni+(i-A->imin); }
    
private:
    size_t nrays_;
    std::vector< std::vector<RayTracingCell*> > cells_;
    std::vector< std::vector<RayTracingInterface*> > interfaces_;
    CellFinder * cf_;
    double dl_lmin_ratio_;
    double dl_min_;
    double E_min_;
    int clustering_;
    int binning_;
    size_t N_bins_;
    std::vector<SpectralBin*> B_;
};

class MonteCarlo : public RadiationTransportModel {
public:
    MonteCarlo( lua_State *L );
    
    ~MonteCarlo();
    
    std::string str() const;
    
    int initialise();
    
    void compute_Q_rad_for_flowfield();
    
private:
    RayTracingRay * create_new_ray_for_cell( RayTracingCell * cell, size_t nrays, size_t ndim, bool planar );
    
    RayTracingRay * create_new_ray_for_interface( RayTracingInterface * interface, size_t nrays, size_t ndim, bool planar );
    
    int trace_ray_standard( RayTracingRay * ray, size_t ib, size_t ic, size_t jc, size_t kc, double nu, double E );
    
    int trace_ray_partitioned_energy( RayTracingRay * ray, size_t ib, size_t ic, size_t jc, size_t kc, double nu, double E );

    size_t get_cell_index( Block * A, size_t i, size_t j, size_t k )
    { return (k-A->kmin)*(A->nnj*A->nni)+(j-A->jmin)*A->nni+(i-A->imin); }
    
private:
    size_t nrays_;
    std::vector< std::vector<RayTracingCell*> > cells_;
    std::vector< std::vector<RayTracingInterface*> > interfaces_;
    CellFinder * cf_;
    double dl_lmin_ratio_;
    double dl_min_;
    double E_min_;
    int clustering_;
    int absorption_;
    CRandomMersenne *rg_;
};

RadiationTransportModel * create_radiation_transport_model( const std::string file_name );

#endif
