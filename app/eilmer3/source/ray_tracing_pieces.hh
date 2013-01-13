/** \file ray_tracing_pieces.hh
 *  \ingroup eilmer3
 *
 *  \author Daniel F. Potter
 *  \version 23-Dec-09: Port from old lib/radiation
 *  \version 2011: Moved to app/eilmer3/source
 *
 **/

#ifndef RAY_TRACING_PIECES_HH
#define RAY_TRACING_PIECES_HH

#include <string>
#include <vector>

#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/geometry2/source/geom.hh"
#include "../../../lib/radiation/source/spectra_pieces.hh"

// NOTE: using face definitions from cell.hh for these
// #define NORTH  0
// #define EAST   1
// #define SOUTH  2
// #define WEST   3
// #define TOP    4
// #define BOTTOM 5
#define INSIDE_GRID    6
#define ERROR          7

class RayTracingPoint {
public:
    RayTracingPoint( Gas_data * Q, double * Q_rE_rad, double s, double vol, CoeffSpectra * X = 0,  BinnedCoeffSpectra * Y = 0, std::vector<double> * Q_rE_rad_temp = 0 );
    
    ~RayTracingPoint();
    
public:
    /* Pointers to CFD cell gas-data and radiative source term */
    Gas_data * Q_;
    double * Q_rE_rad_;
    
    /* 1D spatial location */
    double s_;
    
    /* Cell volume */
    double vol_;
    
    /* Pointer to spectral coefficients */
    CoeffSpectra * X_;
    
    /* Pointer to binned spectral coefficients */
    BinnedCoeffSpectra * Y_;

    /* Point to vector of doubles for threads to store source term contributions */
    std::vector<double> * Q_rE_rad_temp_;
};

class RayTracingRay {
public:
    RayTracingRay( double theta, double phi, double domega, Vector3 ray_origin = Vector3() );
    
    virtual ~RayTracingRay();
    
public:
    virtual Vector3 get_point_on_line( double &L ) = 0;

public:
    /* direction data */
    double theta_;
    double phi_;
    double domega_;
    
    /* cell origin data */
    Vector3 ray_origin_;
    
    /* intersecting point data */
    std::vector<RayTracingPoint*> points_;
    
    /* pointer to q_rad vector element where exit energy is dumped */
    double * q_rad_;
    
    /* interface area of exit surface element */
    double exit_area_;
    
    /* status flags */
    int status_;
    
    /* displacement from origin */
    double L_;
    
    /* remaining energy up grid exit */
    double E_exit_;
};

class RayTracingRay2D : public RayTracingRay {
public:
    RayTracingRay2D( double theta, double phi, double domega, Vector3 ray_origin = Vector3() );
    
    ~RayTracingRay2D();
    
public:
    Vector3 get_point_on_line( double &L );
};

class RayTracingRay3D : public RayTracingRay {
public:
    RayTracingRay3D( double theta, double phi, double domega, Vector3 ray_origin = Vector3() );
    
    ~RayTracingRay3D();
    
public:
    Vector3 get_point_on_line( double &L );
};

class RayTracingInterface;

class RayTracingCell {
public:
    RayTracingCell( Gas_data * Q, double * Q_rE_rad, Vector3 origin, double vol, double area=0.0 );
    
    virtual ~RayTracingCell();
    
public:    
    void recompute_spectra( RadiationSpectralModel * rsm );
    
    void set_CFD_cell_indices( int ii, int jj, int kk );
    
    void get_CFD_cell_indices( int &ii, int &jj, int &kk );
    
    std::string str();
    
    void write_rays_to_file( std::string filename );
    
public:
    /* CFD cell indices */
    int ii_, jj_, kk_;

    /* Pointers to CFD cell gas-data and radiative source term */
    Gas_data * Q_;
    double * Q_rE_rad_;
    
    /* Copy of origin cell coordinates */
    Vector3 origin_;
    
    /* Copy of origin cell volume */
    double vol_;
    
    /* Copy of origin cell area */
    double area_;
    
    /* Spectral coefficients */
    CoeffSpectra * X_;
    
    /* Binned spectral coefficients */
    BinnedCoeffSpectra * Y_;

    /* Ray-tracing rays */
    std::vector<RayTracingRay*> rays_;
    
    /* Ray-tracing interfaces */
    std::vector<RayTracingInterface*> interfaces_;
    
    /* Vector of doubles for threads to store source term contributions */
    std::vector<double> Q_rE_rad_temp_;
    
    /* Total amount of emitted radiant energy */
    double E_rad_;
};

class DiscreteTransferCell : public RayTracingCell {
public:
    DiscreteTransferCell( Gas_data * Q, double * Q_rE_rad, Vector3 origin, double vol, int nrays=0, int ndim=0, bool planar=false );
    
    virtual ~DiscreteTransferCell();
};

class MonteCarloCell : public RayTracingCell {
public:
    MonteCarloCell( Gas_data * Q, double * Q_rE_rad, Vector3 origin, double vol, double area );
    
    virtual ~MonteCarloCell();
};

class RayTracingInterface {
public:
    RayTracingInterface( Gas_data * Q, Vector3 origin, double area, double length );
    
    virtual ~RayTracingInterface();
    
public:    
    virtual void recompute_spectra( RadiationSpectralModel * rsm );
    
    void set_CFD_cell_indices( int ii, int jj, int kk );
    
    void get_CFD_cell_indices( int &ii, int &jj, int &kk );
    
    void write_rays_to_file( std::string filename );
    
public:
    /* Adjacent CFD interface indices */
    int ii_, jj_, kk_;
    
    /* Pointers to CFD interface gas-data and radiative flux term */
    Gas_data * Q_;
    double * q_rad_;
    
    /* Copy of origin interface coordinates */
    Vector3 origin_;
    
    /* Copy of origin interface area */
    double area_;
    
    /* Copy of origin interface length */
    double length_;
    
    /* Spectral coefficients */
    SpectralIntensity * S_;
    
    /* Binned spectral coefficients */
    BinnedSpectralIntensity * U_;

    /* Ray-tracing rays */
    std::vector<RayTracingRay*> rays_;
    
    /* Vector of doubles for threads to store (incident) radiative flux contributions */
    std::vector<double> q_rad_temp_;
    
    /* Total amount of emitted radiant energy */
    double E_rad_;
};

class DiscreteTransferInterface : public RayTracingInterface {
public:
    DiscreteTransferInterface( Gas_data * Q, Vector3 origin, double area, double length, int nrays=0, int ndim=0, bool planar=false );
    
    virtual ~DiscreteTransferInterface();
};

class MonteCarloInterface : public RayTracingInterface {
public:
    MonteCarloInterface( Gas_data * Q, Vector3 origin, double area, double length );
    
    virtual ~MonteCarloInterface();
};

#endif
