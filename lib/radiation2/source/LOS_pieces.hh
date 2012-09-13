/** \file LOS_pieces.hh
 *  \ingroup radiation
 *
 *  \brief Functions, classes and structures for line-of-sight calculations
 *
 *  \author Daniel F. Potter
 *  \version 09-Dec-08: initial implementation
 *           06-Jul-09: port from lib/radiation
 *
 **/

#ifndef LOS_PIECES_HH
#define LOS_PIECES_HH

#include <string>
#include <vector>

#include "../../gas/models/gas_data.hh"

#include "spectral_model.hh"
#include "spectra_pieces.hh"

class RadiatingPoint {
public:
    RadiatingPoint();
    
    RadiatingPoint( Gas_data * Q, double * Q_rE_rad, double s, double ds );
    
    ~RadiatingPoint();
    
    void redefine( Gas_data * Q, double * Q_rE_rad, double s, double ds );
    
public:
    /* Pointers to CFD cell gas-data and radiative source term */
    Gas_data * Q_;
    double * Q_rE_rad_;
    
    /* 1D spatial location */
    double s_;
    double ds_;
    
    /* Spectral coefficients */
    CoeffSpectra * X_;

    /* Binned spectral coefficients */
    BinnedCoeffSpectra * Y_;
};

class LOS_data {
public:
    LOS_data( RadiationSpectralModel * rsm=0, int nrps=0, double T_i=0.0, double T_f=0.0, bool store_spectra_on_disk=false );
    
    virtual ~LOS_data();
    
public:
    virtual void clear();
    
    void set_rad_point( int irp, Gas_data * Q, double * Q_rE_rad, double s, double ds );
    
    void clone_rad_point( int iprp, int irp, double * Q_rE_rad, double s, double ds );
    
    void create_spectral_bins( int bining_type, int N_bins, std::vector<SpectralBin*> & B );

    double integrate_LOS( SpectralIntensity &S );
    
    double integrate_LOS_with_binning( int binning_type, int N_bins );

    double integrate_LOS_MC( SpectralIntensity &S, int nphotons );

    void write_all_points_to_file( void );
    
    void write_point_to_file( int ip, std::string fname );
    
    void write_gas_data_to_file( void );
    
    RadiatingPoint * get_rpoint_pointer( int irp )
    { return rpoints_[irp]; }
    
    void load_spectra( int inu_start, int inu_end );

public:
    RadiationSpectralModel * rsm_;
    
    int nrps_;
    int nnus_;
    double T_i_;
    double T_f_;
    
    bool store_spectra_on_disk_;

    std::vector<RadiatingPoint*> rpoints_;
};

class TS_data : public LOS_data {
public:
    TS_data( RadiationSpectralModel * rsm=0, int nrps=0, double T_i=0.0, double T_f=0.0, bool store_spectra_on_disk=false );
    
    ~TS_data();
    
public:
    //       |                  |#
    // Shock |------------------|# Wall
    //       |                  |#
    
    void clear();

    double quick_solve_for_divq();
    
    double quick_solve_for_divq_in_blocks( int nblks );

    double quick_solve_for_divq_with_binning(  int binning_type, int N_bins );

    double exact_solve_for_divq();
    
    double quick_solve_for_planck_divq( double kappa_nu );

public:
    SpectralFlux * F_;
};

#endif
