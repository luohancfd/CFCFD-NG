/** \file spectra_pieces.hh
 *  \ingroup radiation2
 *
 *  \brief Functions, classes and structures for spectral calculations
 *
 *  \author Daniel F. Potter
 *  \version 15-Sep-09: initial implementation
 *
 **/

#ifndef SPECTRA_PIECES_HH
#define SPECTRA_PIECES_HH

#include <string>
#include <vector>

#define WAVELENGTH 0
#define WAVENUMBER 1
#define FREQUENCY  2
#define ENERGY     3

// Forward declaration of RadiationSpectralModel
class RadiationSpectralModel;

static std::vector<double> zero_vec;

class SpectralContainer {
public:
    /// \brief Minimal Constructor 
    SpectralContainer();
    
    /// \brief Constructor 
    SpectralContainer( RadiationSpectralModel * rsm );
    
    /// \brief Copy constructor 
    SpectralContainer( SpectralContainer &C );
    
    /// \brief Deconstructor
    virtual ~SpectralContainer() = 0;
    
public:    
    virtual double write_to_file( std::string fname, int spectral_units ) = 0;
    
    double write_data_to_file( std::string fname, int spectral_units,
    		    	       std::vector<double> &Y1, std::string Y1_label, std::string Y1_int_label,
    		    	       std::vector<double> &Y2 = zero_vec, std::string Y2_label = "" );
    
public:
    std::vector<double> nu;
};

class CoeffSpectra : public SpectralContainer {
public:
    /// \brief Minimal Constructor 
    CoeffSpectra();
    
    /// \brief Constructor 
    CoeffSpectra( RadiationSpectralModel * rsm );
    
    /// \brief Deconstructor
    ~CoeffSpectra();
    
public:
    /* Spectral properties storage */
    std::vector<double> j_nu;
    std::vector<double> kappa_nu;
    
    /* Cumulative emission vector */
    std::vector<double> j_int;
    
public:
    /// \brief Clone function
    CoeffSpectra * clone();
    
    /// \brief Write CoeffSpectra class data to file
    double write_to_file( std::string fname, int spectral_units=WAVELENGTH );
    
    /// \brief Write CoeffSpectra class data to file in format for TRT_tools interface
    void write_TRT_tools_file( std::string fname );

    /// \brief Integrate the emission coefficient spectrum
    double integrate_emission_spectra();
    
    /// \brief Calculate and store the cumulative emission coefficient
    void calculate_cumulative_emission( bool resize=false );

    /// \brief Determine a random frequency for this spectra via the Monte Carlo method
    double random_frequency( double R );

    /// \brief Determine a random frequency interval for this spectra via the Monte Carlo method
    int random_frequency_interval( double R );

    /// \brief Calculate the absorption coefficient for the given frequency (via interpolation if necessary)
    double kappa_from_nu( double nu );
};

#define NO_BINNING        0
#define FREQUENCY_BINNING 1
#define OPACITY_BINNING   2

class SpectralBin {
public:
    SpectralBin( std::vector<double> & pvec, double p_min, double p_max );
    SpectralBin( std::vector<int> & inus );

public:
    std::vector<int> inu;
};

void create_spectral_bin_vector( std::vector<double> & pvec, int binning_type, int N_bins, std::vector<SpectralBin*> & B );

class BinnedCoeffSpectra {
public:
    /// \brief Minimal Constructor
    BinnedCoeffSpectra( CoeffSpectra * X, std::vector<SpectralBin*> & B );

    /// \brief Deconstructor
    ~BinnedCoeffSpectra();

public:
    double sum_emission();

public:
    std::vector<double> kappa_bin;
    std::vector<double> j_bin;
};

double eval_Gaussian( double x, double delta_x );

class SpectralIntensity : public SpectralContainer {
public:
    /// \brief Minimal Constructor 
    SpectralIntensity();
    
    /// \brief Constructor from rsm
    SpectralIntensity( RadiationSpectralModel * rsm  );
    
    /// \brief Constructor from rsm with plank function evaluated at T
    SpectralIntensity( RadiationSpectralModel * rsm, double T );
    
    /// \brief Copy constructor
    SpectralIntensity( SpectralIntensity &S  );
    
    /// \brief Deconstructor
    ~SpectralIntensity();
    
public:
    double write_to_file( std::string fname, int spectral_units=WAVELENGTH );

    void apply_apparatus_function( double delta_x_ang, int nu_sample=1 );

    void reverse_data_order();
    
    void reset_intensity_vector();
    
    double integrate_intensity_spectra( double lambda_min=-1.0, double lambda_max=-1.0 );
    
    /// \brief Determine a random frequency for this spectra via the Monte Carlo method
    double random_frequency( double R );

    /// \brief Determine a random frequency interval for this spectra via the Monte Carlo method
    int random_frequency_interval( double R );

public:
    std::vector<double> I_nu;
    std::vector<double> I_int;
    int nwidths;
};

class BinnedSpectralIntensity {
public:
    /// \brief Minimal Constructor
    BinnedSpectralIntensity( size_t N_bins );

    BinnedSpectralIntensity( SpectralIntensity * I, std::vector<SpectralBin*> & B );

    BinnedSpectralIntensity( RadiationSpectralModel * rsm, double T, std::vector<SpectralBin*> & B );

    /// \brief Deconstructor
    ~BinnedSpectralIntensity();

public:
    double sum_intensity();

public:
    std::vector<double> I_bin;
};

class SpectralFlux : public SpectralContainer {
public:
    /// \brief Minimal Constructor 
    SpectralFlux();
    
    /// \brief Constructor 
    SpectralFlux( RadiationSpectralModel * rsm  );
    
    /// \brief Constructor from temperature
    SpectralFlux( RadiationSpectralModel * rsm, double T );

    /// \brief Deconstructor
    ~SpectralFlux();
    
public:
    double write_to_file( std::string fname, int spectral_units=WAVELENGTH );
    
public:
    std::vector<double> q_nu;
};

class BinnedSpectralFlux {
public:
    /// \brief Minimal Constructor
    BinnedSpectralFlux( size_t N_bins );

    BinnedSpectralFlux( SpectralFlux * F, std::vector<SpectralBin*> & B );

    BinnedSpectralFlux( RadiationSpectralModel * rsm, double T, std::vector<SpectralBin*> & B );

    /// \brief Deconstructor
    ~BinnedSpectralFlux();

public:
    double sum_flux();

public:
    std::vector<double> q_bin;
};


class IntensityProfile {
public:
    /// \brief Minimal Constructor 
    IntensityProfile();
    
    /// \brief Deconstructor
    ~IntensityProfile();
    
public:
    void add_new_point( double x, double I );
    
    void spatially_smear( double dx_smear );
    
    void spatially_smear_for_varying_dx( double dx_smear );
    
    void write_to_file( std::string fname );
    
public:
    std::vector<double> x_vec;
    std::vector<double> I_vec;
};

class SpectralField {
public:
    /// \brief Minimal Constructor 
    SpectralField();
    
    /// \brief Deconstructor
    ~SpectralField();
    
public:
    void add_new_intensity_spectra( double x, SpectralIntensity * S );
    
    IntensityProfile extract_intensity_profile( double lambda_l, double lambda_u );
    
    SpectralIntensity * last_spectra()
    { return S_vec.back(); }
    
public:
    std::vector<double> x_vec;
    std::vector<SpectralIntensity*> S_vec;
};

/* Some other useful functions */

/// \brief Convert from frequency (Hz) to wavelength (nm)
double nu2lambda(double nu);

/// \brief Convert from wavelength (nm) to frequency (Hz)
double lambda2nu(double lambda_nm);

/// \brief Calculate Planck (blackbody) intensity at frequency nu (Hz) and temperature T (K)
double planck_intensity(const double nu, const double T);

/// \brief Get the frequency index in 'nus' that is just below 'nu'
int get_nu_index( std::vector<double> &nus, double nu );

#endif
