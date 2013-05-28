/** \file spectra_pieces.cxx
 *  \ingroup radiation
 *
 *  \brief Class definitions for line-of-sight calculations
 *
 *  \author Daniel F. Potter
 *  \version 15-Sep-09: initial implementation
 *            
 **/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <algorithm>

#include "../../nm/source/exponential_integrals.hh"
#include "../../util/source/useful.h"

#include "spectral_model.hh"
#include "spectra_pieces.hh"
#include "radiation_constants.hh"
#include "line_shapes.hh"

using namespace std;

/* ------------ SpectralContainer class ------------ */

SpectralContainer::SpectralContainer()
 : nwidths( NWIDTHS ), adaptive( false ) {}

SpectralContainer::SpectralContainer( RadiationSpectralModel * rsm )
{
    // Set nwidths to the default value
    nwidths = NWIDTHS;

    // Set the adaptive flag from the spectral model
    adaptive = rsm->get_adaptive_spectral_grid();

    // Fill out the spectral grid
    // WARNING: if adaptive, this will use the currently set linewidths
    rsm->radiative_spectral_grid( nu );
}

SpectralContainer::SpectralContainer(SpectralContainer &C)
: nu( C.nu ), nwidths ( C.nwidths ), adaptive( C.adaptive ) {}

SpectralContainer::~SpectralContainer()
{
    nu.resize(0);
}

double
SpectralContainer::
write_data_to_file( string fname, int spectral_units,
    		    vector<double> &Y1, string Y1_label, string Y1_int_label,
    		    vector<double> &Y2, string Y2_label )
{
    /* 0. Determine if Y2 data is present */
    bool with_Y2 = ( Y2.size() > 0 ) ? true : false;
    
    /* 1. Setup the output file. */
    ofstream specfile;
    specfile.open(fname.c_str());
    if( specfile.fail() ) {
	cout << "Error opening file: " << fname << endl;
	cout << "Bailing Out!\n";
	exit(FILE_ERROR);
    }
    
    specfile << setprecision(12) << scientific << showpoint;

    specfile << "# " << fname << endl;
    if ( spectral_units==WAVELENGTH )
        specfile << "# Column 1: Wavelength (nm)" << endl;
    else if ( spectral_units==WAVENUMBER )
        specfile << "# Column 1: Wavenumber (cm-1)" << endl;
    else if ( spectral_units==FREQUENCY )
        specfile << "# Column 1: Frequency (Hz)" << endl;
    else if ( spectral_units==ENERGY )
        specfile << "# Column 1: Energy (eV)" << endl;
    else {
        cout << "SpectralContainer::write_data_to_file()" << endl
             << "spectral unit: " << spectral_units << "not recognised" << endl;
        exit( BAD_INPUT_ERROR );
    }
    specfile << "# Column 2: " << Y1_label << endl
             << "# Column 3: " << Y1_int_label << endl;
    if ( with_Y2 ) {
    	specfile << "# Column 4: " << Y2_label << endl;
    }
    /* 2. Write data for each frequency interval to file. */
    double Y1_int = 0.0;
    if ( spectral_units==WAVELENGTH ) {
        for ( int inu=int(nu.size())-1; inu >= 0; --inu ) {
            // Integration of Y1
            if ( inu < int(nu.size())-1 )
                Y1_int += 0.5 * ( Y1[inu] + Y1[inu+1] ) * fabs(nu[inu] - nu[inu+1]);
            // Spectral coordinate
            double SC = nu2lambda( nu[inu] );
            // Write to file
            specfile << setw(20) << SC
                     << setw(20) << Y1[inu] * nu[inu] * nu[inu] / RC_c_SI
                     << setw(20) << Y1_int;
            if ( with_Y2 ) {
                specfile << setw(20) << Y2[inu];
            }
            specfile << endl;
            // move on to next frequency interval
        }
    }
    else {
        for ( int inu=0; inu < int(nu.size()); ++inu ) {
            // Integration of Y1
            if ( inu > 0 )
                Y1_int += 0.5 * ( Y1[inu-1] + Y1[inu] ) * fabs(nu[inu-1] - nu[inu]);
            // Spectral coordinate
            double SC = 0.0;
            if ( spectral_units==WAVENUMBER )
                SC = 1.0 / nu2lambda( nu[inu] ) * 1.0e7;
            else if ( spectral_units==FREQUENCY )
                SC = nu[inu];
            else if ( spectral_units==ENERGY )
                SC = RC_h_SI * nu[inu] / RC_e_SI;
            // Write to file
            specfile << setw(20) << SC
                     << setw(20) << Y1[inu] * nu[inu] / SC
                     << setw(20) << Y1_int;
            if ( with_Y2 ) {
                specfile << setw(20) << Y2[inu];
            }
            specfile << endl;
            // move on to next frequency interval
        }
    }
    
    specfile.close();
    
    return Y1_int;
}

/* ------------ CoeffSpectra class ------------ */

CoeffSpectra::CoeffSpectra() {}

CoeffSpectra::CoeffSpectra( RadiationSpectralModel * rsm )
 : SpectralContainer( rsm )
{
    j_nu.resize( nu.size(), 0.0 );
    kappa_nu.resize( nu.size(), 0.0 );
    j_int.resize( nu.size(), 0.0 );
}

CoeffSpectra::CoeffSpectra( CoeffSpectra * X )
{
    // FIXME: use C++ vector copying functions here
    nu.resize( X->nu.size() );
    j_nu.resize( nu.size() );
    kappa_nu.resize( nu.size() );
    j_int.resize( nu.size() );

    for ( size_t inu=0; inu<X->nu.size(); ++inu) {
        nu[inu] = X->nu[inu];
        j_nu[inu] = X->j_nu[inu];
        kappa_nu[inu] = X->kappa_nu[inu];
        j_int[inu] = X->j_int[inu];
    }

    nwidths = X->nwidths;
    adaptive = X->adaptive;
}

CoeffSpectra::~CoeffSpectra()
{
    j_nu.resize(0);
    kappa_nu.resize(0);
    j_int.resize(0);
}

CoeffSpectra * CoeffSpectra::clone()
{
    CoeffSpectra * X = new CoeffSpectra();
    X->nu.assign( nu.begin(), nu.end() );
    X->j_nu.assign( j_nu.begin(), j_nu.end() );
    X->kappa_nu.assign( kappa_nu.begin(), kappa_nu.end() );
    X->j_int.assign( j_int.begin(), j_int.end() );
    
    X->nwidths = nwidths;
    X->adaptive = adaptive;

    return X;
}

void CoeffSpectra::clear_data()
{
    nu.clear();
    j_nu.clear();
    kappa_nu.clear();
    j_int.clear();
}

double CoeffSpectra::write_to_file( string fname, int spectral_units )
{
    string Y1_label = "Emission coefficient, j_lambda (W/m**3-sr-m)";
    if ( spectral_units==WAVENUMBER )
        Y1_label = "Emission coefficient, j_eta (W/m**3-sr-1/cm)";
    else if ( spectral_units==FREQUENCY )
        Y1_label = "Emission coefficient, j_nu (W/m**3-sr-Hz)";
    else if ( spectral_units==ENERGY )
        Y1_label = "Emission coefficient, j_epsilon (W/m**3-sr-eV)";
    string Y1_int_label = "Integrated emission, j (W/m**3-sr)";
    string Y2_label = "Absorption coefficient, kappa_lambda (1/m)";
    
    return write_data_to_file( fname, spectral_units, j_nu, Y1_label, Y1_int_label, kappa_nu, Y2_label );
}

void CoeffSpectra::read_from_file( string fname, int inu_start, int inu_end )
{
    // make sure this coeffspectra is clear
    this->clear_data();

    string line;
    ifstream specfile (fname.c_str());
    if (specfile.is_open()) {
        int spectral_units = -1;
        getline (specfile,line); // should be the filename
        getline (specfile,line); // should be spectral column descriptor
        if ( line.compare(0,11,"# Column 1:")==0 ) {
            if ( line.compare(12,14,"Frequency (Hz)")==0 ) {
                spectral_units = FREQUENCY;
            }
            else if ( line.compare(12,15,"Wavelength (nm)")==0 ) {
                spectral_units = WAVELENGTH;
            }
            else {
                cout << "SpectralFlux::read_from_file()" << endl
                     << "Only frequency and wavelength units are currently supported" << endl;
                exit(FAILURE);
            }
        }
        else {
            cout << "Expected the string 'Column 1' in this line:" << endl
                 << line << endl;
            exit(FAILURE);
        }
        getline (specfile,line); // should be emission coefficient column descriptor
        getline (specfile,line); // should be abs coefficient column descriptor
        getline (specfile,line); // should be int emission coefficient column descriptor
        int count = 0;
        while ( specfile.good() ) {
            double wav, j_wav, kappa_wav, _j_int;
            specfile >> wav >> j_wav >> _j_int >> kappa_wav;
            if ( count >= inu_start && ( count <= inu_end || inu_end < 0 ) ) {
                double _nu = 0.0, _j_nu = 0.0, _kappa_nu = 0.0;
                if ( spectral_units==FREQUENCY ) {
                    _nu = wav;
                    _j_nu = j_wav;
                    _kappa_nu = kappa_wav;
                }
                else if ( spectral_units==WAVELENGTH ) {
                    _nu = nu2lambda(wav);
                    _j_nu = j_wav * wav * 1.0e-9 / _nu;
                    _kappa_nu = kappa_wav;
                }
                nu.push_back(_nu);
                j_nu.push_back(_j_nu);
                kappa_nu.push_back(_kappa_nu);
                // not storing j_int, calculating later
                // j_int.push_back(_j_int);
            }
            count++;
        }
        specfile.close();
    }
    else {
        cout << "CoeffSpectra::read_from_file()" << endl
             << "Unable to open file with name: " << fname << endl;
        exit(BAD_INPUT_ERROR);
    }

    if ( inu_end<0 ) {
        // remove the last entry which is a double-up
        nu.erase(nu.end()-1);
        j_nu.erase(j_nu.end()-1);
        kappa_nu.erase(kappa_nu.end()-1);
    }

    cout << "Read " << nu.size() << " spectral points from file: " << fname << endl;
    double j_total = this->integrate_emission_spectra();
    cout << "Integral is " << j_total << endl;
}

void CoeffSpectra::write_TRT_tools_file( string fname )
{
    ofstream specfile;
    specfile.open(fname.c_str());
    if( specfile.fail() ) {
        cout << "Error opening file: " << fname << endl;
        cout << "Bailing Out!\n";
        exit(FILE_ERROR);
    }

    specfile << setprecision(12) << scientific << showpoint;

    // Column 1: Wavelength (Ang)
    // Column 2: Emission coefficient, j_lambda (W/m3-m-sr)
    // Column 3: Absorption coefficient, kappa_lambda (1/m)
    for ( int inu=int(nu.size())-1; inu >= 0; --inu ) {
        // Write to file
        specfile << setw(20) << nu2lambda( nu[inu] ) * 10
                 << setw(20) << j_nu[inu] * nu[inu] * nu[inu] / RC_c_SI
                 << setw(20) << kappa_nu[inu]
                 << endl;
        // move on to next frequency interval
    }

    specfile.close();

    return;
}

double CoeffSpectra::integrate_emission_spectra()
{
    double j_total = 0.0;
    for ( int inu=1; inu < int(nu.size()); inu++ ) {
    	j_total += 0.5 * ( j_nu[inu] + j_nu[inu-1] ) * fabs(nu[inu] - nu[inu-1]);
    }

    return j_total;
}

void CoeffSpectra::calculate_cumulative_emission( bool resize )
{
    // 0. (Re)size the vector if requested
    if ( resize )
    	j_int.resize( j_nu.size(), 0.0 );
    
    // 1. Integrate, storing cumulative integrand along the way
    j_int[0] = 0.0;
    for ( int inu=1; inu < int(nu.size()); inu++ ) {
        j_int[inu] = j_int[inu-1] + 0.5 * ( j_nu[inu] + j_nu[inu-1] ) * fabs(nu[inu] - nu[inu-1]);
    }

    return;
}

#define QUADRATIC_EQUATION(a,b,c) ( -b + sqrt( b*b - 4.0 * a * c ) ) / ( 2.0 * a )

double CoeffSpectra::random_frequency( double R )
{
    double j_interval = R * j_int.back();
    // cout << "j_interval = " << j_interval << ", R = " << R << ", j_int.back() = " << j_int.back() << endl;
    int inu_b = 0, inu_u = int( j_int.size() - 1 ), inu_mp = 0;
    double f_b = 0.0, f_mp;
    while ( abs( inu_u - inu_b ) > 1 ) {
        f_b = j_int[inu_b] - j_interval;
        inu_mp = int( ( inu_u + inu_b ) / 2.0 );
        f_mp = j_int[inu_mp] - j_interval;
        ( f_b * f_mp > 0.0 ) ? ( inu_b = inu_mp ) : ( inu_u = inu_mp );
    }

    // now need to solve for the actual frequency
    double m = ( j_nu[inu_u] - j_nu[inu_b] ) / ( nu[inu_u] - nu[inu_b] );
    double nu_interval = QUADRATIC_EQUATION( (0.5*m), (-m*nu[inu_b]+j_nu[inu_b]), (0.5*m*nu[inu_b]*nu[inu_b] - j_nu[inu_b]*nu[inu_b] - ( j_interval - j_int[inu_b] ) ) );

    // cout << "CoeffSpectra::random_frequency()" << endl;
    // cout << "inu_b = " << inu_b << ", inu_u = " << inu_u << endl;
    // cout << "nu = " << nu[inu_b] + dnu << ", nu[inu_b] = " << nu[inu_b] << ", nu[inu_u] = " << nu[inu_u] << endl;
    // cout << "j_interval = " << j_interval << ", j_int[inu_b] = " << j_int[inu_b] << ", j_int[inu_u] = " << j_int[inu_u] << endl;

    // if ( nu_interval > nu[inu_u] ) {
    //     cout << "CoeffSpectra::random_frequency()" << endl
    //          << "nu_interval = " << nu_interval << " > nu[inu_u] = " << nu[inu_u] << endl;
    //     exit( FAILURE );
    // }
    // else if ( nu_interval < nu[inu_b] ) {
    //     cout << "CoeffSpectra::random_frequency()" << endl
    //          << "nu_interval = " << nu_interval << " < nu[inu_b] = " << nu[inu_b] << endl;
    //     exit( FAILURE );
    // }
    return nu_interval;
}

#undef QUADRATIC_EQUATION

int CoeffSpectra::random_frequency_interval( double R )
{
    double j_interval = R * j_int.back();
    // cout << "j_interval = " << j_interval << ", R = " << R << ", j_int.back() = " << j_int.back() << endl;
    int inu_b = 0, inu_u = int( j_int.size() - 1 ), inu_mp = 0;
    double f_b = 0.0, f_mp;
    while ( abs( inu_u - inu_b ) > 1 ) {
        f_b = j_int[inu_b] - j_interval;
        inu_mp = int( ( inu_u + inu_b ) / 2.0 );
        f_mp = j_int[inu_mp] - j_interval;
        ( f_b * f_mp > 0.0 ) ? ( inu_b = inu_mp ) : ( inu_u = inu_mp );
    }
    f_b = j_int[inu_b] - j_interval;
    double f_u = j_int[inu_u] - j_interval;

    return fabs(f_b) < fabs(f_u) ? inu_b : inu_u;
}

double CoeffSpectra::kappa_from_nu( double nu_interval )
{
    if ( nu.size()==1 ) return kappa_nu[0];

    int inu_b = 0, inu_u = int( nu.size() - 1 ), inu_mp = 0;
    double f_b = 0.0, f_mp;
    while ( abs( inu_u - inu_b ) > 1 ) {
        f_b = nu[inu_b] - nu_interval;
        inu_mp = int( ( inu_u + inu_b ) / 2.0 );
        f_mp = nu[inu_mp] - nu_interval;
        ( f_b * f_mp > 0.0 ) ? ( inu_b = inu_mp ) : ( inu_u = inu_mp );
    }

    double dnu = nu_interval - nu[inu_b];
    double m = ( kappa_nu[inu_u] - kappa_nu[inu_b] ) / ( nu[inu_u] - nu[inu_b] );
    double kappa = kappa_nu[inu_b] + m * dnu;

    // cout << "CoeffSpectra::kappa_from_nu()" << endl;
    // cout << "inu_b = " << inu_b << ", inu_u = " << inu_u << endl;
    // cout << "nu_interval = " << nu_interval << ", nu[inu_b] = " << nu[inu_b] << ", nu[inu_u] = " << nu[inu_u] << endl;
    // cout << "kappa = " << kappa << ", kappa_nu[inu_b] = " << kappa_nu[inu_b] << ", kappa_nu[inu_u] = " << kappa_nu[inu_u] << endl;

    return kappa;
}

void CoeffSpectra::coeffs_from_nu( double nu_star, double &j_nu_star, double &kappa_nu_star )
{
    if ( nu.size()==1 ) {
        j_nu_star = j_nu[0];
        kappa_nu_star = kappa_nu[0];
        return;
    }

    int inu_b = 0, inu_u = int( nu.size() - 1 ), inu_mp = 0;
    double f_b = 0.0, f_mp;
    while ( abs( inu_u - inu_b ) > 1 ) {
        f_b = nu[inu_b] - nu_star;
        inu_mp = int( ( inu_u + inu_b ) / 2.0 );
        f_mp = nu[inu_mp] - nu_star;
        ( f_b * f_mp > 0.0 ) ? ( inu_b = inu_mp ) : ( inu_u = inu_mp );
    }

    double dnu = nu_star - nu[inu_b];
    double m;
    // Absorption coefficient
    m = ( kappa_nu[inu_u] - kappa_nu[inu_b] ) / ( nu[inu_u] - nu[inu_b] );
    kappa_nu_star = kappa_nu[inu_b] + m * dnu;
    // Emission coefficient
    m = ( j_nu[inu_u] - j_nu[inu_b] ) / ( nu[inu_u] - nu[inu_b] );
    j_nu_star = j_nu[inu_b] + m * dnu;

    return;
}

void CoeffSpectra::apply_apparatus_function( ApparatusFunction * A  )
{
    // Quick exit if the representative width is too small
    if ( A->gamma_star < 1.0e-10 ) return;

    // Initialise the Apparatus function
    A->initialise();

    // A vector to temporarily hold smeared data
    vector<double> nu_temp;
    int count = 0;
    for( size_t inu=0; inu<nu.size(); inu++) {
	count++;
        if ( count==A->nu_sample ) {
            nu_temp.push_back(nu[inu]);
            count = 0;
        }
    }
    vector<double> j_nu_temp( nu_temp.size() );
    vector<double> kappa_nu_temp( nu_temp.size() );

    int percentage=0;
    for( size_t inu=0; inu<nu_temp.size(); inu++) {
	double nu_val = nu_temp[inu];
	double lambda_ang = 10.0 * nu2lambda( nu_val );
	// convert HWHM's to Hz
	double gamma_star_Hz = A->gamma_star / lambda_ang * nu_val;
	double nu_lower = nu_val - double(nwidths) * gamma_star_Hz;
	double nu_upper = nu_val + double(nwidths) * gamma_star_Hz;
	int jnu_start = get_nu_index(nu,nu_lower,adaptive) + 1;
	int jnu_end = get_nu_index(nu,nu_upper,adaptive) + 1;

	if((double(inu)/double(nu_temp.size())*100.0+1.0)>double(percentage)) {
	    cout << "Smearing spectrum: | " << percentage << "% |, jnu_start = " << jnu_start << ", jnu_end = " << jnu_end << " \r" << flush;
	    percentage += 10;
        }

	// Apply convolution integral over this frequency range with trapezoidal method
	double j_nu_conv = 0.0;
	double kappa_nu_conv = 0.0;
	double AF_integral = 0.0;
	for ( int jnu=jnu_start+1; jnu<jnu_end; jnu++ ) {
            double dnu0 = nu[jnu-1] - nu_val;
            double dnu1 = nu[jnu] - nu_val;
            double f_nu0 = A->eval(nu_val, dnu0);
            double f_nu1 = A->eval(nu_val, dnu1);
            j_nu_conv += 0.5 * ( j_nu[jnu-1]*f_nu0 + j_nu[jnu]*f_nu1 ) * ( nu[jnu] - nu[jnu-1] );
            kappa_nu_conv += 0.5 * ( kappa_nu[jnu-1]*f_nu0 + kappa_nu[jnu]*f_nu1 ) * ( nu[jnu] - nu[jnu-1] );
            AF_integral += 0.5 * ( f_nu0 + f_nu1 ) * ( nu[jnu] - nu[jnu-1] );
	}

	// cout << "AF_integral = " << AF_integral << endl;

	// Make sure a zero value is not returned if nwidths is too small
	if ( jnu_start==(jnu_end-1) ) {
	    cout << "SpectralIntensity::apply_apparatus_function()" << endl
	         << "WARNING: nwidths is too small!" << endl;
	    j_nu_conv = j_nu[inu];
	    kappa_nu_conv = kappa_nu[inu];
	}

	// Save result (and rescale to ensure the integral remains the same)
	j_nu_temp[inu] = j_nu_conv / AF_integral;
	kappa_nu_temp[inu] = kappa_nu_conv / AF_integral;
    }

    cout << endl;

    // Loop again to overwrite old I_nu values in S
    nu.resize(nu_temp.size());
    j_nu.resize(nu_temp.size());
    kappa_nu.resize(nu_temp.size());
    j_int.resize(nu_temp.size());
    for( size_t inu=0; inu<nu_temp.size(); inu++) {
	nu[inu] = nu_temp[inu];
    	j_nu[inu] = j_nu_temp[inu];
    	kappa_nu[inu] = kappa_nu_temp[inu];
    	// recompute cumulative emissivity
    	if ( inu==0.0 )
    	    j_int[inu] = 0.0;
    	else
    	    j_int[inu] = 0.5 * ( j_nu[inu] + j_nu[inu-1] ) * ( nu[inu] - nu[inu-1] );
    }

    return;
}

/* ------------ SpectralBin class ------------ */

SpectralBin::SpectralBin(vector<double> & pvec, double p_min, double p_max )
{
    for ( int ip=0; ip<int(pvec.size()); ++ip ) {
	if ( pvec[ip] >= p_min && pvec[ip] < p_max ) {
	    inu.push_back( ip );
	}
    }
}

SpectralBin::SpectralBin(vector<int> & inus )
 : inu( inus ) {}

int create_spectral_bin_vector( std::vector<double> & pvec, int binning_type, int N_bins, std::vector<SpectralBin*> & B )
{
    if ( N_bins<1 ) {
	cout << "create_spectral_bin_vector()" << endl
	     << "N_bins must be greater than or equal to 1" << endl;
	exit( BAD_INPUT_ERROR );
    }

    if ( binning_type==FREQUENCY_BINNING ) {
	// Simply subdivide spectral range into N_bins of equal size
	int nnus = (int) pvec.size();
	int bin_size = nnus / N_bins;
	int first_bin_size = nnus - ( N_bins - 1 ) * bin_size;
	int inu_start = 0, inu_end = first_bin_size;
	for ( int iB=0; iB<N_bins; ++iB ) {
	    vector<int> inus;
	    for ( int inu=inu_start; inu<inu_end; ++inu )
		inus.push_back( inu );
	    inu_start = inu_end;
	    inu_end += bin_size;
	    B.push_back( new SpectralBin( inus ) );
	}
	// cout << "inu_end = " << inu_end << ", bin_size = " << bin_size << ", first_bin_size = " << first_bin_size << endl;
	// cout << "nnus = " << nnus << ", N_bins = " << N_bins << endl;
    }
    else if ( binning_type==OPACITY_BINNING ) {
	// divide the opacity range into equal segment in log space
	int nps = (int) pvec.size();
	double p_min = 1.0e99, p_max = 0.0;
	for ( int ip=0; ip<nps; ++ip ) {
	    if ( pvec[ip] < p_min && pvec[ip]>0.0 ) p_min = pvec[ip];
	    if ( pvec[ip] > p_max ) p_max = pvec[ip];
	}
	if ( p_min > p_max ) {
	    cout << "create_spectral_bin_vector()" << endl
		 << "p_min and p_max search failed." << endl
		 << "Exiting program." << endl;
	    exit( FAILURE );
	}
	else if ( p_min == p_max && N_bins != 1 ) {
	    cout << "create_spectral_bin_vector()" << endl
		 << "It appears this is a gray gas as the opacity is constant." << endl
		 << "For this type of radiation, spectral binning will only work with a single bin." << endl
		 << "Currently N_bins = " << N_bins << "." << endl
		 << "Exiting program." << endl;
	    exit( BAD_INPUT_ERROR );
	}
	else if ( p_min == p_max && N_bins == 1 ) {
	    // We need to artificially increase p_max so that all spectral points are included in the bin
	    // Any increase is sufficient as all opacities are equal
	    p_max += 1.0;
	}
	vector<double> p_limits;
	double delta_log_p = log(p_max / p_min ) / ( N_bins );
	for ( int iB=0; iB<N_bins+1; ++iB ) {
	   double log_p = log(p_min) + iB*delta_log_p;
	   p_limits.push_back(exp(log_p));
	   // cout << "p_limits = " << p_limits.back() << endl;
	}
	// decrease the first limit, and increase the last limit to ensure all points will be included
        p_limits.front() *= 0.99;
        p_limits.back() *= 1.01;
	int nnu_inc = 0.0;
	for ( int iB=1; iB<N_bins+1; ++iB ) {
	    SpectralBin * SB = new SpectralBin( pvec, p_limits[iB-1], p_limits[iB] );
	    nnu_inc +=  SB->inu.size();
	    if ( SB->inu.size()==0 ) {
	        delete SB;      // is this necessary?
	        continue;
	    }
	    B.push_back( SB );
	    // cout << "Created a new spectral bin with " << B.back()->inu.size() << " entries" << endl;
	}
	// cout << nnu_inc << " spectral points have been included from a total of " << nps << endl;
    }

    return (int) B.size();
}

/* ------------ BinnedCoeffSpectra class ------------ */

BinnedCoeffSpectra::BinnedCoeffSpectra( CoeffSpectra * X, vector<SpectralBin*> & B )
{
    int nnu = X->nu.size();

    // 1. Calculate the kappa_bin and j_bin vectors
    for ( int iB=0; iB<int(B.size()); ++iB ) {
	double j_int = 0.0, S_int = 0.0;
	for ( int j=0; j<int(B[iB]->inu.size()); ++j) {
	    int inu = B[iB]->inu[j];
	    // NOTE: we are taking care here to remain consistent with a trapezoidal discretisation of the spectra
	    //       where the trapezoid heights are taken as the average of the two bounding points
	    double dj_int = 0.0;
	    if      ( inu==0 )     dj_int = X->j_nu[inu] * 0.5 * fabs( X->nu[inu+1] - X->nu[inu] );
	    else if ( inu==nnu-1 ) dj_int = X->j_nu[inu] * 0.5 * fabs( X->nu[inu] - X->nu[inu-1] );
	    else                   dj_int = X->j_nu[inu] * ( 0.5 * fabs( X->nu[inu] - X->nu[inu-1] ) + 0.5 * fabs( X->nu[inu+1] - X->nu[inu] ) );
	    j_int += dj_int;
	    if ( X->kappa_nu[inu] > 0.0 )
	        S_int += dj_int / X->kappa_nu[inu];
	}
	j_bin.push_back( j_int );
	// If the source function is zero then kappa should also be zero
	if ( S_int==0.0 )
	    kappa_bin.push_back( 0.0 );
	else
	    kappa_bin.push_back( j_int / S_int );
    }

}

BinnedCoeffSpectra::~BinnedCoeffSpectra()
{}

double BinnedCoeffSpectra::sum_emission()
{
    double j_total = 0.0;
    for ( size_t ib=0; ib<j_bin.size(); ++ib )
	j_total += j_bin[ib];
    return j_total;
}

int BinnedCoeffSpectra::random_bin( double R )
{
    double j_interval = R * this->sum_emission();
    double j_sum = 0.0;
    for ( size_t ib=0; ib<j_bin.size(); ++ib ) {
        j_sum += j_bin[ib];
        if ( j_sum > j_interval ) return (int) ib;
    }

    // If we get here something has gone wrong (R probably was greater than 1)
    cout << "BinnedCoeffSpectra::random_bin() failed!" << endl;
    exit( FAILURE );
}

void BinnedCoeffSpectra::write_TRT_tools_file( string fname )
{
    ofstream specfile;
    specfile.open(fname.c_str());
    if( specfile.fail() ) {
        cout << "Error opening file: " << fname << endl;
        cout << "Bailing Out!\n";
        exit(FILE_ERROR);
    }

    specfile << setprecision(12) << scientific << showpoint;

    // Column 1: Binned emission coefficient, j_bin (W/m3-sr)
    // Column 2: Binned absorption coefficient, kappa_bin (1/m)
    for ( size_t ib=0; ib<j_bin.size(); ++ib ) {
        // Write to file
        specfile << setw(20) << j_bin[ib]
                 << setw(20) << kappa_bin[ib]
                 << endl;
        // move on to next bin
    }

    specfile.close();

    return;
}

ApparatusFunction::ApparatusFunction( string name, double nu_sample )
 : name( name ), nu_sample( nu_sample ), f_scale( 1.0 ) {}

ApparatusFunction::~ApparatusFunction() {}

void ApparatusFunction::initialise()
{
    /* FIXME: this might not be necessary if we rescale in apply_apparatus_function */
#   if 0
    ofstream ofile;
    ofile.open("profile.txt");
    ofile << setprecision(12) << scientific << showpoint;
    cout << "initialising the apparatus function" << endl;
    // Numerically integrate the eval function to get the unity factor
    double nu_star = nu2lambda(500.0);
    double gamma_star_Hz = gamma_star / ( nu2lambda(nu_star) * 10.0 ) * nu_star;
    double nu_min = -100 * gamma_star_Hz + nu_star;
    double nu_max = 100 * gamma_star_Hz + nu_star;
    double nu = nu_min;
    double dnu = ( nu_max - nu_min ) / 1.0e5;
    double integral = 0.0;
    // FIXME: replace with simpsons Rule
    while ( nu <= nu_max ) {
        integral += this->eval(nu,nu_star-nu) * dnu;
        ofile << nu2lambda(nu) << "\t" << this->eval(nu,nu_star-nu) << endl;
        nu += dnu;
    }

    cout << "integral = " << integral << endl;

    f_scale = 1.0 / integral;
#   endif
}

Voigt::Voigt(double gamma_L, double gamma_G, double nu_sample)
 : ApparatusFunction("Voigt",nu_sample), gamma_L( gamma_L), gamma_G( gamma_G )
{
    gamma_V = calculate_Voigt_width( gamma_L, gamma_G );
    gamma_star = gamma_V;
}

Voigt::~Voigt() {}

double Voigt::eval( double nu, double delta_nu )
{
    double lambda_Ang = nu2lambda(nu) * 10.0;
    double gamma_V_Hz = gamma_V / lambda_Ang * nu;
    double gamma_L_Hz = gamma_L / lambda_Ang * nu;
    double gamma_G_Hz = gamma_G / lambda_Ang * nu;
    double b_V = eval_Voigt_profile( delta_nu, gamma_V_Hz, gamma_L_Hz, gamma_G_Hz );
    return b_V * f_scale;
}

SQRT_Voigt::SQRT_Voigt(double gamma_L, double gamma_G, double nu_sample)
 : ApparatusFunction("SQRT_Voigt",nu_sample), gamma_L( gamma_L), gamma_G( gamma_G )
{
    gamma_V = calculate_Voigt_width( gamma_L, gamma_G );
    /* FIXME: need a more accurate representative width */
    gamma_star = gamma_V;
}

SQRT_Voigt::~SQRT_Voigt() {}

double SQRT_Voigt::eval( double nu, double delta_nu )
{
    double lambda_Ang = nu2lambda(nu) * 10.0;
    double gamma_V_Hz = gamma_V / lambda_Ang * nu;
    double gamma_L_Hz = gamma_L / lambda_Ang * nu;
    double gamma_G_Hz = gamma_G / lambda_Ang * nu;
    double b_V = eval_Voigt_profile( delta_nu, gamma_V_Hz, gamma_L_Hz, gamma_G_Hz );
    return sqrt(b_V) * f_scale;
}

Gaussian_Lorentzian_hybrid::Gaussian_Lorentzian_hybrid(double gamma_L, double gamma_G, double f_pow, double nu_sample)
 : ApparatusFunction("Gaussian_Lorentzian_hybrid",nu_sample), gamma_L( gamma_L), gamma_G( gamma_G ), f_pow( f_pow )
{
    /* FIXME: need a more accurate representative width */
    gamma_star = 0.5 * ( gamma_G + gamma_L );
}

Gaussian_Lorentzian_hybrid::~Gaussian_Lorentzian_hybrid() {}

double Gaussian_Lorentzian_hybrid::eval( double nu, double delta_nu )
{
    cout << "Gaussian_Lorentzian_hybrid::eval()" << endl
	 << "Not implemented yet!" << endl;
    exit( NOT_IMPLEMENTED_ERROR );
}

Gaussian_profile::Gaussian_profile(double gamma_G, double nu_sample)
 : ApparatusFunction("Gaussian_profile",nu_sample), gamma_G( gamma_G )
{
    gamma_star = gamma_G;
}

Gaussian_profile::~Gaussian_profile() {}

double Gaussian_profile::eval( double nu, double delta_nu )
{
    double lambda_Ang = nu2lambda(nu) * 10.0;
    double gamma_G_Hz = gamma_G / lambda_Ang * nu;
    double b_V = eval_Gaussian_profile( delta_nu, gamma_G_Hz );
    return sqrt(b_V);
}

/* ------------ SpectralIntensity class ------------ */

SpectralIntensity::SpectralIntensity() {}

SpectralIntensity::SpectralIntensity( RadiationSpectralModel * rsm )
 : SpectralContainer( rsm )
{
    I_nu.resize( nu.size(), 0.0 );
    I_int.resize( nu.size(), 0.0 );
}

SpectralIntensity::SpectralIntensity( RadiationSpectralModel * rsm, double T )
 : SpectralContainer( rsm )
{
    I_nu.resize( nu.size(), 0.0 );
    I_int.resize( nu.size(), 0.0 );

    for ( size_t inu=0; inu<nu.size(); ++inu ) {
    	I_nu[inu] = planck_intensity( nu[inu], T );
    }
}

SpectralIntensity::SpectralIntensity(SpectralIntensity &S )
 : SpectralContainer( S ), I_nu( S.I_nu ), I_int( S.I_int )
{}

SpectralIntensity::~SpectralIntensity()
{
    I_nu.resize(0);
    I_int.resize(0);
}

double SpectralIntensity::write_to_file( string filename, int spectral_units )
{
    string Y1_label = "Spectral intensity, I_lambda (W/m**2-sr-m)";
    if ( spectral_units==WAVENUMBER )
        Y1_label = "Spectral intensity, I_eta (W/m**2-sr-1/cm)";
    else if ( spectral_units==FREQUENCY )
        Y1_label = "Spectral intensity, I_nu (W/m**2-sr-Hz)";
    else if ( spectral_units==ENERGY )
        Y1_label = "Spectral intensity, I_epsilon (W/m**2-sr-eV)";
    string Y1_int_label = "Integrated intensity, I (W/m**2-sr)";
    
    return write_data_to_file( filename, spectral_units, I_nu, Y1_label, Y1_int_label );
}

void SpectralIntensity::apply_apparatus_function( ApparatusFunction * A  )
{
    // Quick exit if the representative width is too small
    if ( A->gamma_star < 1.0e-10 ) return;

    // Initialise the Apparatus function
    A->initialise();

    // A vector to temporarily hold smeared data
    vector<double> nu_temp;
    int count = 0;
    for( size_t inu=0; inu<nu.size(); inu++) {
	count++;
        if ( count==A->nu_sample ) {
            nu_temp.push_back(nu[inu]);
            count = 0;
        }
    }
    vector<double> I_nu_temp( nu_temp.size() );

    int percentage=0;
    for( size_t inu=0; inu<nu_temp.size(); inu++) {
	double nu_val = nu_temp[inu];
	double lambda_ang = 10.0 * nu2lambda( nu_val );
	// convert HWHM's to Hz
	double gamma_star_Hz = A->gamma_star / lambda_ang * nu_val;
	double nu_lower = nu_val - double(nwidths) * gamma_star_Hz;
	double nu_upper = nu_val + double(nwidths) * gamma_star_Hz;
	int jnu_start = get_nu_index(nu,nu_lower,adaptive) + 1;
	int jnu_end = get_nu_index(nu,nu_upper,adaptive) + 1;

	if((double(inu)/double(nu_temp.size())*100.0+1.0)>double(percentage)) {
	    cout << "Smearing spectrum: | " << percentage << "% |, jnu_start = " << jnu_start << ", jnu_end = " << jnu_end << " \r" << flush;
	    percentage += 10;
        }

	// Apply convolution integral over this frequency range with trapezoidal method
	double I_nu_conv = 0.0;
	double AF_integral = 0.0;
	for ( int jnu=jnu_start+1; jnu<jnu_end; jnu++ ) {
            double dnu0 = nu[jnu-1] - nu_val;
            double dnu1 = nu[jnu] - nu_val;
            double f_nu0 = A->eval(nu_val, dnu0);
            double f_nu1 = A->eval(nu_val, dnu1);
            I_nu_conv += 0.5 * ( I_nu[jnu-1]*f_nu0 + I_nu[jnu]*f_nu1 ) * ( nu[jnu] - nu[jnu-1] );
            AF_integral += 0.5 * ( f_nu0 + f_nu1 ) * ( nu[jnu] - nu[jnu-1] );
	}

	// cout << "AF_integral = " << AF_integral << endl;

	// Make sure a zero value is not returned if nwidths is too small
	if ( jnu_start==(jnu_end-1) ) {
	    cout << "SpectralIntensity::apply_apparatus_function()" << endl
	         << "WARNING: nwidths is too small!" << endl;
	    I_nu_conv = I_nu[inu];
	}

	// Save result (and rescale to ensure the integral remains the same)
	I_nu_temp[inu] = I_nu_conv / AF_integral;
    }

    cout << endl;

    // Loop again to overwrite old I_nu values in S
    nu.resize(nu_temp.size());
    I_nu.resize(nu_temp.size());
    I_int.resize(nu_temp.size());
    for( size_t inu=0; inu<nu_temp.size(); inu++) {
	nu[inu] = nu_temp[inu];
    	I_nu[inu] = I_nu_temp[inu];
    	// recompute cumulative intensity
    	if ( inu==0.0 )
    	    I_int[inu] = 0.0;
    	else
    	    I_int[inu] = 0.5 * ( I_nu[inu] + I_nu[inu-1] ) * ( nu[inu] - nu[inu-1] );
    }

    return;
}

void SpectralIntensity::reverse_data_order()
{
    reverse(I_nu.begin(), I_nu.end());
    reverse(I_int.begin(), I_int.end());
    reverse(nu.begin(), nu.end());
    
    return;
}

void SpectralIntensity::reset_intensity_vectors()
{
    I_nu.resize( nu.size() );
    I_int.resize( nu.size() );
    // FIXME: find the c++ function that sets all values in a vector to zero
    for ( size_t inu=0; inu<nu.size(); ++inu ) {
        I_nu[inu] = 0.0;
        I_int[inu] = 0.0;
    }
    
    return;
}

double SpectralIntensity::integrate_intensity_spectra( double lambda_min, double lambda_max )
{
    // 1. Find spectral indice range to integrate over
    int inu_start = 0;
    if ( lambda_max > 0.0 ) inu_start = get_nu_index(nu, lambda2nu(lambda_max), adaptive) + 1;
    
    int inu_end = nu.size() - 1;
    if ( lambda_min > 0.0 ) inu_end = get_nu_index(nu, lambda2nu(lambda_min), adaptive) + 1;

    double I_total = 0.0;
    for( int inu=inu_start+1; inu<=inu_end; ++inu ) {
     	I_total += 0.5 * ( I_nu[inu] + I_nu[inu-1] ) * ( nu[inu] - nu[inu-1] );
     	I_int[inu] = I_total;
    }

    return I_total;
}

#define QUADRATIC_EQUATION(a,b,c) ( -b + sqrt( b*b - 4.0 * a * c ) ) / ( 2.0 * a )

double SpectralIntensity::random_frequency( double R )
{
    double I_interval = R * I_int.back();
    // cout << "I_int.back() = " << I_int.back() << ", I_interval = " << I_interval << endl;
    int inu_b = 0, inu_u = int( I_int.size() - 1 ), inu_mp = 0;
    double f_b = 0.0, f_mp;
    while ( abs( inu_u - inu_b ) > 1 ) {
        f_b = I_int[inu_b] - I_interval;
        inu_mp = int( ( inu_u + inu_b ) / 2.0 );
        f_mp = I_int[inu_mp] - I_interval;
        ( f_b * f_mp > 0.0 ) ? ( inu_b = inu_mp ) : ( inu_u = inu_mp );
    }

    // now need to solve for the actual frequency
    double m = ( I_nu[inu_u] - I_nu[inu_b] ) / ( nu[inu_u] - nu[inu_b] );
    double nu_interval = QUADRATIC_EQUATION( (0.5*m), (-m*nu[inu_b]+I_nu[inu_b]), (0.5*m*nu[inu_b]*nu[inu_b] - I_nu[inu_b]*nu[inu_b] - ( I_interval - I_int[inu_b])) );

    // cout << "SpectralIntensity::random_frequency()" << endl;
    // cout << "inu_b = " << inu_b << ", inu_u = " << inu_u << endl;
    // cout << "I_interval = " << I_interval << ", I_int[inu_b] = " << I_int[inu_b] << ", I_int[inu_u] = " << I_int[inu_u] << endl;
    // cout << "nu_interval = " << nu_interval << ", nu[inu_b] = " << nu[inu_b] << ", nu[inu_u] = " << nu[inu_u] << endl;
    // cout << "I_nu[inu_b] = " << I_nu[inu_b] << ", I_nu[inu_u] = " << I_nu[inu_u] << endl;

    // if ( nu_interval > nu[inu_u] ) {
    //     cout << "SpectralIntensity::random_frequency()" << endl
    //          << "nu_interval = " << nu_interval << " > nu[inu_u] = " << nu[inu_u] << endl;
    //     exit( FAILURE );
    // }
    // else if ( nu_interval < nu[inu_b] ) {
    //     cout << "SpectralIntensity::random_frequency()" << endl
    //          << "nu_interval = " << nu_interval << " < nu[inu_b] = " << nu[inu_b] << endl;
    //     exit( FAILURE );
    // }
    return nu_interval;
}

#undef QUADRATIC_EQUATION

int SpectralIntensity::random_frequency_interval( double R )
{
    double I_interval = R * I_int.back();
    // cout << "I_interval = " << I_interval << ", R = " << R << ", I_int.back() = " << I_int.back() << endl;
    int inu_b = 0, inu_u = int( I_int.size() - 1 ), inu_mp = 0;
    double f_b = 0.0, f_mp;
    while ( abs( inu_u - inu_b ) > 1 ) {
        f_b = I_int[inu_b] - I_interval;
        inu_mp = int( ( inu_u + inu_b ) / 2.0 );
        f_mp = I_int[inu_mp] - I_interval;
        ( f_b * f_mp > 0.0 ) ? ( inu_b = inu_mp ) : ( inu_u = inu_mp );
    }
    f_b = I_int[inu_b] - I_interval;
    double f_u = I_int[inu_u] - I_interval;

    return fabs(f_b) < fabs(f_u) ? inu_b : inu_u;
}

/* ------------ BinnedSpectralIntensity class ------------ */

BinnedSpectralIntensity::BinnedSpectralIntensity( size_t N_bins )
{
    I_bin.resize(N_bins);
}

BinnedSpectralIntensity::BinnedSpectralIntensity( SpectralIntensity * I, vector<SpectralBin*> & B )
{
    int nnu = I->nu.size();

    // 1. Calculate the I_bin vector
    for ( int iB=0; iB<int(B.size()); ++iB ) {
	double I_int = 0.0;
	for ( int j=0; j<int(B[iB]->inu.size()); ++j) {
	    int inu = B[iB]->inu[j];
	    // NOTE: we are taking care here to remain consistent with a trapezoidal discretisation of the spectra
	    //       where the trapezoid heights are taken as the average of the two bounding points
	    if      ( inu==0 )     I_int += I->I_nu[inu] * 0.5 * fabs( I->nu[inu+1] - I->nu[inu] );
	    else if ( inu==nnu-1 ) I_int += I->I_nu[inu] * 0.5 * fabs( I->nu[inu] - I->nu[inu-1] );
	    else                   I_int += I->I_nu[inu] * ( 0.5 * fabs( I->nu[inu] - I->nu[inu-1] ) + 0.5 * fabs( I->nu[inu+1] - I->nu[inu] ) );
	}
	I_bin.push_back( I_int );
    }

}

BinnedSpectralIntensity::BinnedSpectralIntensity( RadiationSpectralModel * rsm, double T, vector<SpectralBin*> & B )
{
    SpectralIntensity S( rsm, T );

    int nnu = S.nu.size();

    // 1. Calculate the I_bin vector
    for ( int iB=0; iB<int(B.size()); ++iB ) {
	double I_int = 0.0;
	for ( int j=0; j<int(B[iB]->inu.size()); ++j) {
	    int inu = B[iB]->inu[j];
	    // NOTE: we are taking care here to remain consistent with a trapezoidal discretisation of the spectra
	    //       where the trapezoid heights are taken as the average of the two bounding points
	    if      ( inu==0 )     I_int += S.I_nu[inu] * 0.5 * fabs( S.nu[inu+1] - S.nu[inu] );
	    else if ( inu==nnu-1 ) I_int += S.I_nu[inu] * 0.5 * fabs( S.nu[inu] - S.nu[inu-1] );
	    else                   I_int += S.I_nu[inu] * ( 0.5 * fabs( S.nu[inu] - S.nu[inu-1] ) + 0.5 * fabs( S.nu[inu+1] - S.nu[inu] ) );
	}
	I_bin.push_back( I_int );
    }
}

BinnedSpectralIntensity::~BinnedSpectralIntensity()
{}

double BinnedSpectralIntensity::sum_intensity()
{
    double I_total = 0.0;
    for ( size_t ib=0; ib<I_bin.size(); ++ib )
	I_total += I_bin[ib];
    return I_total;
}

/* ------------ SpectralFlux class ------------ */

SpectralFlux::SpectralFlux() {}

SpectralFlux::SpectralFlux( RadiationSpectralModel * rsm )
 : SpectralContainer( rsm )
{
    q_nu.resize( nu.size(), 0.0 );
    q_int.resize( nu.size(), 0.0 );
}

SpectralFlux::SpectralFlux( RadiationSpectralModel * rsm, double T )
 : SpectralContainer( rsm )
{
    q_nu.resize( nu.size() );

    for ( size_t inu=0; inu<nu.size(); ++inu ) {
    	q_nu[inu] = 2.0 * M_PI * planck_intensity( nu[inu], T ) * E_3( 0.0 );
    }
}

SpectralFlux::~SpectralFlux()
{
    q_nu.resize(0);
}

void SpectralFlux::clear_data()
{
    nu.clear();
    q_nu.clear();
    q_int.clear();
}

double SpectralFlux::write_to_file( string filename, int spectral_units )
{
    string Y1_label = "Spectral flux, q_lambda (W/m**2-m)";
    if ( spectral_units==WAVENUMBER )
        Y1_label = "Spectral flux, q_eta (W/m**2-1/cm)";
    else if ( spectral_units==FREQUENCY )
        Y1_label = "Spectral flux, q_nu (W/m**2-Hz)";
    else if ( spectral_units==ENERGY )
        Y1_label = "Spectral flux, q_epsilon (W/m**2-eV)";
    string Y1_int_label = "Integrated flux, q (W/m**2)";
    
    return write_data_to_file( filename, spectral_units, q_nu, Y1_label, Y1_int_label );
}

void SpectralFlux::read_from_file( string fname, int inu_start, int inu_end )
{
    // make sure the vectors are clear
    this->clear_data();

    string line;
    ifstream specfile (fname.c_str());
    if (specfile.is_open()) {
        int spectral_units = -1;
        getline (specfile,line); // should be the filename
        getline (specfile,line); // should be spectral column descriptor
        if ( line.compare(0,11,"# Column 1:")==0 ) {
            if ( line.compare(12,14,"Frequency (Hz)")==0 ) {
                spectral_units = FREQUENCY;
            }
            else if ( line.compare(12,15,"Wavelength (nm)")==0 ) {
                spectral_units = WAVELENGTH;
            }
            else {
                cout << "SpectralFlux::read_from_file()" << endl
                     << "Only frequency and wavelength units are currently supported" << endl;
                exit(FAILURE);
            }

        }
        getline (specfile,line); // should be emission coefficient column descriptor
        getline (specfile,line); // should be abs coefficient column descriptor
        getline (specfile,line); // should be int emission coefficient column descriptor
        int count = 0;
        while ( specfile.good() ) {
            double wav, q_wav, _q_int;
            specfile >> wav >> q_wav >> _q_int;
            // cout << "nu = " << _nu << ", j_nu = " << _j_nu << ", kappa_nu = " << _kappa_nu << ", j_int = " << _j_int << endl;
            if ( count >= inu_start && ( count <= inu_end || inu_end < 0 ) ) {
                double _nu = 0.0, _q_nu = 0.0;
                if ( spectral_units==FREQUENCY ) {
                    _nu = wav;
                    _q_nu = q_wav;
                }
                else if ( spectral_units==WAVELENGTH ) {
                    _nu = nu2lambda(wav);
                    _q_nu = q_wav * wav * 1.0e-9 / _nu;
                }
                nu.push_back(_nu);
                q_nu.push_back(_q_nu);
                // not storing q_int, calculating later
                // q_int.push_back(_q_int);
            }
            count++;
        }
        specfile.close();
    }
    else {
        cout << "SpectralFlux::read_from_file()" << endl
             << "Unable to open file with name: " << fname << endl;
        exit(BAD_INPUT_ERROR);
    }

    if ( inu_end<0 ) {
        // remove the last entry which is a double-up
        nu.erase(nu.end()-1);
        q_nu.erase(q_nu.end()-1);
    }

    // Reverse the order
    reverse(nu.begin(),nu.end());
    reverse(q_nu.begin(),q_nu.end());

    cout << "Read " << nu.size() << " spectral points from file: " << fname << endl;
    q_int.resize(nu.size());
    double q_total = this->integrate_flux_spectra();
    cout << "Integral is " << q_total << endl;
}

void SpectralFlux::apply_apparatus_function( ApparatusFunction * A  )
{
    // Quick exit if the representative width is too small
    if ( A->gamma_star < 1.0e-10 ) return;

    // Initialise the Apparatus function
    A->initialise();

    // A vector to temporarily hold smeared data
    vector<double> nu_temp;
    int count = 0;
    for( size_t inu=0; inu<nu.size(); inu++) {
        count++;
        if ( count==A->nu_sample ) {
            nu_temp.push_back(nu[inu]);
            count = 0;
        }
    }
    vector<double> q_nu_temp( nu_temp.size() );

    int percentage=0;
    for( size_t inu=0; inu<nu_temp.size(); inu++) {
        double nu_val = nu_temp[inu];
        double lambda_ang = 10.0 * nu2lambda( nu_val );
        // convert HWHM's to Hz
        double gamma_star_Hz = A->gamma_star / lambda_ang * nu_val;
        double nu_lower = nu_val - double(nwidths) * gamma_star_Hz;
        double nu_upper = nu_val + double(nwidths) * gamma_star_Hz;
        int jnu_start = get_nu_index(nu,nu_lower,adaptive) + 1;
        int jnu_end = get_nu_index(nu,nu_upper,adaptive) + 1;

        if((double(inu)/double(nu_temp.size())*100.0+1.0)>double(percentage)) {
            cout << "Smearing spectrum: | " << percentage << "% |, jnu_start = " << jnu_start << ", jnu_end = " << jnu_end << " \r" << flush;
            percentage += 10;
        }

        // Apply convolution integral over this frequency range with trapezoidal method
        double q_nu_conv = 0.0;
        double AF_integral = 0.0;
        for ( int jnu=jnu_start+1; jnu<jnu_end; jnu++ ) {
            double dnu0 = nu[jnu-1] - nu_val;
            double dnu1 = nu[jnu] - nu_val;
            double f_nu0 = A->eval(nu_val, dnu0);
            double f_nu1 = A->eval(nu_val, dnu1);
            q_nu_conv += 0.5 * ( q_nu[jnu-1]*f_nu0 + q_nu[jnu]*f_nu1 ) * ( nu[jnu] - nu[jnu-1] );
            AF_integral += 0.5 * ( f_nu0 + f_nu1 ) * ( nu[jnu] - nu[jnu-1] );
        }

        // cout << "AF_integral = " << AF_integral << endl;

        // Make sure a zero value is not returned if nwidths is too small
        if ( jnu_start==(jnu_end-1) ) {
            cout << "SpectralIntensity::apply_apparatus_function()" << endl
                 << "WARNING: nwidths is too small!" << endl;
            q_nu_conv = q_nu[inu];
        }

        // Save result (and rescale to ensure the integral remains the same)
        q_nu_temp[inu] = q_nu_conv / AF_integral;
    }

    cout << endl;

    // Loop again to overwrite old q_nu values in S
    nu.resize(nu_temp.size());
    q_nu.resize(nu_temp.size());
    q_int.resize(nu_temp.size());
    for( size_t inu=0; inu<nu_temp.size(); inu++) {
        nu[inu] = nu_temp[inu];
        q_nu[inu] = q_nu_temp[inu];
        // recompute cumulative intensity
        if ( inu==0.0 )
            q_int[inu] = 0.0;
        else
            q_int[inu] = 0.5 * ( q_nu[inu] + q_nu[inu-1] ) * ( nu[inu] - nu[inu-1] );
    }

    return;
}

double SpectralFlux::integrate_flux_spectra( double lambda_min, double lambda_max )
{
    // 1. Find spectral indice range to integrate over
    int inu_start = 0;
    if ( lambda_max > 0.0 ) inu_start = get_nu_index(nu, lambda2nu(lambda_max), adaptive) + 1;

    int inu_end = nu.size() - 1;
    if ( lambda_min > 0.0 ) inu_end = get_nu_index(nu, lambda2nu(lambda_min), adaptive) + 1;

    double q_total = 0.0;
    for( int inu=inu_start+1; inu<=inu_end; ++inu ) {
        q_total += 0.5 * ( q_nu[inu] + q_nu[inu-1] ) * ( nu[inu] - nu[inu-1] );
        q_int[inu] = q_total;
    }

    return q_total;
}

/* ------------ BinnedSpectralFlux class ------------ */

BinnedSpectralFlux::BinnedSpectralFlux( size_t N_bins )
{
    q_bin.resize(N_bins);
}

BinnedSpectralFlux::BinnedSpectralFlux( SpectralFlux * F, vector<SpectralBin*> & B )
{
    double nnu = F->nu.size();

    // 1. Calculate the q_bin vector
    for ( int iB=0; iB<int(B.size()); ++iB ) {
	double q_int = 0.0;
	for ( int j=0; j<int(B[iB]->inu.size()); ++j) {
	    int inu = B[iB]->inu[j];
	    // NOTE: we are taking care here to remain consistent with a trapezoidal discretisation of the spectra
	    //       where the trapezoid heights are taken as the average of the two bounding points
	    if      ( inu==0 )     q_int += F->q_nu[inu] * 0.5 * fabs( F->nu[inu+1] - F->nu[inu] );
	    else if ( inu==nnu-1 ) q_int += F->q_nu[inu] * 0.5 * fabs( F->nu[inu] - F->nu[inu-1] );
	    else                   q_int += F->q_nu[inu] * ( 0.5 * fabs( F->nu[inu] - F->nu[inu-1] ) + 0.5 * fabs( F->nu[inu+1] - F->nu[inu] ) );
	}
	q_bin.push_back( q_int );
    }

}

BinnedSpectralFlux::BinnedSpectralFlux( RadiationSpectralModel * rsm, double T, vector<SpectralBin*> & B )
{
    SpectralFlux F( rsm, T );

    double nnu = F.nu.size();

    // 1. Calculate the q_bin vector
    for ( int iB=0; iB<int(B.size()); ++iB ) {
	double q_int = 0.0;
	for ( int j=0; j<int(B[iB]->inu.size()); ++j) {
	    int inu = B[iB]->inu[j];
	    // NOTE: we are taking care here to remain consistent with a trapezoidal discretisation of the spectra
	    //       where the trapezoid heights are taken as the average of the two bounding points
	    if      ( inu==0 )     q_int += F.q_nu[inu] * 0.5 * fabs( F.nu[inu+1] - F.nu[inu] );
	    else if ( inu==nnu-1 ) q_int += F.q_nu[inu] * 0.5 * fabs( F.nu[inu] - F.nu[inu-1] );
	    else                   q_int += F.q_nu[inu] * ( 0.5 * fabs( F.nu[inu] - F.nu[inu-1] ) + 0.5 * fabs( F.nu[inu+1] - F.nu[inu] ) );
	}
	q_bin.push_back( q_int );
    }
}

BinnedSpectralFlux::~BinnedSpectralFlux()
{}

double BinnedSpectralFlux::sum_flux()
{
    double q_total = 0.0;
    for ( size_t ib=0; ib<q_bin.size(); ++ib )
	q_total += q_bin[ib];
    return q_total;
}

/* ------------ IntensityProfile class ------------ */

IntensityProfile::IntensityProfile() {}

IntensityProfile::~IntensityProfile() {}

void
IntensityProfile::add_new_point( double x, double I )
{
    x_vec.push_back( x );
    I_vec.push_back( I );

    return;
}

void
IntensityProfile::spatially_smear( double dx_smear )
{
    if ( dx_smear < 0.0 ) return;
    
    if ( x_vec.size()<2 ) {
    	cout << "IntensityProfile::spatially_smear()" << endl
    	     << "Cannot spatially smear as there are only " << x_vec.size()
    	     << " points present.\nBailing out!" << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    double dx = ( x_vec[1] - x_vec[0] );
    int np_smear = int( dx_smear / dx );
    if ( np_smear < 1 ) {
    	cout << "IntensityProfile::spatially_smear()" << endl
    	     << "dx_smear = " << dx_smear << " need to be >> dx = " << dx << endl
    	     << "Returning unsmeared data." << endl;
    	     return;
    }
    int nps = int(x_vec.size());
    vector<double> x_tmp;
    vector<double> I_tmp;
    
    for ( int ip=0; ip < nps+np_smear; ip++ ) {
    	double s = x_vec[0] + dx * ( double( ip - np_smear) );
    	x_tmp.push_back( s );
    	I_tmp.push_back( 0.0 );
    	for ( int jp=ip-np_smear; jp<ip; ++jp ) {
    	    double I = 0.0;
    	    if ( jp < 0 ) I = 0.0;
    	    else if ( jp >= nps ) I = I_vec[nps-1];
	    else I = I_vec[jp];
	    I_tmp.back() += I;
	}
	// Divide by number of averaging points
	I_tmp.back() /= double(np_smear);
    }
    
    // Replace the original data with the smeared data
    x_vec = x_tmp;
    I_vec = I_tmp;
    
    return;
}

void
IntensityProfile::spatially_smear_for_varying_dx( double dx_smear )
{
    if ( dx_smear < 0.0 ) return;
    
    if ( x_vec.size()<2 ) {
    	cout << "IntensityProfile::spatially_smear_for_varying_dx()" << endl
    	     << "Cannot spatially smear as there are only " << x_vec.size()
    	     << " points present.\nBailing out!" << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    int nps = int(x_vec.size());
    vector<double> x_tmp;
    vector<double> I_tmp;
    
    for ( int ip=0; ip < nps; ip++ ) {
    	double s = x_vec[ip] - dx_smear;
    	x_tmp.push_back( s );
    	I_tmp.push_back( 0.0 );
    	double dx_weighting = 0.0;
    	for ( int jp=0; jp < nps; jp++ ) {
    	    if ( x_vec[jp] >= s && x_vec[jp] <= s + dx_smear ) {
    	    	double dx = 0.0;
    	    	if ( jp==0 ) dx = fabs(x_vec[jp+1] - x_vec[jp]);
    	    	else dx = fabs(x_vec[jp] - x_vec[jp-1]);
    	    	I_tmp.back() += I_vec[jp] * dx;
    	    	dx_weighting += dx;
    	    }
	}
	// Divide by number of averaging weighting
	I_tmp.back() /= dx_weighting;
    }
    
    double last_dx = x_vec[nps-1] - x_vec[nps-2];
    double s = x_vec[nps-1] - dx_smear + last_dx;
    while ( s <= x_vec[nps-1] ) {
     	x_tmp.push_back( s );
    	I_tmp.push_back( 0.0 );
    	double dx_weighting = 0.0;
    	for ( int jp=0; jp < nps; jp++ ) {
    	    if ( x_vec[jp] >= s && x_vec[jp] <= s + dx_smear ) {
    	    	double dx = 0.0;
    	    	if ( jp==0 ) dx = fabs(x_vec[jp+1] - x_vec[jp]);
    	    	else dx = fabs(x_vec[jp] - x_vec[jp-1]);
    	    	I_tmp.back() += I_vec[jp] * dx;
    	    	dx_weighting += dx;
    	    }
	}
	// Divide by number of averaging weighting
	I_tmp.back() /= dx_weighting;
	// increment s
	s += last_dx;
    }
    
    // Replace the original data with the smeared data
    x_vec = x_tmp;
    I_vec = I_tmp;
    
    return;
}

void
IntensityProfile::write_to_file( string fname )
{
    /* 1. Setup the output file. */
    ofstream ofile;
    ofile.open(fname.c_str());
    if( ofile.fail() ) {
	cout << "Error opening file: " << fname << endl;
	cout << "Bailing Out!\n";
	exit(FILE_ERROR);
    }
    
    ofile << setprecision(12) << scientific << showpoint;

    ofile << "# " << fname << endl
    	  << "# Column 1: x (m)" << endl
          << "# Column 2: Intensity, I (W/m2-sr)" << endl;
    
    /* 2. Write data to file. */
    for ( int ip=0; ip < int(x_vec.size()); ++ip ) {
	// Write to file
        ofile << setw(20) << x_vec[ip] << setw(20) << I_vec[ip] << endl;
    }
    
    ofile.close();
}

/* ------------ SpectralField class ------------ */

SpectralField::SpectralField() {}

SpectralField::~SpectralField() {}

void
SpectralField::add_new_intensity_spectra( double x, SpectralIntensity * S )
{
    x_vec.push_back( x );
    S_vec.push_back( new SpectralIntensity(*S) );
    
    cout << "S_vec.back()->nu.size() = " << S_vec.back()->nu.size() << endl;
    cout << "S->nu.size() = " << S->nu.size() << endl;

    return;
}

IntensityProfile
SpectralField::
extract_intensity_profile( double lambda_l, double lambda_u )
{
    if ( S_vec.size()==0 ) {
	cout << "SpectralField::extract_intensity_profile()" << endl
	     << "No spectra present, bailing out!" << endl;
	exit( BAD_INPUT_ERROR );
    }
    
    // 1. Find spectral indice range to integrate over
    int inu_start = get_nu_index(S_vec[0]->nu, lambda2nu(lambda_u), S_vec[0]->adaptive) + 1;
    int inu_end = get_nu_index(S_vec[0]->nu, lambda2nu(lambda_l), S_vec[0]->adaptive) + 1;
    
    // 2. Loop over LOS_points and spectrally integrate
    IntensityProfile IvX;
    for ( size_t ip=0; ip<x_vec.size(); ++ip ) {
    	double I_total = 0.0;
    	for( int inu=inu_start+1; inu<inu_end; ++inu )
    	    I_total += 0.5 * ( S_vec[ip]->I_nu[inu-1] + S_vec[ip]->I_nu[inu] ) * ( S_vec[ip]->nu[inu] - S_vec[ip]->nu[inu-1] );
    	IvX.add_new_point( x_vec[ip], I_total );
    }
    
    return IvX;
}

/* ------------- Helper functions -------------- */

double nu2lambda( double nu ) 
{
    double lambda_nm = RC_c_SI / nu * 1.0e9;
    
    return lambda_nm;
}

double lambda2nu( double lambda_nm ) 
{
    double nu = RC_c_SI / lambda_nm * 1.0e9;
    
    return nu;
}

double planck_intensity(const double nu, const double T)
{
    double B_nu = 0.0;
    
    // Impose a lower temperature limit for convenience in initialisations
    if ( T > 50.0 ) {
	double top = 2.0 * RC_h_SI * nu * nu * nu;
	double bottom = RC_c_SI * RC_c_SI * ( exp( RC_h_SI * nu / ( RC_k_SI * T ) ) - 1.0 );
	B_nu = top / bottom;
    }

    return B_nu;
}

int get_nu_index( vector<double> &nu, double nu_star, bool adaptive )
{
    // 0. Firstly check if nu is in range
    int nnu = int ( nu.size() );
    int inu;
    if ( nu_star < nu.front() ) inu=-1;
    else if ( nu_star > nu.back() ) inu=nnu-1;
    // 1. nu is in range, so find the appropriate index
    else if ( adaptive ){
        // The grid is non-uniform
        int inu_b = 0, inu_u = int( nu.size() - 1 ), inu_mp = 0;
        double f_b = 0.0, f_mp;
        while ( abs( inu_u - inu_b ) > 1 ) {
            f_b = nu[inu_b] - nu_star;
            inu_mp = int( ( inu_u + inu_b ) / 2.0 );
            f_mp = nu[inu_mp] - nu_star;
            ( f_b * f_mp > 0.0 ) ? ( inu_b = inu_mp ) : ( inu_u = inu_mp );
        }
        inu = inu_b;
    }
    else {
        // The grid is uniform in frequency space
	double dnu = ( nu.back() - nu.front() ) / double ( nnu - 1 );
	inu = int ( ( nu_star - nu.front() ) / dnu );
    }
    
    return inu;
}

