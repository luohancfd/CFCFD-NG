/** \file gas_data.cxx
 *  \ingroup gas
 *
 *  \author Rowan J Gollan
 *  \version 04-July-2008
 *  \version 10-May-2010: PJ: introduce copying of matrix elements.
 *  \version 19-July-2010: PJ: class-ified
 **/

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "../../util/source/useful.h"
#include "gas_data.hh"
#include "gas-model.hh"

using namespace std;

double get_matrix_element(const matrix &M, int i, int j)
{ 
    return M[i][j];
}

int copy_matrix_elements(const matrix &srcM, matrix &destM)
{ 
    for ( size_t i = 0; i < srcM.size(); ++i ) {
	for ( size_t j = 0; j < srcM[i].size(); ++j ) {
	    destM[i][j] = srcM[i][j];
	}
    }
    return SUCCESS;
}

int average_matrix_elements(const matrix &srcM0, const matrix &srcM1, matrix &destM)
{ 
    for ( size_t i = 0; i < srcM0.size(); ++i ) {
	for ( size_t j = 0; j < srcM0[i].size(); ++j ) {
	    destM[i][j] = 0.5 * (srcM0[i][j] + srcM1[i][j]);
	}
    }
    return SUCCESS;
}

Gas_data::Gas_data(Gas_model *gm)
{
    int nsp = gm->get_number_of_species();
    int nmodes = gm->get_number_of_modes();
    rho = 0.0;
    p = 0.0;
    p_e = 0.0;
    a = 0.0;
    mu = 0.0;
    quality = 1.0;
    massf.resize(nsp, 0.0);
    T.resize(nmodes, 0.0);
    e.resize(nmodes, 0.0);
    k.resize(nmodes, 0.0);
    D_AB.resize(nsp);
    for( size_t i = 0; i < D_AB.size(); ++i ) D_AB[i].resize(nsp, 0.0);
}

Gas_data::Gas_data(const Gas_data &Q)
{
    int nsp = (int) Q.massf.size();
    int nmodes = (int) Q.T.size();
    massf.resize(nsp, 0.0);
    T.resize(nmodes, 0.0);
    e.resize(nmodes, 0.0);
    k.resize(nmodes, 0.0);
    D_AB.resize(nsp);
    for( size_t i = 0; i < D_AB.size(); ++i ) D_AB[i].resize(nsp, 0.0);
    this->copy_values_from(Q);
}

Gas_data::Gas_data(int nsp, int nmodes)
{
    massf.resize(nsp, 0.0);
    T.resize(nmodes, 0.0);
    e.resize(nmodes, 0.0);
    k.resize(nmodes, 0.0);
    D_AB.resize(nsp);
    for( size_t i = 0; i < D_AB.size(); ++i ) D_AB[i].resize(nsp, 0.0);
}

Gas_data::~Gas_data() 
{
    for( size_t i = 0; i < D_AB.size(); ++i ) D_AB[i].resize(0);
    D_AB.resize(0);
    massf.resize(0);
    T.resize(0);
    e.resize(0);
    k.resize(0);
}

Gas_data &
Gas_data::operator=(const Gas_data &Q)
{
    if (this != &Q) { // Avoid aliasing
        int nsp = (int) Q.massf.size();
        int nmodes = (int) Q.T.size();
        massf.resize(nsp, 0.0);
        T.resize(nmodes, 0.0);
        e.resize(nmodes, 0.0);
        k.resize(nmodes, 0.0);
        D_AB.resize(nsp);
        for( size_t i = 0; i < D_AB.size(); ++i ) D_AB[i].resize(nsp, 0.0);
        this->copy_values_from(Q);
    }
    return *this;
}

void Gas_data::print_values(bool print_transport_data) const
{
    cout << showpoint;
    cout << setprecision(12);
    cout << "rho= " << rho << " p= " << p << " p_e= " << p_e << " a= " << a << endl;
    for( size_t i = 0; i < T.size(); ++i ) {
	cout << "T[" << i << "]= " << T[i] << " ";
    }
    cout << endl;
    for( size_t i = 0; i < e.size(); ++i ) {
	cout << "e[" << i << "]= " << e[i] << " ";
    }
    cout << endl;
    if ( print_transport_data==true ) {
	cout << "mu= " << mu << endl;
	for ( size_t i = 0; i < k.size(); ++i ) {
	    cout << "k[" << i << "]= " << k[i] << " ";
	}
	cout << endl;
	for( size_t i = 0; i < D_AB.size(); ++i ) {
	    for( size_t j = 0; j < D_AB[i].size(); ++j ) {
		cout << "D_AB[" << i << "][" << j << "]= " << D_AB[i][j] << " ";
	    }
	    cout << endl;
	}
    }
    for( size_t isp = 0; isp < massf.size(); ++isp ) {
	cout << "massf[" << isp << "]= " << massf[isp] << " ";
    }
    cout << endl;
} // end Gas_data::print_values()

void Gas_data::write_values(string fname, bool print_transport_data) const
{
    /* 1. Setup the output file. */
    ofstream gasfile;
    gasfile.open(fname.c_str());
    if( gasfile.fail() ) {
        cout << "Error opening file: " << fname << endl;
        cout << "Bailing Out!\n";
        exit(FILE_ERROR);
    }

    gasfile << showpoint;
    gasfile << setprecision(12);
    gasfile << "rho= " << rho << " p= " << p << " p_e= " << p_e << " a= " << a << endl;
    for( size_t i = 0; i < T.size(); ++i ) {
        gasfile << "T[" << i << "]= " << T[i] << " ";
    }
    gasfile << endl;
    for( size_t i = 0; i < e.size(); ++i ) {
        gasfile << "e[" << i << "]= " << e[i] << " ";
    }
    gasfile << endl;
    if ( print_transport_data==true ) {
        gasfile << "mu= " << mu << endl;
        for ( size_t i = 0; i < k.size(); ++i ) {
            gasfile << "k[" << i << "]= " << k[i] << " ";
        }
        gasfile << endl;
        for( size_t i = 0; i < D_AB.size(); ++i ) {
            for( size_t j = 0; j < D_AB[i].size(); ++j ) {
                gasfile << "D_AB[" << i << "][" << j << "]= " << D_AB[i][j] << " ";
            }
            gasfile << endl;
        }
    }
    for( size_t isp = 0; isp < massf.size(); ++isp ) {
        gasfile << "massf[" << isp << "]= " << massf[isp] << " ";
    }
    gasfile << endl;

    gasfile.close();
} // end Gas_data::write_values()

void Gas_data::copy_values_from(const Gas_data &src)
{
    // 1. doubles
    rho = src.rho;
    p = src.p;
    p_e = src.p_e;
    a = src.a;
    mu = src.mu;
    quality = src.quality;
    // 2. vectors
    // 2a. loop over thermal modes
    for ( size_t itm=0; itm<src.e.size(); ++itm ) {
    	e[itm] = src.e[itm];
    	T[itm] = src.T[itm];
    	k[itm] = src.k[itm];
    }
    // 2b. loop over mass-fractions
    for ( size_t isp=0; isp<src.massf.size(); ++isp ) {
    	massf[isp] = src.massf[isp];
    }
    // 3. matrices
    copy_matrix_elements(src.D_AB, D_AB);
} // end Gas_data::copy_values_from()

void Gas_data::average_values_from(const Gas_data &src0, const Gas_data &src1, 
				   bool with_diff_coeff)
{
    // 1. doubles
    rho = 0.5 * (src0.rho + src1.rho);
    p = 0.5 * (src0.p + src1.p);
    p_e = 0.5 * (src0.p_e + src1.p_e);
    a = 0.5 * (src0.a + src1.a);
    mu = 0.5 * (src0.mu + src1.mu);
    // 2. vectors
    // 2a. loop over thermal modes
    for ( size_t itm=0; itm<src0.e.size(); ++itm ) {
    	e[itm] = 0.5 * (src0.e[itm] + src1.e[itm]);
    	T[itm] = 0.5 * (src0.T[itm] + src1.T[itm]);
    	k[itm] = 0.5 * (src0.k[itm] + src1.k[itm]);
    }
    // 2b. loop over mass-fractions
    for ( size_t isp=0; isp<src0.massf.size(); ++isp ) {
    	massf[isp] = 0.5 * (src0.massf[isp] + src1.massf[isp]);
    }
    // 3. matrices
    if ( with_diff_coeff ) {
	average_matrix_elements(src0.D_AB, src1.D_AB, D_AB);
    }
} // end Gas_data::average_values_from()

int Gas_data::accumulate_values_from(const Gas_data &src, double alpha)
{
    if ( alpha == 0.0 ) {
	/* This special case avoids the src data contaminating 
	 * the result with NaNs. */
	this->copy_values_from(src);
    } else {
        /* We presume that both sets of data are valid. */
	rho = alpha * rho + src.rho;
	p = alpha * p + src.p;
	p_e = alpha * p_e + src.p_e;
	a = alpha * a + src.a;
	for ( size_t i = 0; i < src.T.size(); ++i ) {
	    e[i] = alpha * e[i] + src.e[i];
	    T[i] = alpha * T[i] + src.T[i];
	    k[i] = alpha * k[i] + src.k[i];
	}
	mu = alpha * mu + src.mu;
	for ( size_t isp = 0; isp < src.massf.size(); ++isp ) {
	    massf[isp] = alpha * massf[isp] + src.massf[isp];
	}
	for ( size_t i = 0; i < src.D_AB.size(); ++i ) {
	    for ( size_t j = 0; j < src.D_AB[i].size(); ++j ) {
		D_AB[i][j] = alpha * D_AB[i][j] + src.D_AB[i][j];
	    }
	}
    }
    return SUCCESS;
} // end Gas_data::accumulate_values_from()

double* Gas_data::copy_values_to_buffer(double *buf) const
{
    *buf++ = rho;
    *buf++ = p;
    *buf++ = p_e;
    *buf++ = a;
    for ( size_t i = 0; i < e.size(); ++i ) *buf++ = e[i];
    for ( size_t i = 0; i < T.size(); ++i ) *buf++ = T[i];
    *buf++ = mu;
    for ( size_t i = 0; i < k.size(); ++i ) *buf++ = k[i];
    for ( size_t i = 0; i < D_AB.size(); ++i ) {
	for ( size_t j = 0; j < D_AB[i].size(); ++j ) {
	    *buf++ = D_AB[i][j];
	}
    }
    for ( size_t isp = 0; isp < massf.size(); ++isp ) *buf++ = massf[isp];
    return buf;
} // end Gas_data::copy_values_to_buffer()

double* Gas_data::copy_values_from_buffer(double *buf)
{
    rho = *buf++;
    p = *buf++;
    p_e = *buf++;
    a = *buf++;
    for ( size_t i = 0; i < e.size(); ++i ) e[i] = *buf++;
    for ( size_t i = 0; i < T.size(); ++i ) T[i] = *buf++;
    mu = *buf++;
    for ( size_t i = 0; i < k.size(); ++i ) k[i] = *buf++;
    for ( size_t i = 0; i < D_AB.size(); ++i ) {
	for ( size_t j = 0; j < D_AB[i].size(); ++j ) {
	    D_AB[i][j] = *buf++;
	}
    }
    for ( size_t isp = 0; isp < massf.size(); ++isp ) massf[isp] = *buf++;
    return buf;
} // end Gas_data::copy_values_from_buffer()

int Gas_data::check_values(bool print_message) const
{
    // Mostly borrowed from check_flow_data_for_cell() in app/eilmer3/source/cell.cxx
    double RHOMIN = 0.0;
    double TMIN = 1.0;
    int data_valid = 1;

    if ( !(finite(rho)) || rho < 1.01 * RHOMIN ) {
	if (print_message) cout << "Density invalid: " << rho << endl;
	data_valid = 0;
    }
    size_t nmodes = e.size();
    for ( size_t imode=0; imode<nmodes; ++imode ) {
    	if ( !(finite(T[imode])) || T[imode] < 1.01 * TMIN ) {
    	    if ( print_message ) 
		cout << "Temperature[" << imode << "] invalid: " << T[imode] << endl;
    	    data_valid = 0;
    	}
    	if ( !(finite(e[imode])) ) {
    	    if ( print_message ) 
		cout << "Energy[" << imode << "] invalid: " << e[imode] << endl;
    	    data_valid = 0;
    	}
    }
    if ( !(finite(p)) ) {
	if ( print_message ) cout << "Total pressure invalid: " << p << endl;
	data_valid = 0;
    }
    if ( !(finite(p_e)) ) {
	if ( print_message ) cout << "Electron pressure invalid: " << p_e << endl;
	data_valid = 0;
    }
    if ( !(finite(a)) ) {
	if ( print_message ) cout << "Sound speed invalid: " << a << endl;
	data_valid = 0;
    }
    size_t nsp = massf.size();
    double f_sum = 0.0;
    for ( size_t isp = 0; isp < nsp; ++isp ) f_sum += massf[isp];
    if ( f_sum < 0.99 || f_sum > 1.01 || !finite(f_sum) ) {
	if ( print_message ) cout << "Mass fraction sum bad: " << f_sum << endl;
	data_valid = 0;
    }
    if ( data_valid == 0 ) {
	this->print_values();
    }
    return data_valid;
} // end of Gas_data::check_values()


double mass_average(const vector<double> &massf, const vector<double> &vec)
{
    double val = 0.0;
    for( size_t isp = 0; isp < massf.size(); ++isp ) {
	val += massf[isp]*vec[isp];
    }
    return val;
}

double mass_average_inverse(const vector<double> &massf, const vector<double> &vec)
{
    double val = 0.0;
    for( size_t isp = 0; isp < massf.size(); ++isp ) {
        val += massf[isp]/vec[isp];
    }
    return val;
}

double mole_average(const vector<double> &molef, const vector<double> &vec)
{
    double val = 0.0;
    for( size_t isp = 0; isp < molef.size(); ++isp ) {
	val += molef[isp]*vec[isp];
    }
    return val;
}

int scale_mass_fractions(vector<double> &massf, double tolerance)
{
    size_t isp;
    int status_flag = SUCCESS;
    double massf_sum = 0.0;
    for ( isp = 0; isp < massf.size(); ++isp ) {
	massf[isp] = massf[isp] >= 0.0 ? massf[isp] : 0.0;
	massf_sum += massf[isp];
    }
    if ( fabs(massf_sum - 1.0) > 0.1 ) {
	/* Something disasterous must have happened to miss by this much. */
	status_flag = FAILURE;
    }
    if ( fabs(massf_sum - 1.0) > tolerance ) {
	for ( isp = 0; isp < massf.size(); ++isp ) massf[isp] /= massf_sum;
    }
    return status_flag;
} // end scale_mass_fractions()


int copy_mass_fractions(const vector<double> &f_src, vector<double> &f_target)
{
    if( f_src.size() != f_target.size() ) {
	cerr << "Mismatched sizes of valarrays for f_src "
	     << " and f_target in copy_mass_fraction()" << endl;
	return (-1);
    }
    for ( size_t isp = 0; isp < f_src.size(); ++isp ) f_target[isp] = f_src[isp];
    return SUCCESS;
} // end copy_mass_fractions()
