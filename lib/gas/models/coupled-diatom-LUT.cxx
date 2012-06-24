// Author: Daniel F. Potter
// Version: 09-Apr-2012

#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "../../util/source/useful.h"
#include "coupled-diatom-LUT.hh"

using namespace std;

NoneqCoupledDiatomicLUT::NoneqCoupledDiatomicLUT( std::string fname, int icol )
 : fname_( fname ), icol_( icol )
{
    cout << "- Creating a new NoneqCoupledDiatomicLUT class" << endl;
    
    // Make the look-up-tables
    ifstream lut (fname.c_str());
    if ( !lut.is_open() ) {
    	cout << "NoneqCoupledDiatomicLUT::NoneqCoupledDiatomicLUT()" << endl
    	     << "Problem opening file: " << fname << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    // Read the file
    string line;
    istringstream iss;
    vector<double> data;
    
    // First line has the temperature list
    getline( lut, line );
    iss.str(line);
    copy( istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(T_list_) );
    ndim_ = int(T_list_.size());
    
    cout << "   * Cubic tabulated data table has " << ndim_ << " elements in each array" << endl;
    
    // initialise the data table
    data_ = new double**[ndim_];
    for ( int i=0; i<ndim_; ++i ) {
        data_[i] = new double*[ndim_];
        for ( int j=0; j<ndim_; ++j ) {
            data_[i][j] = new double[ndim_];
        }
    }
    
    // now the data
    int ie = 0, iv = 0, ir = 0;
    while ( getline( lut, line ) ) {
        iss.clear();
        data.clear();
        iss.str(line);
        copy( istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(data) );
        // cout << "ie = " << ie << ", iv = " << iv << ", ir = " << ir << endl;
        data_[ie][iv][ir] = data[icol_-1];
        ir++;
        if ( ir==ndim_ ) {
            ir=0; iv++;
            if ( iv==ndim_ ) {
                iv=0; ie++;
                // if ( ie==ndim_ ) {
                //    cout << "this should be the last line" << endl;
                // }
            }
        }
    }
}

NoneqCoupledDiatomicLUT::~NoneqCoupledDiatomicLUT()
{
    // clear the T_list vector
    T_list_.clear();
    
    // clear the data table
    for ( int i=0; i<ndim_; ++i ) {
        for ( int j=0; j<ndim_; ++j ) {
            delete data_[i][j];
        }
        delete data_[i];
    }
    delete data_;
}

string NoneqCoupledDiatomicLUT::str()
{
    ostringstream str;
    str << "NoneqCoupledDiatomicLUT instance made from column index: " << icol_ << " of file: " << fname_ << endl;
    
    return str.str();
}

void NoneqCoupledDiatomicLUT::find_bounding_temperature_indices( double T, int &i_l, int &i_u ) 
{
    i_l=-1, i_u=-1;
    if ( T >= T_list_.back() ) { i_l = ndim_-1; i_u = ndim_-1; }
    else if ( T < T_list_.front() ) { i_l = 0; i_u = 0; }
    else {
        for ( int i=0; i<ndim_-1; ++i ) {
            if ( T==T_list_[i] ) {
                i_l = i; i_u = i;
                break;
            }
            else if ( T > T_list_[i] && T < T_list_[i+1] ) {
                i_l = i; i_u = i + 1;
                break;
            }
        }
    }
    
    if ( i_l < 0 || i_u < 0 ) {
        cout << "NoneqCoupledDiatomicLUT::find_bounding_temperature_indices()" << endl
             << "Failed to find data for T = " << T << endl;
        exit( FAILURE );
    }
    
    // cout << "T = " << T << ", i_l = " << i_l << ", i_u = " << i_u << endl;
}

double NoneqCoupledDiatomicLUT::eval( double T_el, double T_vib, double T_rot )
{
    // locate 8 bounding data points for the given temperatures
    int ie_l, ie_u, iv_l, iv_u, ir_l, ir_u;
    find_bounding_temperature_indices(  T_el, ie_l, ie_u );
    find_bounding_temperature_indices( T_vib, iv_l, iv_u );
    find_bounding_temperature_indices( T_rot, ir_l, ir_u );
    
    // if all pairs are doubles, just return the data point
    if ( ie_l==ie_u && iv_l==iv_u && ir_l==ir_u ) return data_[ie_l][iv_l][ir_l];
    
    // trilinear interpolation of the 8 data points
    double c000 = data_[ie_l][iv_l][ir_l];
    double c001 = data_[ie_l][iv_l][ir_u];
    double c010 = data_[ie_l][iv_u][ir_l];
    double c011 = data_[ie_l][iv_u][ir_u];
    double c100 = data_[ie_u][iv_l][ir_l];
    double c101 = data_[ie_u][iv_l][ir_u];
    double c110 = data_[ie_u][iv_u][ir_l];
    double c111 = data_[ie_u][iv_u][ir_u];
    
    // electronic temperature interpolation
    double c00 = ( c100 - c000 ) / ( T_list_[ie_u] - T_list_[ie_l] ) * ( T_el - T_list_[ie_l] ) + c000;
    double c01 = ( c101 - c001 ) / ( T_list_[ie_u] - T_list_[ie_l] ) * ( T_el - T_list_[ie_l] ) + c001;
    double c11 = ( c111 - c011 ) / ( T_list_[ie_u] - T_list_[ie_l] ) * ( T_el - T_list_[ie_l] ) + c011;
    double c10 = ( c110 - c010 ) / ( T_list_[ie_u] - T_list_[ie_l] ) * ( T_el - T_list_[ie_l] ) + c010;
    
    // vibrational temperature interpolation
    double c0 = ( c10 - c00 ) / ( T_list_[iv_u] - T_list_[iv_l] ) * ( T_vib - T_list_[iv_l] ) + c00;
    double c1 = ( c11 - c01 ) / ( T_list_[iv_u] - T_list_[iv_l] ) * ( T_vib - T_list_[iv_l] ) + c01;
    
    // rotational temperature interpolation
    return ( c1 - c0 ) / ( T_list_[ir_u] - T_list_[ir_l] ) * ( T_rot - T_list_[ir_l] ) + c0;
}

EqCoupledDiatomicLUT::EqCoupledDiatomicLUT( std::string fname, int icol )
 : fname_( fname ), icol_( icol )
{
    // Make the look-up-tables
    ifstream lut (fname.c_str());
    if ( !lut.is_open() ) {
    	cout << "EqCoupledDiatomicLUT::EqCoupledDiatomicLUT()" << endl
    	     << "Problem opening file: " << fname << endl;
    	exit( 1 );
    }
    
    // Read the file
    string line;
    istringstream iss;
    vector<double> data;
    
    // First line has the temperature list
    getline( lut, line );
    iss.str(line);
    copy( istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(T_list_) );
    ndim_ = int(T_list_.size());    
    
    // now the data
    while ( getline( lut, line ) ) {
        iss.clear();
        data.clear();
        iss.str(line);
        copy( istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(data) );
        // cout << "ie = " << ie << ", iv = " << iv << ", ir = " << ir << endl;
        data_.push_back( data[icol_-1] );
        if ( int(data.size())==ndim_ ) break;
    }
}

EqCoupledDiatomicLUT::~EqCoupledDiatomicLUT()
{
    T_list_.clear();
    data_.clear();
}

string EqCoupledDiatomicLUT::str()
{
    ostringstream str;
    str << "EqCoupledDiatomicLUT instance made from column index: " << icol_ << " of file: " << fname_ << endl;
    
    return str.str();
}

void EqCoupledDiatomicLUT::find_bounding_temperature_indices( double T, int &i_l, int &i_u ) 
{
    i_l=-1, i_u=-1;
    if ( T >= T_list_.back() ) { i_l = ndim_-1; i_u = ndim_-1; }
    else if ( T < T_list_.front() ) { i_l = 0; i_u = 0; }
    else {
        for ( int i=0; i<ndim_-1; ++i ) {
            if ( T==T_list_[i] ) {
                i_l = i; i_u = i;
                break;
            }
            else if ( T > T_list_[i] && T < T_list_[i+1] ) {
                i_l = i; i_u = i + 1;
                break;
            }
        }
    }
    
    if ( i_l < 0 || i_u < 0 ) {
        cout << "EqCoupledDiatomicLUT::find_bounding_temperature_indices()" << endl
             << "Failed to find data for T = " << T << endl;
        exit( FAILURE );
    }
    
    // cout << "T = " << T << ", i_l = " << i_l << ", i_u = " << i_u << endl;
}

double EqCoupledDiatomicLUT::eval( double T )
{
    // locate nearest neighbours for the given temperatures
    int i_l, i_u;
    find_bounding_temperature_indices(  T, i_l, i_u );
    
    // if the two indices are the same just return the data point
    if ( i_l==i_u ) return data_[i_l];
    
    // linear interpolation of the two data points
    return ( data_[i_u] - data_[i_l] ) / ( T_list_[i_u] - T_list_[i_l] ) * ( T - T_list_[i_l] ) + data_[i_l];
}

