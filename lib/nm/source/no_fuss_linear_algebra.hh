/** \file no_fuss_linear_algebra.hh
 * \brief A no fuss linear algebra module.
 *
 * I wanted/needed a some basic linear algebra and
 * some simple vectors and matrices without the fuss
 * and overheads of the many linear algebra libraries
 * available for C++, so I wrote my own.
 *
 * This also limits the reliance on third-party software.
 *
 * \author Rowan J Gollan
 * \date 22-Apr-2006
 * \version 23-Sep-2006 - implementation changed to valarrays
 **/

#ifndef NF_LIN_ALG_HH
#define NF_LIN_ALG_HH

#include <iostream>
#include <string>
#include <valarray>

#define PVT_DIM 1000  ///< Static working array size for Gaussian elimination method.

void add_valarrays(std::valarray<double> &c, const std::valarray<double> &a, const std::valarray<double> &b );
void subtract_valarrays( std::valarray<double> &c, const std::valarray<double> &a, const std::valarray<double> &b );
void scale_valarray( std::valarray<double> &a, double b );
void scale_valarray2valarray( const std::valarray<double> &a, double b, std::valarray<double> &c);
void copy_valarray( const std::valarray<double> &src, std::valarray<double> &target );

void print_valarray( const std::valarray<double> &a );

double norm2_valarray( const std::valarray<double> &a,
		       int start=0, int finish=-1 );

void merge_valarrays(const std::valarray<double> &src_A, const std::valarray<double> &src_B, std::valarray<double> &target );

void split_valarray( const std::valarray<double> &src, std::valarray<double> &target_A, std::valarray<double> &target_B );

class Valmatrix {
public:
    /// \brief Default constructor
    Valmatrix();

    /// \brief Construct matrix with nrows and ncols
    Valmatrix(size_t nrows, size_t ncols);

    /// \brief Construct matrix (nrows, ncols) and fill with double
    Valmatrix(size_t nrows, size_t ncols, double fill_val);

    /// \brief Copy constructor
    Valmatrix(const Valmatrix &v );

    /// \brief Default destructor
    ~Valmatrix();

    Valmatrix& operator=(const Valmatrix &v);

    /// \brief Resize a matrix
    void resize(size_t nrows, size_t ncols);

    /// \brief string representation
    std::string str() const;

    size_t nrows() const { return nrows_; }
    size_t ncols() const { return ncols_; }

    inline double get(size_t i, size_t j) const
    { return data_[i*ncols_+j]; }

    inline void set(size_t i, size_t j, double val)
    { data_[i*ncols_+j] = val; }
    //--- Some special functions for vectors ---//
    /// \brief provide subscripting as per normal arrays
    //double operator() (size_t i, size_t j) { return data_[i*ncols_+j]; }
    //double& operator() (size_t i, size_t j) { return data_[i*ncols_+j]; }
    
#ifndef SWIG
    //    double& operator() (size_t i, size_t j) const { return &(data_[i*ncols_+j]); }
#endif
    int extract_column(std::valarray<double> &x,
		       size_t start_row, size_t column);
    int insert_column(const std::valarray<double> &x,
		       size_t start_row, size_t column);
    int extract_minor(Valmatrix &x, size_t iA, size_t jA, size_t iB, size_t jB); 
    int insert_minor(const Valmatrix &x, size_t iA, size_t jA, size_t iB, size_t jB); 
    

    int eye(size_t n);

    void scale(double s)
    { scale_valarray(data_, s); }
    
    int gaussian_elimination( std::valarray<double> &x, std::valarray<double> &b, bool with_scaling=false );

private:
    size_t nrows_;
    size_t ncols_;
    std::valarray<double> data_;
    std::valarray<int> row_index_;

};

#ifndef SWIG
std::ostream& operator<<( std::ostream &os, const Valmatrix &v );
#endif


int gaussian_elimination( Valmatrix &A, std::valarray<double> &x,
			  std::valarray<double> &c, bool with_scaling=false );

int householder_transform(const std::valarray<double> &x,
			  std::valarray<double> &v, double &tau);
int householder_hm(double tau, std::valarray<double> &v, Valmatrix &A);
int householder_hv(double tau, const std::valarray<double> &v,
		   std::valarray<double> &w);
int householder_QR(Valmatrix &A, std::valarray<double> &tau);
int least_squares_solve(Valmatrix &A, std::valarray<double> &x, std::valarray<double> &b);
int QR_solve(Valmatrix &A, std::valarray<double> &x, std::valarray<double> &b);
double eval_matrix_determinant( Valmatrix &A, size_t i=0 );
double sub_matrix_det( Valmatrix &A, size_t i_skip, size_t j_skip, size_t row=0 );
#endif
