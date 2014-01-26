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
 * \version 26-Nov-2013 - and then to vector (PJ)
 **/

#ifndef NF_LIN_ALG_HH
#define NF_LIN_ALG_HH

#include <iostream>
#include <string>
#include <vector>

#define PVT_DIM 1000  ///< Static working array size for Gaussian elimination method.

void add_vectors(std::vector<double> &c, const std::vector<double> &a, const std::vector<double> &b );
void subtract_vectors(std::vector<double> &c, const std::vector<double> &a, const std::vector<double> &b );
void elemul_vectors(std::vector<double> &c, const std::vector<double> &a, const std::vector<double> &b);
void elepow_vectors(std::vector<double> &c, const std::vector<double> &b, double a);
void scale_vector(std::vector<double> &a, double b );
void scale_vector2vector(const std::vector<double> &a, double b, std::vector<double> &c);
void copy_vector(const std::vector<double> &src, std::vector<double> &target);

void print_vector(const std::vector<double> &a);

double norm2_vector(const std::vector<double> &a,
		       int start=0, int finish=-1);
void merge_vectors(const std::vector<double> &src_A, const std::vector<double> &src_B, std::vector<double> &target );

void split_vector( const std::vector<double> &src, std::vector<double> &target_A, std::vector<double> &target_B );
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
    int extract_column(std::vector<double> &x,
		       size_t start_row, size_t column);
    int insert_column(const std::vector<double> &x,
		       size_t start_row, size_t column);
    int extract_minor(Valmatrix &x, size_t iA, size_t jA, size_t iB, size_t jB); 
    int insert_minor(const Valmatrix &x, size_t iA, size_t jA, size_t iB, size_t jB); 
    

    int eye(size_t n);
    int print_mat();
    void scale(double s)
    { scale_vector(data_, s); }
    
    // DFP 26-Dec-2014:
    // Although there is another Gaussian elimination function, this implementation needs
    // to be retained for use within programs compiled with OpenMPI (the row_index array
    // used in the function outside of this class is not thread safe)
    int gaussian_elimination( std::vector<double> &x, std::vector<double> &b, bool with_scaling=false );

private:
    size_t nrows_;
    size_t ncols_;
    std::vector<double> data_;
    std::vector<int> row_index_;

};

#ifndef SWIG
std::ostream& operator<<( std::ostream &os, const Valmatrix &v );
#endif
std::vector<double> vector_matmul( Valmatrix &matrix, std::vector<double> &b);
int vector_mul(Valmatrix &A, std::vector<double> &x, std::vector<double> &b);
int gaussian_elimination( Valmatrix &A, std::vector<double> &x,
			  std::vector<double> &c, bool with_scaling=false );

int householder_transform(const std::vector<double> &x,
			  std::vector<double> &v, double &tau);
int householder_hm(double tau, std::vector<double> &v, Valmatrix &A);
int householder_hv(double tau, const std::vector<double> &v,
		   std::vector<double> &w);
int householder_QR(Valmatrix &A, std::vector<double> &tau);
int least_squares_solve(Valmatrix &A, std::vector<double> &x, std::vector<double> &b);
int QR_solve(Valmatrix &A, std::vector<double> &x, std::vector<double> &b);
double eval_matrix_determinant( Valmatrix &A, size_t i=0 );
double sub_matrix_det( Valmatrix &A, size_t i_skip, size_t j_skip, size_t row=0 );
#endif
