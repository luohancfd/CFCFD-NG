/** \file no_fuss_linear_algebra.cxx
 * \brief A no fuss linear algrbra module.
 *
 * \author Rowan J Gollan
 * \date 22-Apr-2006
 *
 **/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include "../../util/source/useful.h"
#include "no_fuss_linear_algebra.hh"

using namespace std;

static int row_index[PVT_DIM];
const double essentially_zero = 1.0e-200;

/// \brief Add vectors a and b, place result in c
///
/// c = a + b
/// This is a fast version for intensive numerical
/// work as no memory allocation is involved.
///
/// \author Rowan J Gollan
/// \version 23-Apr-2006
/// \version 23-Sep-2006 - changed to valarrays
/// \version 26-Nov-2013 - then to vector (PJ)
///
void add_vectors(vector<double> &c, const vector<double> &a, const vector<double> &b )
{
    
    if( (a.size() != b.size()) || (a.size() != c.size()) ) {
	cerr << "Vector dimensions do not agree.\n";
	cerr << "While attempting add_NFVecs,\n";
	cerr << "size of vector a = " << a.size() << endl;
	cerr << "size of vector b = " << b.size() << endl;
	cerr << "size of vector c = " << c.size() << endl;
	cerr << "Bailing out.\n";
	exit(MISMATCHED_DIMENSIONS);
    }
    
    for( size_t i = 0; i < a.size(); ++i ) c[i] = a[i] + b[i];
}

/// \brief Subtract vector b from a, place result in c
///
/// c = a - b
/// This is a fast version for intensive numerical
/// work as no memory allocation is involved.
///
/// \author Rowan J Gollan
/// \version 23-Apr-2006
///
void subtract_vectors(vector<double> &c, const vector<double> &a, const vector<double> &b )
{
    if( (a.size() != b.size()) || (a.size() != c.size()) ) {
	cerr << "Vector dimensions do not agree.\n";
	cerr << "While attempting add_NFVecs,\n";
	cerr << "size of vector a = " << a.size() << endl;
	cerr << "size of vector b = " << b.size() << endl;
	cerr << "size of vector c = " << c.size() << endl;
	cerr << "Bailing out.\n";
	exit(MISMATCHED_DIMENSIONS);
    }
    
    for( size_t i = 0; i < a.size(); ++i ) c[i] = a[i] - b[i];
}

/// brief Multiply two arrays element-wise
// 
/// c = [a]*[b]
// 
/// Follows the same logic as all the code written by Rowan Gollan
/// for adding and subtracting.
// 
// \author Jared Clifford
// \version 11 September 2013
void elemul_vectors(vector<double> &c, const vector<double> &a, const vector<double> &b)
{
  if ( (a.size() != b.size()) || (a.size() != c.size()) ){
    cerr << "Vector dimensions do not agree .\n";
    cerr << "While attempting elemul_vector, \n";
    cerr << "size of vector a = " << a.size() << endl;
    cerr << "size of vector b = " << b.size() << endl;
    cerr << "size of vector c = " << c.size() << endl;
    cerr << "Bailing out. \n";
    exit(MISMATCHED_DIMENSIONS);
  }

  for ( size_t i = 0; i < a.size(); ++i ) c[i] = a[i] * b[i];

}


/// brief Raise each element of an array by a certain power
// 
/// c = [b]**a
// 
/// Follows the same logic as all the code written by Rowan Gollan
/// for adding and subtracting.
// 
// \author Jared Clifford
// \version 11 September 2013
void elepow_vectors(vector<double> &c, const vector<double> &b, double a)
{
    for ( size_t i = 0; i < c.size(); ++i ) c[i] = pow(b[i], a);
}

/// \brief Scale vector a by b
///
/// a = b*a
/// This is a fast version for intensive numerical
/// work as no memory allocation is involved.
///
/// \author Rowan J Gollan
/// \version 23-Apr-2006
///
void scale_vector(vector<double> &a, double b)
{
    for ( size_t i = 0; i < a.size(); ++i ) a[i] *= b;
}

/// \brief Scale vector a by b and place result in c
///
/// c = b*a
/// This is a fast version for intensive numerical
/// work as no memory allocation is involved.
///
/// \author Rowan J Gollan
/// \version 23-Apr-2006
///
void scale_vector2vector( const vector<double> &a, double b, vector<double> &c )
{
    for ( size_t i = 0; i < a.size(); ++i ) c[i] = b*a[i];
}

void copy_vector(const vector<double> &src, vector<double> &target )
{
    if ( src.size() != target.size() ) {
	cerr << "Vector dimensions do not agree.\n";
	cerr << "While attempting to copy two vectors,\n";
	cerr << "size of vector src = " << src.size() << endl;
	cerr << "size of vector target = " << target.size() << endl;
	cerr << "Bailing out.\n";
	exit(MISMATCHED_DIMENSIONS);
    }
	
    for ( size_t i = 0; i < src.size(); ++i ) target[i] = src[i];
    
}

void print_vector( const vector<double> &a )
{
    if ( a.size() == 0 ) return;
    cout << "[ ";
    for ( size_t i = 0; i < a.size()-1; ++i) cout << a[i] << ", ";
    cout << a[a.size()-1] << " ]" << endl;
}

double norm2_vector( const vector<double> &a, int start, int finish )
{
    if ( finish < 0 )
	finish = a.size();

    // Check finish > start
    if ( finish < start ) {
	cout << "Array stride start:finish == " << start << ":" << finish << endl;
	cout << "is not appropriate.\n";
	cout << "Bailing out!\n";
	exit(NUMERICAL_ERROR);
    }

    double norm = 0.0;
    for ( int i = start; i < finish; ++i ) {
	norm += a[i]*a[i];
    }
    return sqrt(norm);
}

void merge_vectors(const vector<double> &src_A, 
    		     const vector<double> &src_B,
		           vector<double> &target )
{
    // Merge two vectors, one after the other
    if ( ( src_A.size() + src_B.size() ) != target.size() ) {
	cerr << "Vector dimensions do not agree.\n";
	cerr << "While attempting to merge two vectors,\n";
	cerr << "size of vector src_A = " << src_A.size() << endl;
	cerr << "size of vector src_B = " << src_B.size() << endl;
	cerr << "size of vector target = " << target.size() << endl;
	cerr << "Bailing out.\n";
	exit(MISMATCHED_DIMENSIONS);
    }
	
    for ( size_t i = 0; i < target.size(); ++i ) {
	if ( i < src_A.size() )
	    target[i] = src_A[i];
	else
	    target[i] = src_B[i-src_A.size()];
    }
    
}

void split_vector( const vector<double> &src, 
                     	   vector<double> &target_A, 
		           vector<double> &target_B )
{
    // Split a vector into 2 parts defined by vector dimensions
    if ( ( target_A.size() + target_B.size() ) != src.size() ) {
	cerr << "Vector dimensions do not agree.\n";
	cerr << "While attempting to merge two vectors,\n";
	cerr << "size of vector target_A = " << target_A.size() << endl;
	cerr << "size of vector target_B = " << target_B.size() << endl;
	cerr << "size of vector src = " << src.size() << endl;
	cerr << "Bailing out.\n";
	exit(MISMATCHED_DIMENSIONS);
    }
	
    for ( size_t i = 0; i < src.size(); ++i ) {
	if ( i < target_A.size() )
	    target_A[i] = src[i];
	else
	    target_B[i-target_A.size()] = src[i];
    }
}

Valmatrix::Valmatrix()
    : nrows_(1), ncols_(1)
{
    data_.resize(1);
    row_index_.resize(1);
}

/// \brief Construct matrix with nrows and ncols
///
/// \author Rowan J Gollan
/// \version 23-Apr-2006
///
Valmatrix::Valmatrix(size_t nrows, size_t ncols )
    : nrows_(nrows), ncols_(ncols)
{
    data_.resize(nrows*ncols);
    row_index_.resize(nrows);
}

/// \brief Construct matrix (nrows, ncols) and fill with double
///
/// \author Rowan J Gollan
/// \version 23-Apr-2006
///
Valmatrix::Valmatrix( size_t nrows, size_t ncols, double fill_val )
    : nrows_(nrows), ncols_(ncols)
{
    data_.resize(nrows*ncols);
    row_index_.resize(nrows);
    for( size_t i = 0; i < data_.size(); ++i ) data_[i] = fill_val;
}


/// \brief Copy constructor
///
/// \author Rowan J Gollan
/// \version 23-Apr-2006
///
Valmatrix::Valmatrix( const Valmatrix &v )
    : nrows_(v.nrows_), ncols_(v.ncols_)
{
    data_.resize(nrows_*ncols_);
    row_index_.resize(nrows_);
    for( size_t i = 0; i < data_.size(); ++i) data_[i] = v.data_[i];
}

///  \brief Default destructor
///
/// \author Rowan J Gollan
/// \version 23-Apr-2006
///
Valmatrix::~Valmatrix() {}


/// \brief Assignment operator
///
/// \author Rowan J. Gollan
/// \version 21-May-2007
Valmatrix&
Valmatrix::
operator=(const Valmatrix &v)
{
    if( this == &v ) { // Avoid aliasing
	return *this;
    }
    // Clear out memory by resizing.
    nrows_ = v.nrows_;
    ncols_ = v.ncols_;
    data_.resize(nrows_*ncols_);
    for( size_t i = 0; i < data_.size(); ++i )data_[i] = v.data_[i];
    return *this;
}

void Valmatrix::resize(size_t nrows, size_t ncols)
{
    data_.resize(nrows*ncols, 0.0);
    row_index_.resize(nrows);
    nrows_ = nrows;
    ncols_ = ncols;
}

/// \brief string representation
///
/// \author Rowan J Gollan
/// \version 23-Apr-2006
///
string Valmatrix::str() const
{
    ostringstream ost;
    ost << setprecision(6) << showpoint;
    ost << "(";
    for( size_t j = 0; j < ncols_; ++j ) ost << "  " << get(0,j);
    for( size_t i = 1; i < nrows_; ++i ) {
	ost << "\n ";
	for( size_t j = 0; j < ncols_; ++j ) ost << "  " << get(i,j);
    }
    ost << "  )\n";
    return ost.str();
}


int
Valmatrix::
extract_column(vector<double> &x, size_t start_row, size_t column)
{
    size_t nr = nrows() - start_row;
    if( nr != x.size() ) {
	cout << "Dimensions do not agree in Valmatrix::extract_column()\n";
	cout << "Bailing Out!\n";
	exit(MISMATCHED_DIMENSIONS);
    }

    for( size_t i = start_row; i < nrows(); ++i ) x[i-start_row] = get(i,column);

    return SUCCESS;
}

int
Valmatrix::
insert_column(const vector<double> &x, size_t start_row, size_t column)
{
    size_t nr = nrows() - start_row;
    if( nr != x.size() ) {
	cout << "Dimensions do not agree in Valmatrix::extract_column()\n";
	cout << "Bailing Out!\n";
	exit(MISMATCHED_DIMENSIONS);
    }

    for( size_t i = start_row; i < nrows(); ++i ) set(i, column, x[i-start_row]);

    return SUCCESS;
}

int
Valmatrix::
extract_minor(Valmatrix &x, size_t iA, size_t jA,
	       size_t iB, size_t jB)
{
    for( size_t i = iA; i <= iB; ++i ) {
	for( size_t j = jA; j <= jB; ++j ) {
	    x.set(i-iA,j-jA, get(i, j));
	}
    }
    return SUCCESS;
}
int
Valmatrix::
insert_minor(const Valmatrix &x, size_t iA, size_t jA,
	       size_t iB, size_t jB)
{
    for( size_t i = iA; i <= iB; ++i ) {
	for( size_t j = jA; j <= jB; ++j ) {
	    set(i, j, x.get(i-iA,j-jA));
	}
    }
    return SUCCESS;
}

int
Valmatrix::
eye(size_t n)
{
    resize(n, n);
    for(size_t i = 0; i < n; ++i ) {
	set(i, i, 1.0); 
    }
    return SUCCESS;
}

// \brief Print out a matrix out on screen in a way that's easy to read.  
// 
// For debugging only
// 
// \author Jared Clifford
// \verion 12 September 2013
// 
int
Valmatrix::
print_mat()
{
    size_t i,j;
    int width = 3;
    cout << endl << "[ ";
    for ( i = 0; i < this->ncols(); i++) {
	if (i != 0)
	    cout << "  ";
	cout << "[ ";
	for (j = 0; j < this->nrows();j++) {
	    cout << setprecision(3) << setw(width) << this->get(i,j) << " ";
	}
	cout << "] ";
	if (i != this->ncols()-1)
	    cout << endl;
    }
    cout << " ]"<< endl;
    return SUCCESS;
}

vector<double> vector_matmul( Valmatrix &matrix, vector<double> &b)
{
  // Ensure the matrix is square and matrix and vector dimensions match
  if ( matrix.ncols() != matrix.nrows() || matrix.ncols() != b.size() ){
		cerr << "vector_mul\n";
		cerr << "ERROR: Dims don't match.\n";
		cerr << "ncols= " << matrix.ncols() << " nrows= " << matrix.nrows() << " b.size() = " << b.size() << endl;
		cerr << "Bailing out.";
		exit(MISMATCHED_DIMENSIONS);
  }

  vector<double> store;
  store.resize(b.size(), 0.0);

  for ( size_t i = 0; i < b.size(); i++){
      for ( size_t j = 0; j < matrix.nrows(); j++){
	  store[i] = store[i] + matrix.get(i,j)*b[j];
      }
  }

  return store;
}

// \brief Multiply a matrix by a vector.  Used to solve the wall equation.
// 
// Solves A.x = b, where A and x are known
// 
// A is a matrix, x is a vector, result is placed in b
// 
// \author Jared Clifford
// \version 12 September 2013
int vector_mul(Valmatrix &A, vector<double> &x, vector<double> &b)
{
    // Ensure the matrix is square and matrix and vector dimensions match
    if ( A.ncols() != A.nrows() || A.ncols() != x.size() || A.ncols() != b.size() ) {
	cerr << "vector_mul\n";
	cerr << "ERROR: Dims don't match.\n";
	cerr << "ncols= " << A.ncols() << " nrows= " << A.nrows() << " x.size() = " << x.size() << " b.size() = " << b.size() << endl;
	cerr << "Bailing out.";
	exit(MISMATCHED_DIMENSIONS);
    }
    
    for ( size_t i = 0; i < b.size(); ++i ) {
	b[i] = 0.0;
	for ( size_t j = 0; j < A.nrows(); j++ ) {
	    b[i] = b[i] + A.get(i,j)*x[j];
	}
    }

    return SUCCESS;
}

// int
// Valmatrix::
// gaussian_elimination( std::vector<double> &x, vector<double> &b, bool with_scaling )
// {
//     int i, j, k, p, tmp;
//     double max, sum, factor;
//     // Ensure the matrix is square
//     if( this->ncols() != this->nrows() ) {
// 	cerr << "gaussian_elimination():\n";
// 	cerr << "ERROR: A square coefficient matrix is not supplied.\n";
// 	cerr << "ncols= " << this->ncols() << " nrows= " << this->nrows() << endl;
// 	cerr << "Bailing out.";
// 	exit(MISMATCHED_DIMENSIONS);
//     }

//     // Setup an index into the rows avoiding the overhead
//     // of swapping rows when required.
//     for( i = 0; i < int(this->nrows()); ++i ) row_index_[i] = i;
    
//     // Elimination and partial-pivoting
//     for( i = 0; i < int(this->nrows()-1); ++i ) {

// 	p = i;
// 	max = fabs(this->get(row_index_[p],i));
// 	for( j = i+1; j < int(this->nrows()); ++j ) {
//  	    if( fabs(this->get(row_index_[j],i)) > max) {
// 		max = fabs(this->get(row_index_[j],i));
// 		p = j;
// 	    }
// 	}

// 	if( fabs(this->get(row_index_[p],i)) < essentially_zero ) {
// 	    cout << "Valmatrix::gaussian_elimination()\n";
// 	    cout << "No unique solution exists when using gaussian_elimination()\n";
// 	    return NUMERICAL_ERROR;
// 	}

// 	if( row_index_[i] != row_index_[p] ) {
// 	    // Perform swap of rows
// 	    tmp = row_index_[i]; row_index_[i] = row_index_[p]; row_index_[p] = tmp;
// 	}

// 	for( j = i+1; j < int(this->nrows()); ++j ) {
// 	    factor = this->get(row_index_[j],i)/this->get(row_index_[i],i);
// 	    for( k = i+1; k < int(this->nrows()); ++k ) {
// 		this->set(row_index_[j],k, this->get(row_index_[j],k) - factor*this->get(row_index_[i],k));
// 	    }
// 	    b[row_index_[j]] = b[row_index_[j]] - factor*b[row_index_[i]];
// 	}
//     }

//     if( fabs(this->get(row_index_[this->nrows()-1],this->nrows()-1)) < essentially_zero ) {
//     	cout << "Valmatrix::gaussian_elimination()\n";
// 	cout << "No unique solution exists when using gaussian_elimination()\n";
// 	return NUMERICAL_ERROR;
//     }

//     // Start backward substution
//     x[x.size()-1] = b[row_index_[b.size()-1]]/this->get(row_index_[this->nrows()-1],this->nrows()-1);

//     for( i = int(this->nrows()-2); i >= 0; --i ) {
// 	sum = 0.0;
// 	for( j = i+1; j < int(this->nrows()); ++j ) {
// 	    sum += this->get(row_index_[i],j)*x[j];
// 	}
// 	x[i] = (b[row_index_[i]] - sum)/this->get(row_index_[i],i);
//     }

//     return SUCCESS;

// }

ostream& operator<<( ostream &os, const Valmatrix &v )
{
    os << v.str();
    return os;
}

int gaussian_elimination( Valmatrix &A, vector<double> &x, vector<double> &b, bool with_scaling )
{
    int i, j, k, p, tmp;
    double max, sum, factor;
    // Ensure the matrix is square
    if( A.ncols() != A.nrows() ) {
	cerr << "gaussian_elimination():\n";
	cerr << "ERROR: A square coefficient matrix is not supplied.\n";
	cerr << "ncols= " << A.ncols() << " nrows= " << A.nrows() << endl;
	cerr << "Bailing out.";
	exit(MISMATCHED_DIMENSIONS);
    }

    // Setup an index into the rows avoiding the overhead
    // of swapping rows when required.
    for( i = 0; i < int(A.nrows()); ++i ) row_index[i] = i;
    
    // Elimination and partial-pivoting
    for( i = 0; i < int(A.nrows()-1); ++i ) {

	p = i;
	max = fabs(A.get(row_index[p],i));
	for( j = i+1; j < int(A.nrows()); ++j ) {
 	    if( fabs(A.get(row_index[j],i)) > max) {
		max = fabs(A.get(row_index[j],i));
		p = j;
	    }
	}

	if( fabs(A.get(row_index[p],i)) < essentially_zero ) {
	    cout << "No unique solution exists when using gaussian_elimination()\n";
	    return NUMERICAL_ERROR;
	}

	if( row_index[i] != row_index[p] ) {
	    // Perform swap of rows
	    tmp = row_index[i]; row_index[i] = row_index[p]; row_index[p] = tmp;
	}

	for( j = i+1; j < int(A.nrows()); ++j ) {
	    factor = A.get(row_index[j],i)/A.get(row_index[i],i);
	    for( k = i+1; k < int(A.nrows()); ++k ) {
		A.set(row_index[j],k, A.get(row_index[j],k) - factor*A.get(row_index[i],k));
	    }
	    b[row_index[j]] = b[row_index[j]] - factor*b[row_index[i]];
	}
    }

    if( fabs(A.get(row_index[A.nrows()-1],A.nrows()-1)) < essentially_zero ) {
	cout << "No unique solution exists when using gaussian_elimination()\n";
	return NUMERICAL_ERROR;
    }

    // Start backward substution
    x[x.size()-1] = b[row_index[b.size()-1]]/A.get(row_index[A.nrows()-1],A.nrows()-1);

    for( i = int(A.nrows()-2); i >= 0; --i ) {
	sum = 0.0;
	for( j = i+1; j < int(A.nrows()); ++j ) {
	    sum += A.get(row_index[i],j)*x[j];
	}
	x[i] = (b[row_index[i]] - sum)/A.get(row_index[i],i);
    }

    return SUCCESS;

}

//------------------------------------------
// The following functions are modifications
// of code appearing in the Gnu Scientific
// Library (GSL).  The GSL is licensed under
// GPLv3 as far as I could ascertain from 
// the website.  As such, I am permitted to 
// make modifications to the code.  There
// are numerous modifications: changes of 
// data structures, and "re-packaging"/
// re-ordering of code.
//
// The present (15-Feb-2008) intended use is
// internal.  Therefore, these modifications
// do not need to be made public.  I have 
// based this on the following item from
// the FAQ for GPLv3.
//
// -----
// Does the GPL require that source code of modified versions
// be posted to the public?
//
// The GPL does not require you to release your modified version,
// or any part of it. You are free to make modifications and use
// them privately, without ever releasing them. This applies to 
// organizations (including companies), too; an organization can
// make a modified version and use it internally without ever
// releasing it outside the organization.
//
// But if you release the modified version to the public in some
// way, the GPL requires you to make the modified source code
// available to the program's users, under the GPL.
//
// Thus, the GPL gives permission to release the modified program
// in certain ways, and not in other ways; but the decision of
// whether to release it is up to you.
// ----
//
// In any case, all of our source code is distributed to our
// users.


int householder_transform(const vector<double> &x,
			  vector<double> &v, double &tau)
{
    if( x.size() != v.size() ) {
	cout << "Vector dimensions do not agree in\n";
	cout << "householder_transform()\n";
	cout << "x.size()= " << x.size() << " v.size()= " << v.size() << endl;
	cout << "Bailing out!\n";
	exit(MISMATCHED_DIMENSIONS);
    }

    size_t n = x.size();
    for(size_t i = 0; i < x.size(); ++i ) v[i] = x[i];
    
    if( n == 1 ) {
	tau = 0.0;
	return SUCCESS;
    }
    else {
	double alpha, beta;
	double xnorm = norm2_vector(x, 1);

	if( xnorm == 0.0 ) {
	    tau = 0.0;
	    return SUCCESS;
	}

	alpha = x[0];
	if( alpha >= 0.0 ) {
	    beta = -hypot(alpha, xnorm);
	}
	else {
	    beta = hypot(alpha, xnorm);
	}
	tau = (beta - alpha) / beta;
	scale_vector(v, 1.0/(alpha - beta));
	v[0] = beta;

	return SUCCESS;
    }
}

int householder_hm(double tau, vector<double> &v, Valmatrix &A)
{
    size_t i, j;
    
    for( j = 0; j < A.ncols(); ++j ) {
	     
	    
	double wj = A.get(0,j);
	
	for( i = 1; i < A.nrows(); ++i ) {
	    /* note, computed for v(0) = 1 above */
	    wj += A.get(i,j) * v[i];
	}
        
	/* Aij = Aij - tau vi wj */
        
	/* i = 0 */

	double A0j = A.get(0,j);
	A.set(0,j, A0j - tau *  wj);
        
	for( i = 1; i < A.nrows(); ++i ) {
	    double Aij = A.get(i, j);
	    double vi = v[i];
	    A.set(i, j, Aij - tau * vi * wj);
	}
    }
    return SUCCESS;
}

int householder_hv(double tau, const vector<double> &v,
		   vector<double> &w)
{
    const size_t N = v.size();
 
    if (tau == 0.0)
	return SUCCESS ;
    
    /* compute d = v'w */

    double d0 = w[0];
    double d1, d;

    vector<double> v1(N-1);
    vector<double> w1(N-1);

    for(size_t i = 1; i < N; ++i ) {
	v1[i-1] = v[i];
	w1[i-1] = w[i];
    }

    // Dot product.
    d1 = 0.0;
    for( size_t i = 0; i < v1.size(); ++i ) {
	d1 += v1[i]*w1[i];
    }
    
    d = d0 + d1;

    /* compute w = w - tau (v) (v'w) */
    double w0 = w[0];
    w[0] = w0 - tau*d;

    for( size_t i = 0; i < v1.size(); ++i ) {
	w1[i] = -tau*d*v1[i] + w1[i];
    }
    // Pack into w.
    for(size_t i = 1; i < N; ++i ) {
	w[i] = w1[i-1];
    }
    
    return SUCCESS;
}


int householder_QR(Valmatrix &A, vector<double> &tau)
{
    const size_t M = A.nrows();
    const size_t N = A.ncols();

    if( tau.size() != MINIMUM(M, N) ){
	cout << "The length of tau vector is not appropriate\n";
	cout << "in householder_QR().  It should be MINIMUM(M,N)\n";
	cout << "Bailing Out!\n";
	exit(MISMATCHED_DIMENSIONS);
    }
    
    size_t i;
    double tau_i;
    vector<double> v;
    vector<double> x;
    Valmatrix m;
    
    for( i = 0; i < MINIMUM(M, N); ++i ) {
	v.resize(M-i);
	x.resize(M-i);
	A.extract_column(x, i, i);
	householder_transform(x, v, tau_i);
	tau[i] = tau_i;
	A.insert_column(v, i, i);

	if( i+1 < N ) {
	    m.resize(M-i,N-i-1);
	    A.extract_minor(m, i, i+1, M-1, N-1);
	    householder_hm(tau_i, v, m);
	    A.insert_minor(m, i, i+1, M-1, N-1);
	}

    }
    return SUCCESS;
}


int QR_QTvec(Valmatrix &QR, vector<double> &tau, vector<double> &b)
{

    const size_t M = QR.nrows();
    const size_t N = QR.ncols();

    if (tau.size() != MINIMUM(M, N) ) {
	cout << "ERROR. Bailing Out!\n";
	exit(MISMATCHED_DIMENSIONS);
    }
    if( b.size() != M ) {
	cout << "ERROR.\n";
	exit(MISMATCHED_DIMENSIONS);
    }
    
    size_t i;

    vector<double> c(M);
    vector<double> h;
    vector<double> w;
    /* compute Q^T v */

    for (i = 0; i < MINIMUM(M, N); i++) {
	QR.extract_column(c, 0, i);

	h.resize(M-i);
	for( size_t ii = i; ii < M; ++ii ) h[ii-i] = c[ii];

	w.resize(M-i);
	for( size_t ii = i; ii < M; ++ii ) w[ii-i] = b[ii];

	double ti = tau[i];

	householder_hv(ti, h, w);

	// Copy back into b.
	for( size_t ii = i; ii < M; ++ii ) b[ii] = w[ii-i];
    }

    return SUCCESS;

}

int least_squares_solve(Valmatrix &A, vector<double> &x, vector<double> &b)
{
    const size_t M = A.nrows();
    const size_t N = A.ncols();

    if( M < N ) {
	cout << "M should be >= N in least_squares_solver()\n";
	cout << "M= " << M << " N= " << N << endl;
	cout << "Bailing Out!\n";
	exit(MISMATCHED_DIMENSIONS);
    }
    if( M != b.size() ) {
	cout << "M should be == b.size() in least_squares_solver()\n";
	cout << "M= " << M << " b.size()= " << b.size() << endl;
	cout << "Bailing Out!\n";
	exit(MISMATCHED_DIMENSIONS);
    }
    if( N != x.size() ) {
	cout << "N should be == x.size() in least_squares_solver()\n";
	cout << "N= " << N << " x.size()= " << x.size() << endl;
	cout << "Bailing Out!\n";
	exit(MISMATCHED_DIMENSIONS);
    }

    vector<double> tau(MINIMUM(M, N));

    // First find QR-factorization
    householder_QR(A, tau);
    // A is now in QR-form

    Valmatrix R(N, N);
    A.extract_minor(R, 0, 0, N-1, N-1);

    QR_QTvec(A, tau, b);

    for( size_t i = 0; i < x.size(); ++i ) x[i] = b[i];
    
    // Now just use back-substituion...

    x[N-1] = b[N-1]/R.get(N-1,N-1);

    for( int i = int(N-2); i >= 0; --i ) {
	double sum = 0.0;
	for( size_t j = i+1; j < N; ++j ) {
	    sum += R.get(i,j)*x[j];
	}
	x[i] = (b[i] - sum)/R.get(i,i);
    }
    return SUCCESS;
}


int QR_solve(Valmatrix &A, vector<double> &x, vector<double> &b)
{
    const size_t N = A.ncols();
    // Ensure the matrix is square
    if( A.ncols() != A.nrows() ) {
	cerr << "QR_solve():\n";
	cerr << "ERROR: A square coefficient matrix is not supplied.\n";
	cerr << "ncols= " << A.ncols() << " nrows= " << A.nrows() << endl;
	cerr << "Bailing out.";
	exit(MISMATCHED_DIMENSIONS);
    }

    vector<double> tau(b.size());

    // First find QR-factorization
    householder_QR(A, tau);

    QR_QTvec (A, tau, b);
    
    // And then back-substitute...
    x[N-1] = b[N-1]/A.get(N-1,N-1);

    for( int i = int(N-2); i >= 0; --i ) {
	double sum = 0.0;
	for( size_t j = i+1; j < N; ++j ) {
	    sum += A.get(i,j)*x[j];
	}
	x[i] = (b[i] - sum)/A.get(i,i);
    }
    return SUCCESS;
}

// A lazy determinent solver, not meant for large numerical usage
// DFP 21.10.2009

double eval_matrix_determinant( Valmatrix &A, size_t i )
{
    const size_t N = A.ncols();
    // Ensure the matrix is square
    if( A.ncols() != A.nrows() ) {
	cerr << "eval_matrix_determinant():\n";
	cerr << "ERROR: A square coefficient matrix is not supplied.\n";
	cerr << "ncols= " << A.ncols() << " nrows= " << A.nrows() << endl;
	cerr << "Bailing out.";
	exit(MISMATCHED_DIMENSIONS);
    }
    double det = 0.0;
    
    for ( size_t j=0; j<N; ++j )
    	det += ( ( 2*(i+j)/2==(i+j) ) ? 1.0 : -1.0 ) * A.get(i,j)  * sub_matrix_det( A, i, j );
    
    return det;
}

double sub_matrix_det( Valmatrix &A, size_t i_skip, size_t j_skip, size_t row )
{
    Valmatrix Asub(A.nrows()-1,A.ncols()-1);
    
    size_t i_sub = 0;
    
    for ( size_t i=0; i<A.nrows(); ++i ) {
    	if ( i==i_skip ) continue;
    	else {
    	    size_t j_sub = 0;
	    for ( size_t j=0; j<A.ncols(); ++j ) {
		if ( j==j_skip ) continue;
		else {
		    Asub.set(i_sub,j_sub,A.get(i,j));
		    j_sub++;
		}
	    }
	    i_sub++;
	}
    }
    
    double det = 0.0;
    
    size_t i = row;
    for ( size_t j=0; j<Asub.ncols(); ++j )
    	det += ( ( 2*(i+j)/2==(i+j) ) ? 1.0 : -1.0 ) * Asub.get(i,j) * sub_matrix_det( Asub, i, j );
    
    return det;
}
