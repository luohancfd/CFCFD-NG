/**
 * bbla.d
 * Bare-bones linear algebra functions.
 *
 * Author: Peter J.
 * Version: 2014-06-28, just enough to get our CFD code going.
 */

import std.conv;
import std.algorithm;
import std.math;

class Matrix {
    size_t _nrows;
    size_t _ncols;
    double[][] _data;

    this(size_t n) {
	_nrows = n;
	_ncols = n;
	_data.length = n;
	foreach(i; 0 .. n) _data[i].length = n;
    }

    this(size_t nrows, size_t ncols) {
	_nrows = nrows;
	_ncols = ncols;
	_data.length = nrows;
	foreach(i; 0 .. nrows) _data[i].length = ncols;
    }

    this(in Matrix other) {
	this(other._nrows, other._ncols);
	foreach(row; 0 .. _nrows)
	    foreach(col; 0 .. _ncols)
		_data[row][col] = other._data[row][col];
    }

    this(in double[] vec) {
	this(vec.length, 1);
	foreach(row; 0 .. _nrows) _data[row][0] = vec[row];
    }

    this(in double[][] other) {
	this(other.length, other[0].length);
	foreach(row; 0 .. _nrows)
	    foreach(col; 0 .. _ncols)
		_data[row][col] = other[row][col];
    }

    @property const size_t nrows() { return _nrows; }
    @property const size_t ncols() { return _ncols; }

    const double opIndex(size_t row, size_t col) {
	return _data[row][col];
    }

    ref double opIndexAssign(double c, size_t row, size_t col) {
	_data[row][col] = c;
	return _data[row][col];
    }

    ref double opIndexOpAssign(string op)(double c, size_t row, size_t col)
	if ( op == "+" )
    {
	_data[row][col] += c;
	return _data[row][col];
    }

    ref double opIndexOpAssign(string op)(double c, size_t row, size_t col)
	if ( op == "-" )
    {
	_data[row][col] -= c;
	return _data[row][col];
    }

    ref double opIndexOpAssign(string op)(double c, size_t row, size_t col)
	if ( op == "*" )
    {
	_data[row][col] *= c;
	return _data[row][col];
    }

    ref double opIndexOpAssign(string op)(double c, size_t row, size_t col)
	if ( op == "/" )
    {
	_data[row][col] /= c;
	return _data[row][col];
    }

    override string toString() {
	string s = "Matrix[";
	foreach(row; 0 .. _nrows) {
	    s ~= "[";
	    foreach(col; 0 .. _ncols) {
		s ~= to!string(_data[row][col]);
		if ( col < _ncols-1 ) s ~= ",";
	    }
	    s ~= "]";
	    if ( row < _nrows-1 ) s ~= ",";
	}
	s ~= "]";
	return s;
    }

    void swap_rows(size_t i1, size_t i2)
    {
	swap(_data[i1], _data[i2]);
    }
}

Matrix zeros(size_t rows, size_t cols)
{
    Matrix my_matrix = new Matrix(rows, cols);
    foreach(row; 0 .. rows) {
	foreach(col; 0 .. cols) {
	    my_matrix[row,col] = 0.0;
	}
    }
    return my_matrix;
}

Matrix eye(size_t n)
{
    Matrix ident_matrix = new Matrix(n);
    foreach(row; 0 .. n) {
	foreach(col; 0 .. n) {
	    ident_matrix[row,col] = (row == col) ? 1.0 : 0.0;
	}
    }
    return ident_matrix;
}

Matrix transpose(in Matrix other)
{
    Matrix my_matrix = new Matrix(other.ncols, other.nrows);
    foreach(row; 0 .. other.nrows) {
	foreach(col; 0 .. other.ncols) {
	    my_matrix[col,row] = other[row,col];
	}
    }
    return my_matrix;
}

Matrix dot(in Matrix a, in Matrix b)
{
    if ( a.ncols != b.nrows ) {
	throw new Exception("incompatible matrices for dot product");
    }
    size_t nrows = a.nrows;
    size_t ncols = b.ncols;
    Matrix c = zeros(nrows, ncols);
    foreach(row; 0 .. nrows) {
	foreach(col; 0 .. ncols) {
	    foreach(i; 0 .. a.ncols) {
		c[row,col] += a[row,i] * b[i,col];
	    }
	}
    }
    return c;
}

/**
 * Perform Gauss-Jordan elimination on an augmented matrix.
 * c = [A|b] such that the mutated matrix becomes [I|x]
 * where x is the solution vector(s) to A.x = b
 */
void gauss_jordan_elimination(ref Matrix c, double very_small_value=1.0e-16)
{
    if (c.ncols < c.nrows) {
	throw new Exception("too few columns supplied");
    }
    foreach(j; 0 .. c.nrows) {
	// Select pivot.
	size_t p = j;
	foreach(i; j+1 .. c.nrows) {
	    if ( abs(c[i,j]) > abs(c[p,j]) ) p = i;
	}
	if ( abs(c[p,j]) < very_small_value ) {
	    throw new Exception("matrix is essentially singular");
	}
	c.swap_rows(p,j);
	// Scale row j to get unity on the diagonal.
	double cjj = c[j,j];
	foreach(col; 0 .. c.ncols) c[j,col] /= cjj;
	// Do the elimination to get zeros in all off diagonal values in column j.
	foreach(i; 0 .. c.nrows) {
	    if ( i == j ) continue;
	    double cij = c[i,j];
	    foreach(col; 0 .. c.ncols) c[i,col] -= cij * c[j,col]; 
	}
    } // end foreach j
} // end gauss_jordan_elimination()
