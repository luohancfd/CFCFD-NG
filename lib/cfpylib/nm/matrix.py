## \file matrix.py
## \ingroup nm
##
## \brief Simple matrix functions.
##
## These functions provide basic matrix initialization and operations.
## The storage scheme is a list-of-lists containing the matrix elements
## in row-major form.  For example:   \n
##     A = [[ a_00 a_01 a_02 ]        \n
##          [ a_10 a_11 a_12 ]        \n
##          [ a_20 a_21 a_22 ]]       \n
## Note that indices that start from zero as is usual for Python arrays.
## The first index references the row and the second the column
## (i.e. the element within the row).
## Vectors are stored in the same way, as an n-by-1 column vector
## or a 1-by-n row vector.
##
## See chapter 2 of the text by Schilling & Harris (Applied Numerical Methods
## for Engineers) for detailed descriptions of the Gauss-Jordan elimination
## and determinant evaluation algorithms.
##
## These functions are intended to be simple to study and use rather than
## be very efficient.  Generally, copies of the user's matrices are manipulated
## and returned rather than mutating the originals. Although, this increases
## memory use and reduces speed, the functions involve fewer (unexpected)
## side effects.
##    
## \author PA Jacobs.
##         School of Engineering, The University of Queensland
##
## \version 27-Oct-03
##

#----------------------------------------------------------------------------
# House-keeping functions...

def new_matrix(nr, nc=1, fillValue=0.0):
    """Returns a new nr-by-nc matrix with all elements bound to the same value.

    nr : number of rows
    nc : number of columns
    """
    dmatrix = []
    for ir in range(nr):
        dmatrix.append([fillValue] * nc)
    return dmatrix

def zeros(nr, nc=1):
    "Returns a new nr-by-nc matrix of all 0.0 values"
    return new_matrix(nr, nc, 0.0)

def ones(nr, nc=1):
    "Returns a new nr-by-nc matrix of all 1.0 values"
    return new_matrix(nr, nc, 1.0)

def eye(n):
    "Returns the n-by-n identity matrix."
    imatrix = zeros(n, n)
    for i in range(n): imatrix[i][i] = 1.0
    return imatrix

from copy import copy

def copy_matrix(A):
    "Returns a (new) copy of the original matrix."
    B = []
    for ir in range( len(A) ):
        B.append( copy(A[ir]) )
    return B

def make_column_vector(lst):
    "Returns a column vector (nx1 matrix) created from the given list of numbers."
    assert type(lst) is list
    vc = []
    for num in lst:
        vc.append([num,])
    return vc

def get_column(A, jc):
    "Returns the elements of column j of matrix A."
    B = []
    for ir in range( len(A) ):
        B.append( A[ir][jc] )
    return B

def transpose(A):
    "Returns the transpose of the supplied matrix."
    B = []
    for jc in range( len(A[0]) ):
        B.append( get_column(A, jc) )
    return B

#-----------------------------------------------------------------------------
# Simple arithmetic functions...
# These functions return new matrices and leave the originals unchanged.

def mmult(A, B):
    "Matrix multiply"
    nrA,ncA = len(A),len(A[0])
    nrB,ncB = len(B),len(B[0])
    if ncA != nrB: raise IndexError, "Matrices are not compatible."
    C = []
    for ir in range(nrA):
        C.append([0.0] * ncB)
        for jc in range(ncB):
            sum = 0.0
            for k in range(ncA):
                sum += A[ir][k] * B[k][jc]
            C[ir][jc] = sum
    return C

def madd(A, B):
    "Matrix add (element by element)"
    nrA,ncA = len(A),len(A[0])
    nrB,ncB = len(B),len(B[0])
    if ncA != ncB or nrA != nrB:
        raise IndexError, "Matrices are not compatible."
    C = []
    for ir in range(nrA):
        C.append([0.0] * ncA)
        for jc in range(ncA):
            C[ir][jc] = A[ir][jc] + B[ir][jc]
    return C

def mscale(A, alpha):
    "Scale all elements of A by factor alpha."
    A = copy_matrix(A)
    for ir in range(len(A)):
        scale_row(A, ir, alpha)
    return A

def msub(A, B):
    "Matrix subtract (element by element) B from A"
    return madd( A, mscale(B, -1.0) )

#-----------------------------------------------------------------------------
# Solution of linear systems of equations...

# First, elementary row operations
# These functions alter the matrices.

def scale_row(A, k, alpha):
    "Scale the elements of the k-th row of matrix A by factor alpha"
    ncA = len(A[k])
    for jc in range(ncA): A[k][jc] *= alpha

def interchange_rows(A, k, i):
    "Swap the contents of rows k and i. (Alters matrix A)"
    A[k], A[i] = A[i], A[k]

def locate_pivot_row(A, k):
    """Locate the pivot element.

    This is done by searching column k for the largest element below the
    current row, k
    """
    pivot_mag = abs(A[k][k])
    pivot_row = k
    for ir in range(k+1, len(A)):
        this_mag = abs(A[ir][k])
        if this_mag > pivot_mag:
            pivot_row = ir
            pivot_mag = this_mag
    return pivot_row

def madd_rows(A, k, i, alpha):
    "Add alpha times row i to row k. (Alters matrix A)"
    nc = len(A[k])
    for jc in range(nc):
        A[k][jc] += alpha * A[i][jc]

def equilibriate_equations(A, b):
    "Scale coefficients of the equations to get a maximum magnitude of 1.0"
    nr, ncA = len(A), len(A[0])
    for ir in range(nr):
        big = abs(A[ir][0])
        for jc in range(ncA):
            big = max(big, abs(A[ir][jc]))
        scale_row(A, ir, 1.0/big)
        scale_row(b, ir, 1.0/big)

# then, the elimination procedure itself.

def gauss_jordan(A, b, preserve_original=1):
    """Solves the set of linear equations Ax = b by Gauss-Jordan elimination.

    A : square matrix
    b : right-hand-side column vector (or matrix with several columns)
        on return, this is also the solution, x
        
    Be aware that both A and b are mutated by this function
    if preserve_original is set to 0.
    """
    nr,nc = len(A),len(A[0])
    if nr != nc: raise IndexError, "Matrix is not square."

    if preserve_original:
        A = copy_matrix(A); b = copy_matrix(b)
    equilibriate_equations(A, b)

    # Start work on the elimination proper.
    for j in range(nr):
        # Select pivot row and make diagonal element unity.
        jpivot = locate_pivot_row(A, j)
        if jpivot > j:
            interchange_rows(A, j, jpivot); interchange_rows(b, j, jpivot)
        pivot = A[j][j]
        if abs(pivot) < 1.0e-15:
            print "Warning: very small pivot magnitude:", pivot
        scale_row(A, j, 1.0/pivot); scale_row(b, j, 1.0/pivot)
        # Now, eliminate all off diagonal elements in the j-th column.
        for i in range(nr):
            if i == j: continue
            alpha = -A[i][j]
            madd_rows(A, i, j, alpha); madd_rows(b, i, j, alpha)
    return b

#-----------------------------------------------------------------------------

def inverse(A, preserve_original=1):
    "Returns the inverse of (square matrix) A."
    return gauss_jordan(A, eye(len(A)), preserve_original)
    
#-----------------------------------------------------------------------------

def det(A):
    "Returns the determinant of matrix A computed via Gaussian elimination."
    nr,nc = len(A),len(A[0])
    if nr != nc: raise IndexError, "Matrix is not square."
    A = copy_matrix(A)
    sign = 1
    # Start work on the elimination proper.
    for j in range(nr):
        # Select pivot row and sub-diagonal elements.
        jpivot = locate_pivot_row(A, j)
        if jpivot > j:
            interchange_rows(A, j, jpivot)
            sign = -1 * sign
        pivot = A[j][j]
        if abs(pivot) < 1.0e-15:
            print "Warning: looks like matrix is singular; pivot=", pivot
            return 0.0
        # Now, eliminate lower elements in the j-th column.
        for i in range(j+1,nr):
            alpha = -A[i][j]/pivot
            madd_rows(A, i, j, alpha)
        # Form determinant from diagonal elements.
        d = sign
        for i in range(nr):
            d *= A[i][i]
    return d
    
#-----------------------------------------------------------------------------

def norm(A):
    "Returns the row-sum norm (for a matrix) or infinity norm (for a vector)."
    nr, nc = len(A), len(A[0])
    v = []
    if nr > 1 and nc > 1:
        # Compute the row-sums of the matrix.
        for ir in range(nr):
            sum = 0.0
            for jc in range(nc):
                sum += abs( A[ir][jc] )
            v.append(sum)
    elif nc > 1:
        # row vector
        for jc in range(nc): v.append( abs(A[0][jc]) )
    else:
        # column vector (or a single element)
        for ir in range(nr): v.append( abs(A[ir][0]) )
    return max( v )


def cond(A):
    "Returns the condition number of matrix A, assuming A is nonsingular."
    return norm(A) * norm(inverse(A))
    
#-----------------------------------------------------------------------------

if __name__ == "__main__":
    print "Begin matrix test..."

    print "Create matrices"
    a = eye(3)
    print "a=", a
    b = new_matrix(3, 1, 4.5)
    print "b=", b

    print "Multiply matrices"
    c = mmult(a, b)
    print "c=a*b=", c
    # c2 = mmult(b, a)

    print "Try Gauss-Jordan elimination (Schilling & Harris Example 2.2.2)"
    a = [[1.0, -1.0, 0.0], [-2.0, 2.0, -1.0], [0.0, 1.0, -2.0]]
    b = [[2.0], [-1.0], [6.0]]
    print "before elimination:"
    print "a=", a
    print "b=", b
    x = gauss_jordan(a, b)
    print "after elimination:"
    print "a=", a
    print "b=", b
    print "x=", x
    print "residual=", msub(mmult(a,x), b)

    print "Transpose matrices and vectors"
    print "a=", a
    print "a^T=", transpose(a)
    print "x^T=", transpose(x)
    print "x^T * a^T=", mmult(transpose(x),transpose(a))

    print "Try to compute inverse (Schilling & Harris Example 2.2.5)"
    a = [[1.0, -1.0, 0.0], [2.0, 0.0, 4.0], [0.0, 2.0, -1.0]]
    ainv = inverse(a)
    print "a=", a
    print "ainv=", ainv
    print "ainv * a=", mmult(ainv,a)

    print "Determinant (Schilling & Harris Example 2.3.2)"
    a = [[1.0, 4.0, 0.0], [0.0, 2.0, 6.0], [-1.0, 0.0, 1.0]]
    print "a=", a
    print "det(a)=", det(a)

    print "Norms (Schilling & Harris Example 2.5.2)"
    a = [[4.0, -6.0, 9.0], [5.0, 8.0, -2.0], [7.0, -3.0, 1.0]]
    b = transpose([[-1.0, 12.0, 13.0]])
    x = mmult(inverse(a), b)
    print "a=", a, "norm(a)=", norm(a)
    print "b=", b, "norm(b)=", norm(b)
    print "x=", x, "norm(x)=", norm(x)

    print "Condition number (Schilling & Harris Example 2.5.3)"
    a = [[100.0, -200.0], [-200.0, 401.0]]
    print "a=", a, " cond(a)=", cond(a)
    
    print "Done."
