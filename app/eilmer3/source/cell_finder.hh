/** \file cell_finder.hh
 *  \ingroup eilmer3
 *  \brief Header file for the cell finder class used for ray-tracing.
 **/

#ifndef CELL_FINDER_HH
#define CELL_FINDER_HH

#include <vector>

#include "block.hh"
#include "cell.hh"

class CellFinder {
public:
    CellFinder( size_t nvertices );
    
    virtual ~CellFinder();
    
public:
    virtual int find_cell( const Vector3 &p, size_t &ib, size_t &ic, size_t &jc, size_t &kc ) = 0;
    
    virtual void test_cell( const FV_Cell * cell, const Vector3 &p, int *dc ) = 0;
    
protected:
    size_t nvertices_;
};


class CellFinder2D : public CellFinder {
public:
    CellFinder2D( size_t nvertices = 4 );
    
    ~CellFinder2D();

public:
    int find_cell( const Vector3 &p, size_t &ib, size_t &ic, size_t &jc, size_t &kc );
    
    void test_cell( const FV_Cell * cell, const Vector3 &p, int *dc );
    
private:
    // NOTE: we have a vectors here to allow for multiple threads
    std::vector< std::vector<Vector3> > rp_;
    std::vector< std::vector<Vector3> > a_;
    std::vector< int* > dc_;
};

class CellFinder3D : public CellFinder {
public:
    CellFinder3D( size_t nvertices = 8 );
    
    ~CellFinder3D();

public:
    int find_cell( const Vector3 &p, size_t &ib, size_t &ic, size_t &jc, size_t &kc );
    
    void test_cell( const FV_Cell * cell, const Vector3 &p, int *dc );
    
private:
    // NOTE: we have a vectors here to allow for multiple threads
    std::vector< std::vector<Vector3> > vp_;
    std::vector< std::vector<double> > a_;
    std::vector< int* > dc_;
};

#endif

