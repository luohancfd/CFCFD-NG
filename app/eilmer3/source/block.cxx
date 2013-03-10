/// \file block.cxx
/// \ingroup eilmer3
/// \brief Functions that apply to the block as a whole.
///
/// \author PJ
/// \version 23-Jun-2006 Abstracted from cns_tstp.cxx
///

#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <unistd.h>
extern "C" {
#include <zlib.h>
}
#include "cell.hh"
#include "kernel.hh"
#include "block.hh"
#include "bc.hh"

//-----------------------------------------------------------------------------

// Set some values for debugging a particular cell in a particular block
// These may be defined externally at compile time. Otherwise they are set here.
// If PRNT_DBG_CELL is not set, then nothing is compiled.
#ifndef PRNT_DBG_CELL
#define PRNT_DBG_CELL 0
#endif

#ifndef DBG_BLK
#define DBG_BLK 0
#endif

#ifndef DBG_CELL_I
#define DBG_CELL_I 2
#endif

#ifndef DBG_CELL_J
#define DBG_CELL_J 2
#endif

#ifndef DBG_CELL_K
#define DBG_CELL_K 0
#endif


Block::Block()
    : bcp(NI,NULL) // Let everything else default initialize.
{}

Block::Block(const Block &b)
    : id(b.id), active(b.active), omegaz(b.omegaz),
      dt_allow(b.dt_allow),
      cfl_min(b.cfl_min), cfl_max(b.cfl_max),
      mass_residual(b.mass_residual),
      energy_residual(b.energy_residual),
      mass_residual_loc(b.mass_residual_loc),
      energy_residual_loc(b.energy_residual_loc),
      hncell(b.hncell),
      hicell(b.hicell), hjcell(b.hjcell), hkcell(b.hkcell),
      nidim(b.nidim), njdim(b.njdim), nkdim(b.nkdim),
      nni(b.nni), nnj(b.nnj), nnk(b.nnk),
      imin(b.imin), imax(b.imax),
      jmin(b.jmin), jmax(b.jmax),
      kmin(b.kmin), kmax(b.kmax),
      baldwin_lomax_iturb(b.baldwin_lomax_iturb),
      bcp(b.bcp),
      ctr_(b.ctr_),
      ifi_(b.ifi_), ifj_(b.ifj_), ifk_(b.ifk_),
      vtx_(b.vtx_),
      sifi_(b.sifi_), sifj_(b.sifj_), sifk_(b.sifk_)
{}

Block & Block::operator=(const Block &b)
{
    if ( this != &b ) { // Avoid aliasing
	id = b.id; active = b.active; omegaz = b.omegaz;
	dt_allow = b.dt_allow;
	cfl_min = b.cfl_min; cfl_max = b.cfl_max;
	mass_residual = b.mass_residual;
	energy_residual = b.energy_residual;
	mass_residual_loc = b.mass_residual_loc;
	energy_residual_loc = b.energy_residual_loc;
	hncell = b.hncell;
	hicell = b.hicell; hjcell = b.hjcell; hkcell = b.hkcell;
	nidim = b.nidim; njdim = b.njdim; nkdim = b.nkdim;
	nni = b.nni; nnj = b.nnj; nnk = b.nnk;
	imin = b.imin; imax = b.imax;
	jmin = b.jmin; jmax = b.jmax;
	kmin = b.kmin; kmax = b.kmax;
	baldwin_lomax_iturb = b.baldwin_lomax_iturb;
	bcp = b.bcp;
	ctr_ = b.ctr_;
	ifi_ = b.ifi_; ifj_ = b.ifj_; ifk_ = b.ifk_;
	vtx_ = b.vtx_;
	sifi_ = b.sifi_; sifj_ = b.sifj_; sifk_ = b.sifk_;
    }
    return *this;
}

Block::~Block() 
{
    for ( size_t i = 0; i < bcp.size(); ++i )
	delete bcp[i];
    // for ( size_t i = 0; i < shock_iface_pos.size(); ++i )
    //	delete shock_iface_pos[i];
}

/// \brief Allocate memory for the internal arrays of the block.
/// Returns 0 if successful, 1 otherwise.
int Block::array_alloc(int dimensions)
{
    if ( get_verbose_flag() ) cout << "array_alloc(): Begin for block " <<  id << endl;
    // Check for obvious errors.
    if ( nidim <= 0 || njdim <= 0 || nkdim <= 0 ) {
        cerr << "array_alloc(): Declared dimensions are zero or negative: " 
	     << nidim << ", " << njdim << ", " << nkdim << endl;
	exit( VALUE_ERROR );
    }
    Gas_model *gm = get_gas_model_ptr();

    // Allocate vectors of pointers for the entire block.
    size_t ntot = nidim * njdim * nkdim;
    ctr_.resize(ntot); 
    ifi_.resize(ntot);
    ifj_.resize(ntot);
    if ( dimensions == 3 ) ifk_.resize(ntot);
    vtx_.resize(ntot);
    sifi_.resize(ntot);
    sifj_.resize(ntot);
    if ( dimensions == 3 ) sifk_.resize(ntot);
    // Now, create the actual objects.
    for (size_t ijk = 0; ijk < ntot; ++ijk) {
        ctr_[ijk] = new FV_Cell(gm);
        ifi_[ijk] = new FV_Interface(gm);
        ifj_[ijk] = new FV_Interface(gm);
        if ( dimensions == 3 ) ifk_[ijk] = new FV_Interface(gm);
        vtx_[ijk] = new FV_Vertex(gm);
        sifi_[ijk] = new FV_Interface(gm);
        sifj_[ijk] = new FV_Interface(gm);
        if ( dimensions == 3 ) sifk_[ijk] = new FV_Interface(gm);
    } // ijk loop

    if ( get_verbose_flag() || id == 0 ) {
	cout << "Block " << id << ": finished creating " << ntot << " cells." << endl;
    }
    return SUCCESS;
} // end of array_alloc()


int Block::array_cleanup(int dimensions)
{
    if ( get_verbose_flag() || id == 0 ) {
	cout << "array_cleanup(): Begin for block " <<  id << endl;
    }
    // Need to clean up allocated memory.
    size_t ntot = nidim * njdim * nkdim;
    for (size_t ijk = 0; ijk < ntot; ++ijk) {
	delete ctr_[ijk];
	delete ifi_[ijk];
	delete ifj_[ijk];
	if ( dimensions == 3 ) delete ifk_[ijk];
	delete vtx_[ijk];
	delete sifi_[ijk];
	delete sifj_[ijk];
	if ( dimensions == 3 ) delete sifk_[ijk];
    } // ijk loop
    ctr_.clear(); 
    ifi_.clear();
    ifj_.clear();
    ifk_.clear();
    vtx_.clear();
    sifi_.clear();
    sifj_.clear();
    sifk_.clear();
    return SUCCESS;
} // end of array_cleanup()

//-----------------------------------------------------------------------------

/// \brief Apply the function (one extra double parameter) to all cells.
int Block::apply(FV_Cell_MemberFunction_void f, string failure_message_header)
{
    int result_flag;
    FV_Cell *cellp;
    for ( int k = kmin; k <= kmax; ++k ) {
        for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "About to apply: " << failure_message_header << endl;
			cout << "to cell[" << DBG_CELL_I << "," << DBG_CELL_J << "," << DBG_CELL_K << "]";
			cout << " in block: " << DBG_BLK << endl;
			cout << "Data before call: " << endl;
			cellp->print();
		    }
		}
		#endif
		result_flag = CALL_MEMBER_FN(*cellp,f)();
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "Data after call: " << endl;
			cellp->print();
		    }
		}
		#endif
		if ( result_flag != 0 ) {
		    cout << failure_message_header << endl;
		    printf("Block %d: cell[%d][%d][%d] \n", id, i, j, k);
		    cellp->print();
		    exit( BAD_CELLS_ERROR );
		}
	    }
	}
    }
    return SUCCESS;
} // end of apply(f)

/// \brief Apply the function (one extra double parameter) to all cells.
int Block::apply(FV_Cell_MemberFunction_double f, double param1, string failure_message_header)
{
    int result_flag;
    FV_Cell *cellp;
    for ( int k = kmin; k <= kmax; ++k ) {
        for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "About to apply: " << failure_message_header << endl;
			cout << "to cell[" << DBG_CELL_I << "," << DBG_CELL_J << "," << DBG_CELL_K << "]";
			cout << " in block: " << DBG_BLK << endl;
			cout << "Data before call: " << endl;
			cellp->print();
		    }
		}
		#endif
		result_flag = CALL_MEMBER_FN(*cellp,f)(param1);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "Data after call: " << endl;
			cellp->print();
		    }
		}
		#endif
		if ( result_flag != 0 ) {
		    cout << failure_message_header << endl;
		    printf("Block %d: cell[%d][%d][%d] \n", id, i, j, k);
		    cellp->print();
		    exit( BAD_CELLS_ERROR );
		}
	    }
	}
    }
    return SUCCESS;
} // end of apply(f, p1)

/// \brief Apply the function (two extra double parameter) to all cells.
int Block::apply(FV_Cell_MemberFunction_double_double f, double param1, double param2, 
		 string failure_message_header)
{
    int result_flag;
    FV_Cell *cellp;
    for ( int k = kmin; k <= kmax; ++k ) {
        for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "About to apply: " << failure_message_header << endl;
			cout << "to cell[" << DBG_CELL_I << "," << DBG_CELL_J << "," << DBG_CELL_K << "]";
			cout << " in block: " << DBG_BLK << endl;
			cout << "Data before call: " << endl;
			cellp->print();
		    }
		}
		#endif
		result_flag = CALL_MEMBER_FN(*cellp,f)(param1,param2);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "Data after call: " << endl;
			cellp->print();
		    }
		}
		#endif
		if ( result_flag != 0 ) {
		    cout << failure_message_header << endl;
		    printf("Block %d: cell[%d][%d][%d] \n", id, i, j, k);
		    cellp->print();
		    exit( BAD_CELLS_ERROR );
		}
	    }
	}
    }
    return SUCCESS;
} // end of apply(f, p1, p2)

/// \brief Apply the function (one extra int one extra double parameter) to all cells.
int Block::apply(FV_Cell_MemberFunction_int_double f, int param1, double param2, 
		 string failure_message_header)
{
    int result_flag;
    FV_Cell *cellp;
    for ( int k = kmin; k <= kmax; ++k ) {
        for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "About to apply: " << failure_message_header << endl;
			cout << "to cell[" << DBG_CELL_I << "," << DBG_CELL_J << "," << DBG_CELL_K << "]";
			cout << " in block: " << DBG_BLK << endl;
			cout << "Data before call: " << endl;
			cellp->print();
		    }
		}
		#endif
		result_flag = CALL_MEMBER_FN(*cellp,f)(param1,param2);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "Data after call: " << endl;
			cellp->print();
		    }
		}
		#endif
		if ( result_flag != 0 ) {
		    cout << failure_message_header << endl;
		    printf("Block %d: cell[%d][%d][%d] \n", id, i, j, k);
		    cellp->print();
		    exit( BAD_CELLS_ERROR );
		}
	    }
	}
    }
    return SUCCESS;
} // end of apply(f, p1, p2)

/// \brief Apply the function (one extra int parameter) to all cells.
int Block::apply(FV_Cell_MemberFunction_int f, int param1, string failure_message_header)
{
    int result_flag;
    FV_Cell *cellp;
    for ( int k = kmin; k <= kmax; ++k ) {
        for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "About to apply: " << failure_message_header << endl;
			cout << "to cell[" << DBG_CELL_I << "," << DBG_CELL_J << "," << DBG_CELL_K << "]";
			cout << " in block: " << DBG_BLK << endl;
			cout << "Data before call: " << endl;
			cellp->print();
		    }
		}
		#endif
		result_flag = CALL_MEMBER_FN(*cellp,f)(param1);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "Data after call: " << endl;
			cellp->print();
		    }
		}
		#endif
		if ( result_flag != 0 ) {
		    cout << failure_message_header << endl;
		    printf("Block %d: cell[%d][%d][%d] \n", id, i, j, k);
		    cellp->print();
		    exit( BAD_CELLS_ERROR );
		}
	    }
	}
    }
    return SUCCESS;
} // end of apply(f, p1)

/// \brief Apply the function (two extra int parameters) to all cells.
int Block::apply(FV_Cell_MemberFunction_int_int f, int param1, int param2, string failure_message_header)
{
    int result_flag;
    FV_Cell *cellp;
    for ( int k = kmin; k <= kmax; ++k ) {
        for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "About to apply: " << failure_message_header << endl;
			cout << "to cell[" << DBG_CELL_I << "," << DBG_CELL_J << "," << DBG_CELL_K << "]";
			cout << " in block: " << DBG_BLK << endl;
			cout << "Data before call: " << endl;
			cellp->print();
		    }
		}
		#endif
		result_flag = CALL_MEMBER_FN(*cellp,f)(param1,param2);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "Data after call: " << endl;
			cellp->print();
		    }
		}
		#endif
		if ( result_flag != 0 ) {
		    cout << failure_message_header << endl;
		    printf("Block %d: cell[%d][%d][%d] \n", id, i, j, k);
		    cellp->print();
		    exit( BAD_CELLS_ERROR );
		}
	    }
	}
    }
    return SUCCESS;
} // end of apply(f, p1, p2)

/// \brief Apply the function (no extra parameters) to all cells.
int Block::apply(int (*f)(FV_Cell *cellp), 
		 string failure_message_header)
{
    int result_flag;
    FV_Cell *cellp;
    for ( int k = kmin; k <= kmax; ++k ) {
        for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "About to apply: " << failure_message_header << endl;
			cout << "to cell[" << DBG_CELL_I << "," << DBG_CELL_J << "," << DBG_CELL_K << "]";
			cout << " in block: " << DBG_BLK << endl;
			cout << "Data before call: " << endl;
			cellp->print();
		    }
		}
		#endif
		result_flag = (*f)( cellp );
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "Data after call: " << endl;
			cellp->print();
		    }
		}
		#endif
		if ( result_flag != 0 ) {
		    cout << failure_message_header << endl;
		    printf("Block %d: cell[%d][%d][%d] \n", id, i, j, k);
		    cellp->print();
		    exit( BAD_CELLS_ERROR );
		}
	    }
	}
    }
    return SUCCESS;
} // end of apply(f)

/// \brief Apply the function (one extra double parameter) to all cells.
int Block::apply(int (*f)(FV_Cell *cellp, double param1), 
		 double param1, string failure_message_header)
{
    int result_flag;
    FV_Cell *cellp;
    for ( int k = kmin; k <= kmax; ++k ) {
        for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "About to apply: " << failure_message_header << endl;
			cout << "to cell[" << DBG_CELL_I << "," << DBG_CELL_J << "," << DBG_CELL_K << "]";
			cout << " in block: " << DBG_BLK << endl;
			cout << "Data before call: " << endl;
			cellp->print();
		    }
		}
		#endif
		result_flag = (*f)( cellp, param1 );
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "Data after call: " << endl;
			cellp->print();
		    }
		}
		#endif
		if ( result_flag != 0 ) {
		    cout << failure_message_header << endl;
		    printf("Block %d: cell[%d][%d][%d] \n", id, i, j, k);
		    cellp->print();
		    exit( BAD_CELLS_ERROR );
		}
	    }
	}
    }
    return SUCCESS;
} // end of apply(f, p1)

/// \brief Apply the function (with one extra int parameter) to all cells.
int Block::apply(int (*f)(FV_Cell *cellp, int param1), 
		 int param1, string failure_message_header)
{
    int result_flag;
    FV_Cell *cellp;
    for ( int k = kmin; k <= kmax; ++k ) {
        for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "About to apply: " << failure_message_header << endl;
			cout << "to cell[" << DBG_CELL_I << "," << DBG_CELL_J << "," << DBG_CELL_K << "]";
			cout << " in block: " << DBG_BLK << endl;
			cout << "Data before call: " << endl;
			cellp->print();
		    }
		}
		#endif
		result_flag = (*f)( cellp, param1 );
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "Data after call: " << endl;
			cellp->print();
		    }
		}
		#endif
		if ( result_flag != 0 ) {
		    cout << failure_message_header << endl;
		    printf("Block %d: cell[%d][%d][%d] \n", id, i, j, k);
		    cellp->print();
		    exit( BAD_CELLS_ERROR );
		}
	    }
	}
    }
    return SUCCESS;
} // end of apply(f, p1)

/// \brief Apply the function (with two extra int parameters) to all cells.
int Block::apply(int (*f)(FV_Cell *cellp, int param1, int param2), 
		 int param1, int param2, string failure_message_header)
{
    int result_flag;
    FV_Cell *cellp;
    for ( int k = kmin; k <= kmax; ++k ) {
        for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "About to apply: " << failure_message_header << endl;
			cout << "to cell[" << DBG_CELL_I << "," << DBG_CELL_J << "," << DBG_CELL_K << "]";
			cout << " in block: " << DBG_BLK << endl;
			cout << "Data before call: " << endl;
			cellp->print();
		    }
		}
		#endif
		result_flag = (*f)( cellp, param1, param2 );
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "Data after call: " << endl;
			cellp->print();
		    }
		}
		#endif
		if ( result_flag != 0 ) {
		    cout << failure_message_header << endl;
		    printf("Block %d: cell[%d][%d][%d] \n", id, i, j, k);
		    cellp->print();
		    exit( BAD_CELLS_ERROR );
		}
	    }
	}
    }
    return SUCCESS;
} // end of apply(f, p1, p2)

/// \brief Apply the function (one extra int one extra double parameter) to all cells.
int Block::apply(int (*f)(FV_Cell *cellp, int param1, double param2), 
		 int param1, double param2, string failure_message_header)
{
    int result_flag;
    FV_Cell *cellp;
    for ( int k = kmin; k <= kmax; ++k ) {
        for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "About to apply: " << failure_message_header << endl;
			cout << "to cell[" << DBG_CELL_I << "," << DBG_CELL_J << "," << DBG_CELL_K << "]";
			cout << " in block: " << DBG_BLK << endl;
			cout << "Data before call: " << endl;
			cellp->print();
		    }
		}
		#endif
		result_flag = (*f)( cellp, param1, param2 );
		#if PRNT_DBG_CELL
		if ( id == DBG_BLK) {
		    if ( i == DBG_CELL_I && j == DBG_CELL_J && k == DBG_CELL_K ) {
			cout << "Data after call: " << endl;
			cellp->print();
		    }
		}
		#endif
		if ( result_flag != 0 ) {
		    cout << failure_message_header << endl;
		    printf("Block %d: cell[%d][%d][%d] \n", id, i, j, k);
		    cellp->print();
		    exit( BAD_CELLS_ERROR );
		}
	    }
	}
    }
    return SUCCESS;
} // end of apply(f, p1, p2)


int Block::bind_interfaces_to_cells( int dimensions )
{
    FV_Cell *cellp;
    int kstart, kend;

    if ( dimensions == 3 ) {
	kstart = kmin-1;
	kend = kmax+1;
    } else {
	kstart = 0;
	kend = 0;
    }
    for ( int k = kstart; k <= kend; ++k ) {
	for ( int j = jmin-1; j <= jmax+1; ++j ) {
	    for ( int i = imin-1; i <= imax+1; ++i ) {
		cellp = get_cell(i,j,k);
		cellp->iface[NORTH] = get_ifj(i,j+1,k);
		cellp->iface[EAST] = get_ifi(i+1,j,k);
		cellp->iface[SOUTH] = get_ifj(i,j,k);
		cellp->iface[WEST] = get_ifi(i,j,k);
		cellp->vtx[0] = get_vtx(i,j,k);
		cellp->vtx[1] = get_vtx(i+1,j,k);
		cellp->vtx[2] = get_vtx(i+1,j+1,k);
		cellp->vtx[3] = get_vtx(i,j+1,k);
		if ( dimensions == 3 ) {
		    cellp->iface[TOP] = get_ifk(i,j,k+1);
		    cellp->iface[BOTTOM] = get_ifk(i,j,k);
		    cellp->vtx[4] = get_vtx(i,j,k+1);
		    cellp->vtx[5] = get_vtx(i+1,j,k+1);
		    cellp->vtx[6] = get_vtx(i+1,j+1,k+1);
		    cellp->vtx[7] = get_vtx(i,j+1,k+1);
		} else {
		    cellp->iface[TOP] = NULL;
		    cellp->iface[BOTTOM] = NULL;
		    cellp->vtx[4] = NULL;
		    cellp->vtx[5] = NULL;
		    cellp->vtx[6] = NULL;
		    cellp->vtx[7] = NULL;
		} // end if
	    }
	}
    }
    return SUCCESS;
} // end bind_interfaces_to_cells()


/// \brief Set the base heat source values for this block.
int Block::set_base_qdot( global_data &gd )
{
    double total_qdot_for_block = 0.0;
    CHeatZone *hzp;
    FV_Cell *cellp;

    total_qdot_for_block = 0.0;
    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int i = imin; i <= imax; ++i ) {
	    for ( int j = jmin; j <= jmax; ++j ) {
		cellp = get_cell(i,j,k);
		cellp->base_qdot = 0.0;
		for ( int indx = 0; indx < gd.n_heat_zone; ++indx ) {
		    hzp = &(gd.heat_zone[indx]);
		    if ( cellp->pos.x >= hzp->x0 && cellp->pos.x <= hzp->x1 &&
			 cellp->pos.y >= hzp->y0 && cellp->pos.y <= hzp->y1 &&
			 (gd.dimensions == 2 || 
			  (cellp->pos.z >= hzp->z0 && cellp->pos.z <= hzp->z1)) ) {
			cellp->base_qdot += hzp->qdot;
		    }
		} // for indx
		total_qdot_for_block += cellp->base_qdot * cellp->volume;
	    } // for j
	} // for i
    } // for k
    if ( total_qdot_for_block > 1.0e-10 ) {
	cout << "set_base_qdot(): block " << id
	     << " base_qdot= " << total_qdot_for_block << " Watts" << endl;
    }
    return SUCCESS;
} // end set_base_qdot()


/// \brief Set the reactions-allowed flag for cells in this block.
int Block::identify_reaction_zones( global_data &gd )
{
    int total_cells_in_reaction_zones = 0;
    int total_cells = 0;
    CReactionZone *rzp;
    FV_Cell *cellp;

    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int i = imin; i <= imax; ++i ) {
	    for ( int j = jmin; j <= jmax; ++j ) {
		cellp = get_cell(i,j,k);
		if ( gd.n_reaction_zone > 0 ) {
		    // User-specified reaction zones; mask off reacting/nonreacting zones.
		    cellp->fr_reactions_allowed = 0;
		    for ( int indx = 0; indx < gd.n_reaction_zone; ++indx ) {
			rzp = &(gd.reaction_zone[indx]);
			if ( cellp->pos.x >= rzp->x0 && cellp->pos.x <= rzp->x1 &&
			     cellp->pos.y >= rzp->y0 && cellp->pos.y <= rzp->y1 &&
			     (gd.dimensions == 2 || 
			      (cellp->pos.z >= rzp->z0 && cellp->pos.z <= rzp->z1)) ) {
			    cellp->fr_reactions_allowed = 1;
			}
		    } // end for( indx
		} else {
		    // No user-specified zones; always allow reactions.
		    cellp->fr_reactions_allowed = 1;
		}
		total_cells_in_reaction_zones += cellp->fr_reactions_allowed;
		total_cells += 1;
	    } // for j
	} // for i
    } // for k
    if ( get_reacting_flag() ) {
	cout << "identify_reaction_zones(): block " << id
	     << " cells inside zones = " << total_cells_in_reaction_zones 
	     << " out of " << total_cells << endl;
	if ( gd.n_reaction_zone == 0 ) {
	    cout << "Note that for no user-specified zones,"
		 << " the whole domain is allowed to be reacting." << endl;
	}
    }
    return SUCCESS;
} // end identify_reaction_zones()


/// \brief Set the in-turbulent-zone flag for cells in this block.
int Block::identify_turbulent_zones( global_data &gd )
{
    int total_cells_in_turbulent_zones = 0;
    int total_cells = 0;
    CTurbulentZone *tzp;
    FV_Cell *cellp;
    
    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int i = imin; i <= imax; ++i ) {
	    for ( int j = jmin; j <= jmax; ++j ) {
		cellp = get_cell(i,j,k);
		if ( gd.n_turbulent_zone > 0 ) {
		    cellp->in_turbulent_zone = 0;
		    for ( int indx = 0; indx < gd.n_turbulent_zone; ++indx ) {
			tzp = &(gd.turbulent_zone[indx]);
			if ( cellp->pos.x >= tzp->x0 && cellp->pos.x <= tzp->x1 &&
			     cellp->pos.y >= tzp->y0 && cellp->pos.y <= tzp->y1 &&
			     (gd.dimensions == 2 || 
			      (cellp->pos.z >= tzp->z0 && cellp->pos.z <= tzp->z1)) ) {
			    cellp->in_turbulent_zone = 1;
			}
		    } // for indx
		} else {
		    cellp->in_turbulent_zone = 1;
		}
		total_cells_in_turbulent_zones += cellp->in_turbulent_zone;
		total_cells += 1;
	    } // for j
	} // for i
    } // for k
    if ( get_turbulence_flag() ) {
	cout << "identify_turbulent_zones(): block " << id
	     << " cells inside zones = " << total_cells_in_turbulent_zones 
	     << " out of " << total_cells << endl;
	if ( gd.n_turbulent_zone == 0 ) {
	    cout << "Note that for no user-specified zones,"
		 << " the whole domain is allowed to be turbulent." << endl;
	}
    }
    return SUCCESS;
} // end identify_turbulent_zones()


int Block::clear_fluxes_of_conserved_quantities( int dimensions )
{
    FV_Interface *IFace;

    for ( int k = kmin; k <= kmax; ++k ) {
	for (int j = jmin; j <= jmax; ++j) {
	    for (int i = imin; i <= imax+1; ++i) {
		IFace = get_ifi(i,j,k);
		IFace->F->clear_values();
	    } // for i
	} // for j
    } // for k
    for ( int k = kmin; k <= kmax; ++k ) {
	for (int j = jmin; j <= jmax+1; ++j) {
	    for (int i = imin; i <= imax; ++i) {
		IFace = get_ifj(i,j,k);
		IFace->F->clear_values();
	    } // for i
	} // for j
    } // for k
    if ( dimensions == 3 ) {
	for ( int k = kmin; k <= kmax+1; ++k ) {
	    for (int j = jmin; j <= jmax; ++j) {
		for (int i = imin; i <= imax; ++i) {
		    IFace = get_ifk(i,j,k);
		    IFace->F->clear_values();
		} // for i
	    } // for j
	} // for k
    } // end if G.dimensions == 3
    return SUCCESS;
}

int Block::propagate_data_west_to_east( int dimensions )
// Propagate data from the west ghost cell, right across the block.
// This is a useful starting state for the block-sequenced calculation
// where the final flow is expected to be steady-state.
{
    FV_Cell *src, *dest;
    Gas_model *gm = get_gas_model_ptr();
    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int j = jmin; j <= jmax; ++j) {
	    src = get_cell(imin-1,j);
	    for ( int i = imin; i <= imax; ++i ) {
		dest = get_cell(i,j,k);
		dest->copy_values_from(*src, COPY_FLOW_STATE);
		Gas_data *gas = dest->fs->gas;
		if ( gm->eval_thermo_state_pT(*gas) != SUCCESS ||
		     gm->eval_transport_coefficients(*gas) != SUCCESS ) {
		    printf( "propagate_data_west_to_east(): Duff call to thermo model.\n" );
		    printf( "   i=%d, j=%d, k=%d\n", i, j, k );
		    gas->print_values();
		    exit(DUFF_EOS_ERROR);  /* Might as well quit early. */
		}
	    } // for i
	} // for j
    } // for k
    return SUCCESS;
} // end propagate_data_west_to_east()

int Block::compute_initial_primary_cell_geometric_data( int dimensions )
// Compute cell and interface geometric properties.
{
    int i, j, k;
    FV_Cell *cell, *cell_1, *cell_2, *ghost_cell;
    Vector3 dummy;
    FV_Interface *iface;
    Vector3 *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7;
    Vector3 ds;

    if ( dimensions == 2 ) {
	calc_initial_volumes_2D();
	calc_initial_faces_2D();
	calc_initial_ghost_cell_geom_2D();
	return SUCCESS;
    }

    // Cell properties of volume and position.
    // Estimates of cross-cell distances for use in high-order reconstruction.
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		cell = get_cell(i,j,k);
		p0 = &(get_vtx(i,j,k)->pos);
		p1 = &(get_vtx(i+1,j,k)->pos);
		p2 = &(get_vtx(i+1,j+1,k)->pos);
		p3 = &(get_vtx(i,j+1,k)->pos);
		p4 = &(get_vtx(i,j,k+1)->pos);
		p5 = &(get_vtx(i+1,j,k+1)->pos);
		p6 = &(get_vtx(i+1,j+1,k+1)->pos);
		p7 = &(get_vtx(i,j+1,k+1)->pos);
		hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				     cell->pos, cell->volume, cell->iLength,
				     cell->jLength, cell->kLength);
		cell->L_min = cell->iLength;
		if ( cell->jLength < cell->L_min ) cell->L_min = cell->jLength;
		if ( cell->kLength < cell->L_min ) cell->L_min = cell->kLength;
	    }
	}
    }

    // work on ifi face as a WEST face
    // t1 in the j-ordinate direction
    // t2 in the k-ordinate direction
    for ( i = imin; i <= imax + 1; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		iface = get_ifi(i,j,k);
		p0 = &(get_vtx(i,j,k)->pos);
		p3 = &(get_vtx(i,j+1,k)->pos);
		p7 = &(get_vtx(i,j+1,k+1)->pos);
		p4 = &(get_vtx(i,j,k+1)->pos);
		quad_properties( *p0, *p3, *p7, *p4,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area );
	    }
	}
    }

    // work on ifj face as a SOUTH face
    // t1 in the k-ordinate direction
    // t2 in the i-ordinate direction
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax + 1; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		iface = get_ifj(i,j,k);
		p0 = &(get_vtx(i,j,k)->pos);
		p4 = &(get_vtx(i,j,k+1)->pos);
		p5 = &(get_vtx(i+1,j,k+1)->pos);
		p1 = &(get_vtx(i+1,j,k)->pos);
		quad_properties( *p0, *p4, *p5, *p1,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area );
	    }
	}
    }

    // work on ifk face as a BOTTOM face
    // t1 in the i-ordinate direction
    // t2 in the j-ordinate direction
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax + 1; ++k ) {
		iface = get_ifk(i,j,k);
		p0 = &(get_vtx(i,j,k)->pos);
		p1 = &(get_vtx(i+1,j,k)->pos);
		p2 = &(get_vtx(i+1,j+1,k)->pos);
		p3 = &(get_vtx(i,j+1,k)->pos);
		quad_properties( *p0, *p1, *p2, *p3,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area );
	    }
	}
    }

    /* Propagate cross-cell lengths into the ghost cells. */
    for ( j = jmin; j <= jmax; ++j ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    i = imin;
	    cell = get_cell(i,j,k);
	    get_cell(i-1,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i-2,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    i = imax;
	    cell = get_cell(i,j,k);
	    get_cell(i+1,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i+2,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    j = jmin;
	    cell = get_cell(i,j,k);
	    get_cell(i,j-1,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i,j-2,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    j = jmax;
	    cell = get_cell(i,j,k);
	    get_cell(i,j+1,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i,j+2,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    k = kmin;
	    cell = get_cell(i,j,k);
	    get_cell(i,j,k-1)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i,j,k-2)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    k = kmax;
	    cell = get_cell(i,j,k);
	    get_cell(i,j,k+1)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i,j,k+2)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	}
    }

    /* Extrapolate (with first-order) cell positions and volumes to ghost cells. */
    for ( j = jmin; j <= jmax; ++j ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    i = imin;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i+1,j,k);
	    ghost_cell = get_cell(i-1,j,k);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i-2,j,k);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_2->volume - cell_2->volume;
	    i = imax;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i-1,j,k);
	    ghost_cell = get_cell(i+1,j,k);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i+2,j,k);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    j = jmin;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j+1,k);
	    ghost_cell = get_cell(i,j-1,k);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j-2,k);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	    j = jmax;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j-1,k);
	    ghost_cell = get_cell(i,j+1,k);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j+2,k);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    k = kmin;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j,k+1);
	    ghost_cell = get_cell(i,j,k-1);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j,k-2);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	    k = kmax;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j,k-1);
	    ghost_cell = get_cell(i,j,k+1);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j,k+2);
	    ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	    ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	}
    }
    
    return SUCCESS;
} // end compute_initial_primary_cell_geometric_data()


int Block::compute_primary_cell_geometric_data( int dimensions, int time_level )
// Compute cell and interface geometric properties.
{
    int i, j, k;
    FV_Cell *cell, *cell_1, *cell_2, *ghost_cell;
    Vector3 dummy;
    FV_Interface *iface;
    Vector3 *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7;
    Vector3 ds;

    if ( dimensions == 2 ) {
	calc_volumes_2D(time_level);
	calc_faces_2D(time_level);
	calc_ghost_cell_geom_2D(time_level);
	return SUCCESS;
    }

    // Cell properties of volume and position.
    // Estimates of cross-cell distances for use in high-order reconstruction.
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		cell = get_cell(i,j,k);
		p0 = &(get_vtx(i,j,k)->position[time_level]);
		p1 = &(get_vtx(i+1,j,k)->position[time_level]);
		p2 = &(get_vtx(i+1,j+1,k)->position[time_level]);
		p3 = &(get_vtx(i,j+1,k)->position[time_level]);
		p4 = &(get_vtx(i,j,k+1)->position[time_level]);
		p5 = &(get_vtx(i+1,j,k+1)->position[time_level]);
		p6 = &(get_vtx(i+1,j+1,k+1)->position[time_level]);
		p7 = &(get_vtx(i,j+1,k+1)->position[time_level]);
		hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				     cell->position[time_level], cell->vol[time_level], cell->iLength,
				     cell->jLength, cell->kLength);
		cell->L_min = cell->iLength;
		if ( cell->jLength < cell->L_min ) cell->L_min = cell->jLength;
		if ( cell->kLength < cell->L_min ) cell->L_min = cell->kLength;
	    }
	}
    }

    // work on ifi face as a WEST face
    // t1 in the j-ordinate direction
    // t2 in the k-ordinate direction
    for ( i = imin; i <= imax + 1; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		iface = get_ifi(i,j,k);
		p0 = &(get_vtx(i,j,k)->position[time_level]);
		p3 = &(get_vtx(i,j+1,k)->position[time_level]);
		p7 = &(get_vtx(i,j+1,k+1)->position[time_level]);
		p4 = &(get_vtx(i,j,k+1)->position[time_level]);
		quad_properties( *p0, *p3, *p7, *p4,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->ar[time_level] );
	    }
	}
    }

    // work on ifj face as a SOUTH face
    // t1 in the k-ordinate direction
    // t2 in the i-ordinate direction
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax + 1; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		iface = get_ifj(i,j,k);
		p0 = &(get_vtx(i,j,k)->position[time_level]);
		p4 = &(get_vtx(i,j,k+1)->position[time_level]);
		p5 = &(get_vtx(i+1,j,k+1)->position[time_level]);
		p1 = &(get_vtx(i+1,j,k)->position[time_level]);
		quad_properties( *p0, *p4, *p5, *p1,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->ar[time_level] );
	    }
	}
    }

    // work on ifk face as a BOTTOM face
    // t1 in the i-ordinate direction
    // t2 in the j-ordinate direction
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax + 1; ++k ) {
		iface = get_ifk(i,j,k);
		p0 = &(get_vtx(i,j,k)->position[time_level]);
		p1 = &(get_vtx(i+1,j,k)->position[time_level]);
		p2 = &(get_vtx(i+1,j+1,k)->position[time_level]);
		p3 = &(get_vtx(i,j+1,k)->position[time_level]);
		quad_properties( *p0, *p1, *p2, *p3,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->ar[time_level] );
	    }
	}
    }

    /* Propagate cross-cell lengths into the ghost cells. */
    for ( j = jmin; j <= jmax; ++j ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    i = imin;
	    cell = get_cell(i,j,k);
	    get_cell(i-1,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i-2,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    i = imax;
	    cell = get_cell(i,j,k);
	    get_cell(i+1,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i+2,j,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    j = jmin;
	    cell = get_cell(i,j,k);
	    get_cell(i,j-1,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i,j-2,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    j = jmax;
	    cell = get_cell(i,j,k);
	    get_cell(i,j+1,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i,j+2,k)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    k = kmin;
	    cell = get_cell(i,j,k);
	    get_cell(i,j,k-1)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i,j,k-2)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    k = kmax;
	    cell = get_cell(i,j,k);
	    get_cell(i,j,k+1)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	    get_cell(i,j,k+2)->copy_values_from(*cell, COPY_CELL_LENGTHS);
	}
    }

    /* Extrapolate (with first-order) cell positions and volumes to ghost cells. */
    for ( j = jmin; j <= jmax; ++j ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    i = imin;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i+1,j,k);
	    ghost_cell = get_cell(i-1,j,k);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i-2,j,k);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_2->vol[time_level] - cell_2->vol[time_level];
	    i = imax;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i-1,j,k);
	    ghost_cell = get_cell(i+1,j,k);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i+2,j,k);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    j = jmin;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j+1,k);
	    ghost_cell = get_cell(i,j-1,k);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j-2,k);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	    j = jmax;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j-1,k);
	    ghost_cell = get_cell(i,j+1,k);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j+2,k);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    k = kmin;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j,k+1);
	    ghost_cell = get_cell(i,j,k-1);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j,k-2);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	    k = kmax;
	    cell_1 = get_cell(i,j,k);
	    cell_2 = get_cell(i,j,k-1);
	    ghost_cell = get_cell(i,j,k+1);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	    cell_2 = cell_1;
	    cell_1 = ghost_cell;
	    ghost_cell = get_cell(i,j,k+2);
	    ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	    ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	}
    }
    
    return SUCCESS;
} // end compute_primary_cell_geometric_data()

int Block::compute_distance_to_nearest_wall_for_all_cells( int dimensions )
{
    int i, j, k;
    FV_Cell *cell, *cell_at_wall[6];
    double dx, dy, dz, dist[6], half_width[6];
    FV_Interface *face_at_wall;

    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		cell = get_cell(i,j,k);
		// Step 1: get distances to all boundaries along all index directions.
		// If the block is not too distorted, these directions should take us
		// straight to the bounding walls.
		// North
		face_at_wall = get_ifj(i,jmax+1,k);
	        dx = cell->pos.x - face_at_wall->pos.x;
	        dy = cell->pos.y - face_at_wall->pos.y;
	        dz = cell->pos.z - face_at_wall->pos.z;
	        dist[NORTH] = sqrt(dx*dx + dy*dy + dz*dz);
		cell_at_wall[NORTH] = get_cell(i,jmax,k);
	        dx = cell_at_wall[NORTH]->pos.x - face_at_wall->pos.x;
	        dy = cell_at_wall[NORTH]->pos.y - face_at_wall->pos.y;
	        dz = cell_at_wall[NORTH]->pos.z - face_at_wall->pos.z;
	        half_width[NORTH] = sqrt(dx*dx + dy*dy + dz*dz);
		// East
		face_at_wall = get_ifi(imax+1,j,k);
	        dx = cell->pos.x - face_at_wall->pos.x;
	        dy = cell->pos.y - face_at_wall->pos.y;
	        dz = cell->pos.z - face_at_wall->pos.z;
	        dist[EAST] = sqrt(dx*dx + dy*dy + dz*dz);
		cell_at_wall[EAST] = get_cell(imax,j,k);
	        dx = cell_at_wall[EAST]->pos.x - face_at_wall->pos.x;
	        dy = cell_at_wall[EAST]->pos.y - face_at_wall->pos.y;
	        dz = cell_at_wall[EAST]->pos.z - face_at_wall->pos.z;
	        half_width[EAST] = sqrt(dx*dx + dy*dy + dz*dz);
		// South
		face_at_wall = get_ifj(i,jmin,k);
	        dx = cell->pos.x - face_at_wall->pos.x;
	        dy = cell->pos.y - face_at_wall->pos.y;
	        dz = cell->pos.z - face_at_wall->pos.z;
	        dist[SOUTH] = sqrt(dx*dx + dy*dy + dz*dz);
		cell_at_wall[SOUTH] = get_cell(i,jmin,k);
	        dx = cell_at_wall[SOUTH]->pos.x - face_at_wall->pos.x;
	        dy = cell_at_wall[SOUTH]->pos.y - face_at_wall->pos.y;
	        dz = cell_at_wall[SOUTH]->pos.z - face_at_wall->pos.z;
	        half_width[SOUTH] = sqrt(dx*dx + dy*dy + dz*dz);
		// West
		face_at_wall = get_ifi(imin,j,k);
	        dx = cell->pos.x - face_at_wall->pos.x;
	        dy = cell->pos.y - face_at_wall->pos.y;
	        dz = cell->pos.z - face_at_wall->pos.z;
	        dist[WEST] = sqrt(dx*dx + dy*dy + dz*dz);
		cell_at_wall[WEST] = get_cell(imin,j,k);
	        dx = cell_at_wall[WEST]->pos.x - face_at_wall->pos.x;
	        dy = cell_at_wall[WEST]->pos.y - face_at_wall->pos.y;
	        dz = cell_at_wall[WEST]->pos.z - face_at_wall->pos.z;
	        half_width[WEST] = sqrt(dx*dx + dy*dy + dz*dz);
		if ( dimensions == 3 ) {
		    // Top
		    face_at_wall = get_ifk(i,j,kmax+1);
		    dx = cell->pos.x - face_at_wall->pos.x;
		    dy = cell->pos.y - face_at_wall->pos.y;
		    dz = cell->pos.z - face_at_wall->pos.z;
		    dist[TOP] = sqrt(dx*dx + dy*dy + dz*dz);
		    cell_at_wall[TOP] = get_cell(i,j,kmax);
		    dx = cell_at_wall[TOP]->pos.x - face_at_wall->pos.x;
		    dy = cell_at_wall[TOP]->pos.y - face_at_wall->pos.y;
		    dz = cell_at_wall[TOP]->pos.z - face_at_wall->pos.z;
		    half_width[TOP] = sqrt(dx*dx + dy*dy + dz*dz);
		    // Bottom
		    face_at_wall = get_ifk(i,j,kmin);
		    dx = cell->pos.x - face_at_wall->pos.x;
		    dy = cell->pos.y - face_at_wall->pos.y;
		    dz = cell->pos.z - face_at_wall->pos.z;
		    dist[BOTTOM] = sqrt(dx*dx + dy*dy + dz*dz);
		    cell_at_wall[BOTTOM] = get_cell(i,j,kmin);
		    dx = cell_at_wall[BOTTOM]->pos.x - face_at_wall->pos.x;
		    dy = cell_at_wall[BOTTOM]->pos.y - face_at_wall->pos.y;
		    dz = cell_at_wall[BOTTOM]->pos.z - face_at_wall->pos.z;
		    half_width[BOTTOM] = sqrt(dx*dx + dy*dy + dz*dz);
		} else {
		    cell_at_wall[TOP] = NULL;
		    dist[TOP] = 0.0;
		    half_width[TOP] = 0.0;
		    cell_at_wall[BOTTOM] = NULL;
		    dist[BOTTOM] = 0.0;
		    half_width[BOTTOM] = 0.0;
		}

		// Step 2: Just in case there are no real walls for this block...
		//
		// We'll start by storing the largest distance and 
		// corresponding wall-cell half-width so that we have 
		// a relatively large distance in case there are no walls
		// on the boundary of the block.
		cell->distance_to_nearest_wall = dist[0];
		cell->half_cell_width_at_wall = half_width[0];
		for ( int iface = 1; iface < 6; ++iface ) {
		    if ( cell_at_wall[iface] != NULL &&
			 dist[iface] > cell->distance_to_nearest_wall ) {
			cell->distance_to_nearest_wall = dist[iface];
			cell->half_cell_width_at_wall = half_width[iface];
			cell->cell_at_nearest_wall = cell_at_wall[iface];
		    }
		}
		
		// Step 3: find the closest real wall.
		for ( int iface = 0; iface < 6; ++iface ) {
		    if ( dimensions == 2 && iface >= 4 ) break; // Top,Bottom are 4,5
		    if ( bcp[iface]->is_wall() && 
			 bcp[iface]->type_code != SLIP_WALL &&
			 dist[iface] < cell->distance_to_nearest_wall && 
			 cell_at_wall[iface] != NULL ) {
			cell->distance_to_nearest_wall = dist[iface];
			cell->half_cell_width_at_wall = half_width[iface];
			cell->cell_at_nearest_wall = cell_at_wall[iface];
		    }
		}
	    } // for k
	} // for j
    } // for i
    
    return SUCCESS;
} // end compute_distance_to_nearest_wall_for_all_cells()

int Block::compute_secondary_cell_geometric_data( int dimensions )
// Compute secondary-cell and interface geometric properties.
// Will be used for computing gradients for viscous terms.
//
// To do: The code for the 3D cells has been ported from eilmer2 without
//        taking advantage of the eilmer3 structure where ghost-cells have
//        been filled in with useful geometric information.
//        We should make use of this information.
{
    int i, j, k;
    Vector3 dummy;
    double iLen, jLen, kLen;
    FV_Vertex *vertex;
    FV_Interface *iface;
    Vector3 *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7;
    Vector3 ds;

    if ( dimensions == 2 ) {
	secondary_areas_2D();
	return SUCCESS;
    }

    /*
     * Internal secondary cell geometry information
     */
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    for ( k = kmin; k <= kmax-1; ++k ) {
		vertex = get_vtx(i+1,j+1,k+1);
		p0 = &(get_cell(i,j,k)->pos);
		p1 = &(get_cell(i+1,j,k)->pos);
		p2 = &(get_cell(i+1,j+1,k)->pos);
		p3 = &(get_cell(i,j+1,k)->pos);
		p4 = &(get_cell(i,j,k+1)->pos);
		p5 = &(get_cell(i+1,j,k+1)->pos);
		p6 = &(get_cell(i+1,j+1,k+1)->pos);
		p7 = &(get_cell(i,j+1,k+1)->pos);
		hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				     dummy, vertex->volume, iLen, jLen, kLen );
	    }
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    for ( k = kmin; k <= kmax-1; ++k ) {
		iface = get_sifi(i,j,k);
		p0 = &(get_cell(i,j,k)->pos);
		p3 = &(get_cell(i,j+1,k)->pos);
		p7 = &(get_cell(i,j+1,k+1)->pos);
		p4 = &(get_cell(i,j,k+1)->pos);
		quad_properties( *p0, *p3, *p7, *p4,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area );
	    }
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( k = kmin; k <= kmax-1; ++k ) {
		iface = get_sifj(i,j,k);
		p0 = &(get_cell(i,j,k)->pos);
		p4 = &(get_cell(i,j,k+1)->pos);
		p5 = &(get_cell(i+1,j,k+1)->pos);
		p1 = &(get_cell(i+1,j,k)->pos);
		quad_properties( *p0, *p4, *p5, *p1,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area );
	    }
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    for ( k = kmin; k <= kmax; ++k ) {
		iface = get_sifk(i,j,k);
		p0 = &(get_cell(i,j,k)->pos);
		p1 = &(get_cell(i+1,j,k)->pos);
		p2 = &(get_cell(i+1,j+1,k)->pos);
		p3 = &(get_cell(i,j+1,k)->pos);
		quad_properties( *p0, *p1, *p2, *p3,
				 iface->pos, iface->n, iface->t1, iface->t2,
				 iface->area );
	    }
	}
    }

    /*
     * East boundary secondary cell geometry information
     */
    i = imax;
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_cell(i,j,k)->pos);
	    p1 = &(get_ifi(i+1,j,k)->pos);
	    p2 = &(get_ifi(i+1,j+1,k)->pos);
	    p3 = &(get_cell(i,j+1,k)->pos);
	    p4 = &(get_cell(i,j,k+1)->pos);
	    p5 = &(get_ifi(i+1,j,k+1)->pos);
	    p6 = &(get_ifi(i+1,j+1,k+1)->pos);
	    p7 = &(get_cell(i,j+1,k+1)->pos);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifi(i+1,j,k);
	    p1 = &(get_ifi(i+1,j,k)->pos);
	    p2 = &(get_ifi(i+1,j+1,k)->pos);
	    p6 = &(get_ifi(i+1,j+1,k+1)->pos);
	    p5 = &(get_ifi(i+1,j,k+1)->pos);
	    quad_properties( *p1, *p2, *p6, *p5,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( j = jmin; j <= jmax; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifj(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos);
	    p4 = &(get_cell(i,j,k+1)->pos);
	    p5 = &(get_ifi(i+1,j,k+1)->pos);
	    p1 = &(get_ifi(i+1,j,k)->pos);
	    quad_properties( *p0, *p4, *p5, *p1,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    iface = get_sifk(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos);
	    p1 = &(get_ifi(i+1,j,k)->pos);
	    p2 = &(get_ifi(i+1,j+1,k)->pos);
	    p3 = &(get_cell(i,j+1,k)->pos);
	    quad_properties( *p0, *p1, *p2, *p3,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }

    /*
     * West boundary secondary cell geometry information
     */
    i = imin - 1;
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_ifi(i+1,j,k)->pos);
	    p1 = &(get_cell(i+1,j,k)->pos);
	    p2 = &(get_cell(i+1,j+1,k)->pos);
	    p3 = &(get_ifi(i+1,j+1,k)->pos);
	    p4 = &(get_ifi(i+1,j,k+1)->pos);
	    p5 = &(get_cell(i+1,j,k+1)->pos);
	    p6 = &(get_cell(i+1,j+1,k+1)->pos);
	    p7 = &(get_ifi(i+1,j+1,k+1)->pos);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifi(i,j,k);
	    p0 = &(get_ifi(i+1,j,k)->pos);
	    p3 = &(get_ifi(i+1,j+1,k)->pos);
	    p7 = &(get_ifi(i+1,j+1,k+1)->pos);
	    p4 = &(get_ifi(i+1,j,k+1)->pos);
	    quad_properties( *p0, *p3, *p7, *p4,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( j = jmin; j <= jmax; ++j ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifj(i,j,k);
	    p0 = &(get_ifi(i+1,j,k)->pos);
	    p4 = &(get_ifi(i+1,j,k+1)->pos);
	    p5 = &(get_cell(i+1,j,k+1)->pos);
	    p1 = &(get_cell(i+1,j,k)->pos);
	    quad_properties( *p0, *p4, *p5, *p1,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( j = jmin; j <= jmax-1; ++j ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    iface = get_sifk(i,j,k);
	    p0 = &(get_ifi(i+1,j,k)->pos);
	    p1 = &(get_cell(i+1,j,k)->pos);
	    p2 = &(get_cell(i+1,j+1,k)->pos);
	    p3 = &(get_ifi(i+1,j+1,k)->pos);
	    quad_properties( *p0, *p1, *p2, *p3,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }

    /*
     * North boundary secondary cell geometry information
     */
    j = jmax;
    for ( i = imin; i <= imax-1; ++i ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_cell(i,j,k)->pos);
	    p1 = &(get_cell(i+1,j,k)->pos);
	    p2 = &(get_ifj(i+1,j+1,k)->pos);
	    p3 = &(get_ifj(i,j+1,k)->pos);
	    p4 = &(get_cell(i,j,k+1)->pos);
	    p5 = &(get_cell(i+1,j,k+1)->pos);
	    p6 = &(get_ifj(i+1,j+1,k+1)->pos);
	    p7 = &(get_ifj(i,j+1,k+1)->pos);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for (k = kmin; k <= kmax-1; ++k) {
	    iface = get_sifi(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos);
	    p3 = &(get_ifj(i,j+1,k)->pos);
	    p7 = &(get_ifj(i,j+1,k+1)->pos);
	    p4 = &(get_cell(i,j,k+1)->pos);
	    quad_properties( *p0, *p3, *p7, *p4,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for (k = kmin; k <= kmax-1; ++k) {
	    iface = get_sifj(i,j+1,k);
	    p3 = &(get_ifj(i,j+1,k)->pos);
	    p7 = &(get_ifj(i,j+1,k+1)->pos);
	    p6 = &(get_ifj(i+1,j+1,k+1)->pos);
	    p2 = &(get_ifj(i+1,j+1,k)->pos);
	    quad_properties( *p3, *p7, *p6, *p2,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    iface = get_sifk(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos);
	    p1 = &(get_cell(i+1,j,k)->pos);
	    p2 = &(get_ifj(i+1,j+1,k)->pos);
	    p3 = &(get_ifj(i,j+1,k)->pos);
	    quad_properties( *p0, *p1, *p2, *p3,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }

    /*
     * South boundary secondary cell geometry information
     */
    j = jmin - 1;
    for ( i = imin; i <= imax-1; ++i ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_ifj(i,j+1,k)->pos);
	    p1 = &(get_ifj(i+1,j+1,k)->pos);
	    p2 = &(get_cell(i+1,j+1,k)->pos);
	    p3 = &(get_cell(i,j+1,k)->pos);
	    p4 = &(get_ifj(i,j+1,k+1)->pos);
	    p5 = &(get_ifj(i+1,j+1,k+1)->pos);
	    p6 = &(get_cell(i+1,j+1,k+1)->pos);
	    p7 = &(get_cell(i,j+1,k+1)->pos);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifi(i,j,k);
	    p0 = &(get_ifj(i,j+1,k)->pos);
	    p3 = &(get_cell(i,j+1,k)->pos);
	    p7 = &(get_cell(i,j+1,k+1)->pos);
	    p4 = &(get_ifj(i,j+1,k+1)->pos);
	    quad_properties( *p0, *p3, *p7, *p4,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( k = kmin; k <= kmax-1; ++k ) {
	    iface = get_sifj(i,j,k);
	    p0 = &(get_ifj(i,j+1,k)->pos);
	    p4 = &(get_ifj(i,j+1,k+1)->pos);
	    p5 = &(get_ifj(i+1,j+1,k+1)->pos);
	    p1 = &(get_ifj(i+1,j+1,k)->pos);
	    quad_properties( *p0, *p4, *p5, *p1,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( k = kmin; k <= kmax; ++k ) {
	    iface = get_sifk(i,j,k);
	    p0 = &(get_ifj(i,j+1,k)->pos);
	    p1 = &(get_ifj(i+1,j+1,k)->pos);
	    p2 = &(get_cell(i+1,j+1,k)->pos);
	    p3 = &(get_cell(i,j+1,k)->pos);
	    quad_properties( *p0, *p1, *p2, *p3,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }

    /*
     * Top boundary secondary cell geometry information
     */
    k = kmax;
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_cell(i,j,k)->pos);
	    p1 = &(get_cell(i+1,j,k)->pos);
	    p2 = &(get_cell(i+1,j+1,k)->pos);
	    p3 = &(get_cell(i,j+1,k)->pos);
	    p4 = &(get_ifk(i,j,k+1)->pos);
	    p5 = &(get_ifk(i+1,j,k+1)->pos);
	    p6 = &(get_ifk(i+1,j+1,k+1)->pos);
	    p7 = &(get_ifk(i,j+1,k+1)->pos);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    iface = get_sifi(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos);
	    p3 = &(get_cell(i,j+1,k)->pos);
	    p7 = &(get_ifk(i,j+1,k+1)->pos);
	    p4 = &(get_ifk(i,j,k+1)->pos);
	    quad_properties( *p0, *p3, *p7, *p4,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    iface = get_sifj(i,j,k);
	    p0 = &(get_cell(i,j,k)->pos);
	    p4 = &(get_ifk(i,j,k+1)->pos);
	    p5 = &(get_ifk(i+1,j,k+1)->pos);
	    p1 = &(get_cell(i+1,j,k)->pos);
	    quad_properties( *p0, *p4, *p5, *p1,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    iface = get_sifk(i,j,k+1);
	    p4 = &(get_ifk(i,j,k+1)->pos);
	    p5 = &(get_ifk(i+1,j,k+1)->pos);
	    p6 = &(get_ifk(i+1,j+1,k+1)->pos);
	    p7 = &(get_ifk(i,j+1,k+1)->pos);
	    quad_properties( *p4, *p5, *p6, *p7,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }

    /*
     * Bottom boundary secondary cell geometry information
     */
    k = kmin - 1;
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    vertex = get_vtx(i+1,j+1,k+1);
	    p0 = &(get_ifk(i,j,k+1)->pos);
	    p1 = &(get_ifk(i+1,j,k+1)->pos);
	    p2 = &(get_ifk(i+1,j+1,k+1)->pos);
	    p3 = &(get_ifk(i,j+1,k+1)->pos);
	    p4 = &(get_cell(i,j,k+1)->pos);
	    p5 = &(get_cell(i+1,j,k+1)->pos);
	    p6 = &(get_cell(i+1,j+1,k+1)->pos);
	    p7 = &(get_cell(i,j+1,k+1)->pos);
	    hex_cell_properties( *p0, *p1, *p2, *p3, *p4, *p5, *p6, *p7, 
				 dummy, vertex->volume, iLen, jLen, kLen );
	}
    }
    for ( i = imin; i <= imax; ++i) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    iface = get_sifi(i,j,k);
	    p0 = &(get_ifk(i,j,k+1)->pos);
	    p3 = &(get_ifk(i,j+1,k+1)->pos);
	    p7 = &(get_cell(i,j+1,k+1)->pos);
	    p4 = &(get_cell(i,j,k+1)->pos);
	    quad_properties( *p0, *p3, *p7, *p4,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    iface = get_sifj(i,j,k);
	    p0 = &(get_ifk(i,j,k+1)->pos);
	    p4 = &(get_cell(i,j,k+1)->pos);
	    p5 = &(get_cell(i+1,j,k+1)->pos);
	    p1 = &(get_ifk(i+1,j,k+1)->pos);
	    quad_properties( *p0, *p4, *p5, *p1,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }
    for ( i = imin; i <= imax-1; ++i ) {
	for ( j = jmin; j <= jmax-1; ++j ) {
	    iface = get_sifk(i,j,k);
	    p0 = &(get_ifk(i,j,k+1)->pos);
	    p1 = &(get_ifk(i+1,j,k+1)->pos);
	    p2 = &(get_ifk(i+1,j+1,k+1)->pos);
	    p3 = &(get_ifk(i,j+1,k+1)->pos);
	    quad_properties( *p0, *p1, *p2, *p3,
			     iface->pos, iface->n, iface->t1, iface->t2,
			     iface->area );
	}
    }

    return SUCCESS;
} // end compute_secondary_cell_geometric_data_3D()

int Block::calc_initial_volumes_2D( void )
/// \brief Compute the PRIMARY cell volumes, areas, and centers 
///        from the vertex positions.
///
/// For 2D planar, assume unit length in the Z-direction.
/// For axisymmetry, compute a volume per radian.
///
/// Determine minimum length and aspect ratio, also.
{
    int i, j;
    double xA, yA, xB, yB, xC, yC, xD, yD;
    double xN, yN, xS, yS, xE, yE, xW, yW;
    double vol, max_vol, min_vol, xyarea;
    double dx, dy, dxN, dyN, dxE, dyE;
    double lengthN, lengthE, length_max, length_min, length_cross;
    double max_aspect, aspect_ratio;
    FV_Cell *cell, *source_cell, *target_cell;

    // Cell layout
    // C-----B     3-----2
    // |     |     |     |
    // |  c  |     |  c  |
    // |     |     |     |
    // D-----A     0-----1
    
    max_vol = 0.0;
    min_vol = 1.0e6;    /* arbitrarily large */
    max_aspect = 0.0;
    for ( i = imin; i <= imax; ++i ) {
        for ( j = jmin; j <= jmax; ++j ) {
            cell = get_cell(i,j);
	    // These are the corners.
	    xA = cell->vtx[1]->pos.x;
	    yA = cell->vtx[1]->pos.y;
	    xB = cell->vtx[2]->pos.x;
	    yB = cell->vtx[2]->pos.y;
	    xC = cell->vtx[3]->pos.x;
	    yC = cell->vtx[3]->pos.y;
	    xD = cell->vtx[0]->pos.x;
	    yD = cell->vtx[0]->pos.y;
	    // Cell area in the (x,y)-plane.
            xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                            (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
	    // Cell Centroid.
            cell->pos.x = 1.0 / (xyarea * 6.0) * 
                      ((yB - yA) * (xA * xA + xA * xB + xB * xB) + 
                       (yC - yB) * (xB * xB + xB * xC + xC * xC) +
                       (yD - yC) * (xC * xC + xC * xD + xD * xD) + 
                       (yA - yD) * (xD * xD + xD * xA + xA * xA));
            cell->pos.y = -1.0 / (xyarea * 6.0) * 
                     ((xB - xA) * (yA * yA + yA * yB + yB * yB) + 
                      (xC - xB) * (yB * yB + yB * yC + yC * yC) +
                      (xD - xC) * (yC * yC + yC * yD + yD * yD) + 
                      (xA - xD) * (yD * yD + yD * yA + yA * yA));
	    cell->pos.z = 0.0;
	    // Cell Volume.
            if (get_axisymmetric_flag() == 1) {
                // Volume per radian = centroid y-ordinate * cell area
                vol = xyarea * cell->pos.y;
            } else {
                // Assume unit depth in the z-direction.
                vol = xyarea;
            }
            if (vol < 0.0) {
                printf("Negative cell volume: vol[%d][%d] = %e\n", i, j, vol);
                return 1;
            }
            if (vol > max_vol) max_vol = vol;
            if (vol < min_vol) min_vol = vol;
            cell->volume = vol;
            cell->area = xyarea;
	    // Check cell length scale using North and East boundaries.
	    // Also, save the minimum length for later use in the CFL
	    // checking routine.  Record max aspect ratio over all cells.
            dxN = xC - xB;
            dyN = yC - yB;
            dxE = xA - xB;
            dyE = yA - yB;
            lengthN = sqrt(dxN * dxN + dyN * dyN);
            lengthE = sqrt(dxE * dxE + dyE * dyE);

            length_max = lengthN;
            if (lengthE > length_max) length_max = lengthE;
            length_cross = xyarea / length_max; 
            // estimate of minimum width of cell
            length_min = lengthN;
            if (lengthE < length_min) length_min = lengthE;
            if (length_cross < length_min) length_min = length_cross;
            cell->L_min = length_min;
            aspect_ratio = length_max / length_min;
            if (aspect_ratio > max_aspect) max_aspect = aspect_ratio;

	    // Record the cell widths in the i- and j-index directions.
	    // The widths are measured between corresponding midpoints of
            // the bounding interfaces.
            // This data is later used by the high-order reconstruction.
            xN = 0.5 * (xC + xB);
            yN = 0.5 * (yC + yB);
            xS = 0.5 * (xD + xA);
            yS = 0.5 * (yD + yA);
            xE = 0.5 * (xA + xB);
            yE = 0.5 * (yA + yB);
            xW = 0.5 * (xD + xC);
            yW = 0.5 * (yD + yC);
            dx = xN - xS;
            dy = yN - yS;
            cell->jLength = sqrt(dx * dx + dy * dy);
            dx = xE - xW;
            dy = yE - yW;
            cell->iLength = sqrt(dx * dx + dy * dy);
	    cell->kLength = 0.0;
        } // j loop
    } // i loop

    // We now need to mirror the cell iLength and jLength
    // around the boundaries.
    // Those boundaries that are adjacent to another block
    // will be updated later with the other-block's cell lengths.
    for ( i = imin; i <= imax; ++i ) {
        // North boundary
        j = jmax;
        source_cell = get_cell(i,j);
	target_cell = get_cell(i,j+1);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
	target_cell = get_cell(i,j+2);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        // South boundary
        j = jmin;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i,j-1);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i,j-2);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
    } // end for i

    for ( j = jmin; j <= jmax; ++j ) {
        // East boundary
        i = imax;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i+1,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i+2,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        // West boundary
        i = imin;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i-1,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i-2,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
    } // end for j
	
    return SUCCESS;
} // end calc_initial_volumes_2D()

int Block::calc_volumes_2D( int time_level )
/// \brief Compute the PRIMARY cell volumes, areas, and centers 
///        from the vertex positions.
///
/// For 2D planar, assume unit length in the Z-direction.
/// For axisymmetry, compute a volume per radian.
///
/// Determine minimum length and aspect ratio, also.
{
    int i, j;
    double xA, yA, xB, yB, xC, yC, xD, yD;
    double xN, yN, xS, yS, xE, yE, xW, yW;
    double vol, max_vol, min_vol, xyarea;
    double dx, dy, dxN, dyN, dxE, dyE;
    double lengthN, lengthE, length_max, length_min, length_cross;
    double max_aspect, aspect_ratio;
    FV_Cell *cell, *source_cell, *target_cell;

    // Cell layout
    // C-----B     3-----2
    // |     |     |     |
    // |  c  |     |  c  |
    // |     |     |     |
    // D-----A     0-----1
    
    max_vol = 0.0;
    min_vol = 1.0e6;    /* arbitrarily large */
    max_aspect = 0.0;
    for ( i = imin; i <= imax; ++i ) {
        for ( j = jmin; j <= jmax; ++j ) {
            cell = get_cell(i,j);
	    // These are the corners.
	    xA = cell->vtx[1]->position[time_level].x;
	    yA = cell->vtx[1]->position[time_level].y;
	    xB = cell->vtx[2]->position[time_level].x;
	    yB = cell->vtx[2]->position[time_level].y;
	    xC = cell->vtx[3]->position[time_level].x;
	    yC = cell->vtx[3]->position[time_level].y;
	    xD = cell->vtx[0]->position[time_level].x;
	    yD = cell->vtx[0]->position[time_level].y;
	    // Cell area in the (x,y)-plane.
            xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                            (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
	    // Cell Centroid.
            cell->position[time_level].x = 1.0 / (xyarea * 6.0) * 
                      ((yB - yA) * (xA * xA + xA * xB + xB * xB) + 
                       (yC - yB) * (xB * xB + xB * xC + xC * xC) +
                       (yD - yC) * (xC * xC + xC * xD + xD * xD) + 
                       (yA - yD) * (xD * xD + xD * xA + xA * xA));
            cell->position[time_level].y = -1.0 / (xyarea * 6.0) * 
                     ((xB - xA) * (yA * yA + yA * yB + yB * yB) + 
                      (xC - xB) * (yB * yB + yB * yC + yC * yC) +
                      (xD - xC) * (yC * yC + yC * yD + yD * yD) + 
                      (xA - xD) * (yD * yD + yD * yA + yA * yA));
	    cell->position[time_level].z = 0.0;
	    // Cell Volume.
            if (get_axisymmetric_flag() == 1) {
                // Volume per radian = centroid y-ordinate * cell area
                vol = xyarea * cell->position[time_level].y;
            } else {
                // Assume unit depth in the z-direction.
                vol = xyarea;
            }
            if (vol < 0.0) {
                printf("Negative cell volume: vol[%d][%d] = %e\n", i, j, vol);
                return 1;
            }
            if (vol > max_vol) max_vol = vol;
            if (vol < min_vol) min_vol = vol;
            cell->vol[time_level] = vol;
            cell->ar[time_level] = xyarea;
	    // Check cell length scale using North and East boundaries.
	    // Also, save the minimum length for later use in the CFL
	    // checking routine.  Record max aspect ratio over all cells.
            dxN = xC - xB;
            dyN = yC - yB;
            dxE = xA - xB;
            dyE = yA - yB;
            lengthN = sqrt(dxN * dxN + dyN * dyN);
            lengthE = sqrt(dxE * dxE + dyE * dyE);

            length_max = lengthN;
            if (lengthE > length_max) length_max = lengthE;
            length_cross = xyarea / length_max; 
            // estimate of minimum width of cell
            length_min = lengthN;
            if (lengthE < length_min) length_min = lengthE;
            if (length_cross < length_min) length_min = length_cross;
            cell->L_min = length_min;
            aspect_ratio = length_max / length_min;
            if (aspect_ratio > max_aspect) max_aspect = aspect_ratio;

	    // Record the cell widths in the i- and j-index directions.
	    // The widths are measured between corresponding midpoints of
            // the bounding interfaces.
            // This data is later used by the high-order reconstruction.
            xN = 0.5 * (xC + xB);
            yN = 0.5 * (yC + yB);
            xS = 0.5 * (xD + xA);
            yS = 0.5 * (yD + yA);
            xE = 0.5 * (xA + xB);
            yE = 0.5 * (yA + yB);
            xW = 0.5 * (xD + xC);
            yW = 0.5 * (yD + yC);
            dx = xN - xS;
            dy = yN - yS;
            cell->jLength = sqrt(dx * dx + dy * dy);
            dx = xE - xW;
            dy = yE - yW;
            cell->iLength = sqrt(dx * dx + dy * dy);
	    cell->kLength = 0.0;
        } // j loop
    } // i loop

    // We now need to mirror the cell iLength and jLength
    // around the boundaries.
    // Those boundaries that are adjacent to another block
    // will be updated later with the other-block's cell lengths.
    for ( i = imin; i <= imax; ++i ) {
        // North boundary
        j = jmax;
        source_cell = get_cell(i,j);
	target_cell = get_cell(i,j+1);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
	target_cell = get_cell(i,j+2);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        // South boundary
        j = jmin;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i,j-1);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i,j-2);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
    } // end for i

    for ( j = jmin; j <= jmax; ++j ) {
        // East boundary
        i = imax;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i+1,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i+2,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        // West boundary
        i = imin;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i-1,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
        source_cell = get_cell(i,j);
        target_cell = get_cell(i-2,j);
        target_cell->iLength = source_cell->iLength;
        target_cell->jLength = source_cell->jLength;
	target_cell->kLength = 0.0;
    } // end for j
	
    return SUCCESS;
} // end calc_volumes_2D()


/// \brief Compute the secondary cell cell areas in the (x,y)-plane.
///
/// The secondary cells are centred on the vertices of the 
/// primary cells and have primary cell centres as their corners.
/// For this particular secondary cell, centred on a vertex v (i,j),
/// the near-by primary cell i,j is centred on B'.
///
///          +-----+
///          |     |
///       C'-+--B' |
///       |  |  |  |
///       |  v--+--+
///       |     |
///       D'----A'
///
int Block::secondary_areas_2D( void )
{
    int i, j;
    double xA, yA, xB, yB, xC, yC, xD, yD;
    double xyarea, max_area, min_area;

    max_area = 0.0;
    min_area = 1.0e6;   // arbitrarily large
    // First, do all of the internal secondary cells.
    // i.e. The ones centred on primary vertices which 
    // are not on a boundary.
    for (i = imin+1; i <= imax; ++i) {
        for (j = jmin+1; j <= jmax; ++j) {
	    // These are the corners.
            xA = get_cell(i,j-1)->pos.x;
            yA = get_cell(i,j-1)->pos.y;
            xB = get_cell(i,j)->pos.x;
            yB = get_cell(i,j)->pos.y;
            xC = get_cell(i-1,j)->pos.x;
            yC = get_cell(i-1,j)->pos.y;
            xD = get_cell(i-1,j-1)->pos.x;
            yD = get_cell(i-1,j-1)->pos.y;
	    // Cell area in the (x,y)-plane.
            xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                            (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
            if (xyarea < 0.0) {
                printf("Negative secondary-cell area: Block %d, vtx[%d,%d] = %e\n", 
		       id, i, j, xyarea);
                return BAD_CELLS_ERROR;
            }
            if (xyarea > max_area) max_area = xyarea;
            if (xyarea < min_area) min_area = xyarea;
            get_vtx(i,j)->area = xyarea;
        } // j loop
    } // i loop

    // Note that the secondary cells along block boundaries are HALF cells.
    //
    // East boundary.
    i = imax+1;
    for (j = jmin+1; j <= jmax; ++j) {
	xA = get_ifi(i,j-1)->pos.x;
	yA = get_ifi(i,j-1)->pos.y;
	xB = get_ifi(i,j)->pos.x;
	yB = get_ifi(i,j)->pos.y;
        xC = get_cell(i-1,j)->pos.x;
        yC = get_cell(i-1,j)->pos.y;
        xD = get_cell(i-1,j-1)->pos.x;
        yD = get_cell(i-1,j-1)->pos.y;
	// Cell area in the (x,y)-plane.
        xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                        (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
        if (xyarea < 0.0) {
            printf("Negative secondary-cell area: Block %d, vtx[%d,%d] = %e\n", 
		   id, i, j, xyarea);
            return BAD_CELLS_ERROR;
        }
        if (xyarea > max_area) max_area = xyarea;
        if (xyarea < min_area) min_area = xyarea;
        get_vtx(i,j)->area = xyarea;
    } // j loop 

    // Fudge corners -- not expecting to use this data.
    get_vtx(i,jmin)->area = 0.5 * get_vtx(i,jmin+1)->area;
    get_vtx(i,jmax+1)->area = 0.5 * get_vtx(i,jmax)->area;
    
    // West boundary.
    i = imin;
    for (j = jmin+1; j <= jmax; ++j) {
        xA = get_cell(i,j-1)->pos.x;
        yA = get_cell(i,j-1)->pos.y;
        xB = get_cell(i,j)->pos.x;
        yB = get_cell(i,j)->pos.y;
	xC = get_ifi(i,j)->pos.x;
	yC = get_ifi(i,j)->pos.y;
	xD = get_ifi(i,j-1)->pos.x;
	yD = get_ifi(i,j-1)->pos.y;
	// Cell area in the (x,y)-plane.
        xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                        (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
        if (xyarea < 0.0) {
	    printf("Negative secondary-cell area: Block %d, vtx[%d,%d] = %e\n", 
		   id, i, j, xyarea);
            return BAD_CELLS_ERROR;
        }
        if (xyarea > max_area) max_area = xyarea;
        if (xyarea < min_area) min_area = xyarea;
        get_vtx(i,j)->area = xyarea;
    } // j loop 

    // Fudge corners.
    get_vtx(i,jmin)->area = 0.5 * get_vtx(i,jmin+1)->area;
    get_vtx(i,jmax+1)->area = 0.5 * get_vtx(i,jmax)->area;

    // North boundary.
    j = jmax+1;
    for (i = imin+1; i <= imax; ++i) {
	// These are the corners.
        xA = get_cell(i,j-1)->pos.x;
        yA = get_cell(i,j-1)->pos.y;
	xB = get_ifj(i,j)->pos.x;
	yB = get_ifj(i,j)->pos.y;
	xC = get_ifj(i-1,j)->pos.x;
	yC = get_ifj(i-1,j)->pos.y;
        xD = get_cell(i-1,j-1)->pos.x;
        yD = get_cell(i-1,j-1)->pos.y;
	// Cell area in the (x,y)-plane.
        xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                        (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
        if (xyarea < 0.0) {
            printf("Negative secondary-cell area: Block %d, vtx[%d,%d] = %e\n", 
		   id, i, j, xyarea);
            return BAD_CELLS_ERROR;
        }
        if (xyarea > max_area) max_area = xyarea;
        if (xyarea < min_area) min_area = xyarea;
        get_vtx(i,j)->area = xyarea;
    } // i loop 

    // Fudge corners.
    get_vtx(imin,j)->area = 0.5 * get_vtx(imin+1,j)->area;
    get_vtx(imax+1,j)->area = 0.5 * get_vtx(imax,j)->area;

    // South boundary.
    j = jmin;
    for (i = imin+1; i <= imax; ++i) {
	xA = get_ifj(i,j)->pos.x;
	yA = get_ifj(i,j)->pos.y;
        xB = get_cell(i,j)->pos.x;
        yB = get_cell(i,j)->pos.y;
        xC = get_cell(i-1,j)->pos.x;
        yC = get_cell(i-1,j)->pos.y;
	xD = get_ifj(i-1,j)->pos.x;
	yD = get_ifj(i-1,j)->pos.y;
	// Cell area in the (x,y)-plane.
        xyarea = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                        (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
        if (xyarea < 0.0) {
            printf("Negative secondary-cell area: Block %d, vtx[%d,%d] = %e\n", 
		   id, i, j, xyarea);
            return BAD_CELLS_ERROR;
        }
        if (xyarea > max_area) max_area = xyarea;
        if (xyarea < min_area) min_area = xyarea;
        get_vtx(i,j)->area = xyarea;
    } // i loop

    // Fudge corners.
    get_vtx(imin,j)->area = 0.5 * get_vtx(imin+1,j)->area;
    get_vtx(imax+1,j)->area = 0.5 * get_vtx(imax,j)->area;

    //printf("Max area = %e, Min area = %e\n", max_area, min_area);

    return SUCCESS;
} // end secondary_areas_2D()

/// \brief Compute the interface lengths and direction cosines for interfaces.
///
int Block::calc_initial_faces_2D( void )
{
    FV_Interface *IFace;
    int i, j;
    double xA, xB, yA, yB, xC, yC;
    double LAB, LBC;

    // East-facing interfaces.
    for (i = imin; i <= imax+1; ++i) {
        for (j = jmin; j <= jmax; ++j) {
            IFace = get_ifi(i,j);
	    // These are the corners.
            xA = get_vtx(i,j)->pos.x; 
	    yA = get_vtx(i,j)->pos.y;
            xB = get_vtx(i,j+1)->pos.x; 
	    yB = get_vtx(i,j+1)->pos.y;
            LAB = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA));
            if (LAB < 1.0e-9) {
                printf("Zero length ifi[%d,%d]: %e\n", i, j, LAB);
            }
            // Direction cosines for the unit normal.
            IFace->n.x = (yB - yA) / LAB;
            IFace->n.y = -(xB - xA) / LAB;
	    IFace->n.z = 0.0;  // 2D plane
	    IFace->t2 = Vector3(0.0, 0.0, 1.0);
	    IFace->t1 = cross(IFace->n, IFace->t2);
            // Length in the XY-plane.
            IFace->length = LAB;
            // Mid-point and area.
            IFace->Ybar = 0.5 * (yA + yB);
            if (get_axisymmetric_flag() == 1) {
                // Interface area per radian.
                IFace->area = LAB * IFace->Ybar;
            } else {
                // Assume unit depth in the Z-direction.
                IFace->area = LAB;
            }
	    IFace->pos = (get_vtx(i,j)->pos + get_vtx(i,j+1)->pos)/2.0;
	    
        } // j loop
    } // i loop
    
    // North-facing interfaces.
    for (i = imin; i <= imax; ++i) {
        for (j = jmin; j <= jmax+1; ++j) {
            IFace = get_ifj(i,j);
	    // These are the corners.
            xB = get_vtx(i+1,j)->pos.x;
            yB = get_vtx(i+1,j)->pos.y;
            xC = get_vtx(i,j)->pos.x;
            yC = get_vtx(i,j)->pos.y;
            LBC = sqrt((xC - xB) * (xC - xB) + (yC - yB) * (yC - yB));
            if (LBC < 1.0e-9) {
                printf("Zero length ifj[%d,%d]: %e\n", i, j, LBC);
            }
            // Direction cosines for the unit normal.
            IFace->n.x = (yC - yB) / LBC;
            IFace->n.y = -(xC - xB) / LBC;
	    IFace->n.z = 0.0;  // 2D plane
	    IFace->t2 = Vector3(0.0, 0.0, 1.0);
	    IFace->t1 = cross(IFace->n, IFace->t2);
            // Length in the XY-plane.
            IFace->length = LBC;
            // Mid-point and area.
            IFace->Ybar = 0.5 * (yC + yB);
            if (get_axisymmetric_flag() == 1) {
                // Interface area per radian.
                IFace->area = LBC * IFace->Ybar;
            } else {
                // Assume unit depth in the Z-direction.
                IFace->area = LBC;
            }
	    IFace->pos = (get_vtx(i+1,j)->pos + get_vtx(i,j)->pos)/2.0;
        } // j loop
    } // i loop
    return SUCCESS;
} // end calc_initial_faces_2D()

/// \brief Compute the interface lengths and direction cosines for interfaces.
///
int Block::calc_faces_2D( int time_level )
{
    FV_Interface *IFace;
    int i, j;
    double xA, xB, yA, yB, xC, yC;
    double LAB, LBC;

    // East-facing interfaces.
    for (i = imin; i <= imax+1; ++i) {
        for (j = jmin; j <= jmax; ++j) {
            IFace = get_ifi(i,j);
	    // These are the corners.
            xA = get_vtx(i,j)->position[time_level].x; 
	    yA = get_vtx(i,j)->position[time_level].y;
            xB = get_vtx(i,j+1)->position[time_level].x; 
	    yB = get_vtx(i,j+1)->position[time_level].y;
            LAB = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA));
            if (LAB < 1.0e-9) {
                printf("Zero length ifi[%d,%d]: %e\n", i, j, LAB);
            }
            // Direction cosines for the unit normal.
            IFace->n.x = (yB - yA) / LAB;
            IFace->n.y = -(xB - xA) / LAB;
	    IFace->n.z = 0.0;  // 2D plane
	    IFace->t2 = Vector3(0.0, 0.0, 1.0);
	    IFace->t1 = cross(IFace->n, IFace->t2);
            // Length in the XY-plane.
            IFace->length = LAB;
            // Mid-point and area.
            IFace->Ybar = 0.5 * (yA + yB);
            if (get_axisymmetric_flag() == 1) {
                // Interface area per radian.
                IFace->ar[time_level] = LAB * IFace->Ybar;
            } else {
                // Assume unit depth in the Z-direction.
                IFace->ar[time_level] = LAB;
            }
	    IFace->pos = (get_vtx(i,j)->position[time_level] + get_vtx(i,j+1)->position[time_level])/2.0;
	    
        } // j loop
    } // i loop
    
    // North-facing interfaces.
    for (i = imin; i <= imax; ++i) {
        for (j = jmin; j <= jmax+1; ++j) {
            IFace = get_ifj(i,j);
	    // These are the corners.
            xB = get_vtx(i+1,j)->position[time_level].x;
            yB = get_vtx(i+1,j)->position[time_level].y;
            xC = get_vtx(i,j)->position[time_level].x;
            yC = get_vtx(i,j)->position[time_level].y;
            LBC = sqrt((xC - xB) * (xC - xB) + (yC - yB) * (yC - yB));
            if (LBC < 1.0e-9) {
                printf("Zero length ifj[%d,%d]: %e\n", i, j, LBC);
            }
            // Direction cosines for the unit normal.
            IFace->n.x = (yC - yB) / LBC;
            IFace->n.y = -(xC - xB) / LBC;
	    IFace->n.z = 0.0;  // 2D plane
	    IFace->t2 = Vector3(0.0, 0.0, 1.0);
	    IFace->t1 = cross(IFace->n, IFace->t2);
            // Length in the XY-plane.
            IFace->length = LBC;
            // Mid-point and area.
            IFace->Ybar = 0.5 * (yC + yB);
            if (get_axisymmetric_flag() == 1) {
                // Interface area per radian.
                IFace->ar[time_level] = LBC * IFace->Ybar;
            } else {
                // Assume unit depth in the Z-direction.
                IFace->ar[time_level] = LBC;
            }
	    IFace->pos = (get_vtx(i+1,j)->position[time_level] + get_vtx(i,j)->position[time_level])/2.0;
        } // j loop
    } // i loop
    return SUCCESS;
} // end calc_faces_2D()

int Block::calc_initial_ghost_cell_geom_2D( void )
/// \brief Compute the ghost cell positions and volumes.
///
/// 'Compute' is a bit too strong to describe what we do here.
///  Rather this is a first-order extrapolation
/// from interior cells to estimate the position
/// and volume of the ghost cells.
{
    int i, j;
    FV_Cell *cell_1, *cell_2, *ghost_cell;
    // East boundary
    i = imax;
    for ( j = jmin; j <= jmax; ++j ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i-1,j);
	ghost_cell = get_cell(i+1,j);
	ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i+2,j);
	ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
    }
    // West boundary
    i = imin;
    for ( j = jmin; j <= jmax; ++j ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i+1,j);
	ghost_cell = get_cell(i-1,j);
	ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i-2,j);
	ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
    }
    // North boundary
    j = jmax;
    for ( i = imin; i <= imax; ++i ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i,j-1);
	ghost_cell = get_cell(i,j+1);
	ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i,j+2);
	ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
    }
    // South boundary
    j = jmin;
    for ( i = imin; i <= imax; ++i ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i,j+1);
	ghost_cell = get_cell(i,j-1);
	ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i,j-2);
	ghost_cell->pos = 2.0*cell_1->pos - cell_2->pos;
	ghost_cell->volume = 2.0*cell_1->volume - cell_2->volume;
    }
    return SUCCESS;
}

int Block::calc_ghost_cell_geom_2D( int time_level )
/// \brief Compute the ghost cell positions and volumes.
///
/// 'Compute' is a bit too strong to describe what we do here.
///  Rather this is a first-order extrapolation
/// from interior cells to estimate the position
/// and volume of the ghost cells.
{
    int i, j;
    FV_Cell *cell_1, *cell_2, *ghost_cell;
    // East boundary
    i = imax;
    for ( j = jmin; j <= jmax; ++j ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i-1,j);
	ghost_cell = get_cell(i+1,j);
	ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i+2,j);
	ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
    }
    // West boundary
    i = imin;
    for ( j = jmin; j <= jmax; ++j ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i+1,j);
	ghost_cell = get_cell(i-1,j);
	ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i-2,j);
	ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
    }
    // North boundary
    j = jmax;
    for ( i = imin; i <= imax; ++i ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i,j-1);
	ghost_cell = get_cell(i,j+1);
	ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i,j+2);
	ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
    }
    // South boundary
    j = jmin;
    for ( i = imin; i <= imax; ++i ) {
	cell_1 = get_cell(i,j);
	cell_2 = get_cell(i,j+1);
	ghost_cell = get_cell(i,j-1);
	ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
	cell_2 = cell_1;
	cell_1 = ghost_cell;
	ghost_cell = get_cell(i,j-2);
	ghost_cell->position[time_level] = 2.0*cell_1->position[time_level] - cell_2->position[time_level];
	ghost_cell->vol[time_level] = 2.0*cell_1->vol[time_level] - cell_2->vol[time_level];
    }
    return SUCCESS;
}

/// \brief Perform Euler step using calculated vertex
/// velocity and set the vertex position for time-level 1.
///
int Block::predict_vertex_positions( int dimensions, double dt )
{
    int i, j, k, krangemax;
    double gamma_1;
    FV_Vertex *vtx;
    if (get_Torder_flag() == 3) {
	/* 3rd order Runge-Kutta */
	gamma_1 = 8.0 / 15.0;
    } else {
	/* Normal Predictor-Corrector or Euler */
	gamma_1 = 1.0;
    }
    if ( dimensions == 2 ) {
	krangemax = kmax;
    } else {
	krangemax = kmax+1;
    }
    for (k = kmin; k <= krangemax; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
		vtx = get_vtx(i,j,k);
		vtx->position[1] = vtx->position[0] + dt * gamma_1 * vtx->velocity[0];
	    }
	}
    }
    return SUCCESS;
}

/// \brief Perform corrector step using re-calculated vertex
/// velocity and set the vertex position for time-level 2.
///
int Block::correct_vertex_positions( int dimensions, double dt )
{
    int i, j, k, krangemax;
    FV_Vertex *vtx;
    double th, th_inv;
    int tl_old;
    if (get_Torder_flag() == 3) {
	/* 3rd order Runge-Kutta */
	th = 5.0 / 12.0;
	th_inv = -17.0 / 60.0;
	tl_old = 1;
    } else {
	/* Normal Predictor-Corrector or Euler */
	th = 0.5;
	th_inv = 0.5;
	tl_old = 0;
    }
    if ( dimensions == 2 ) {
	krangemax = kmax;
    } else {
	krangemax = kmax+1;
    }
    for (k = kmin; k <= krangemax; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
		vtx = get_vtx(i,j,k);
		vtx->position[2] = vtx->position[tl_old] + dt * (th_inv * vtx->velocity[0] + 
								 th * vtx->velocity[1]);
	    }
	}
    }
    return SUCCESS;
}

/// \brief Perform corrector step using re-calculated vertex
/// velocity and set the vertex position for time-level 2.
///
int Block::rk3_vertex_positions( int dimensions, double dt )
{
    int i, j, k, krangemax;
    FV_Vertex *vtx;
    double gamma_3 = 3.0 / 4.0;
    double psi_2 = -5.0 / 12.0;
    if ( dimensions == 2 ) {
	krangemax = kmax;
    } else {
	krangemax = kmax+1;
    }
    for (k = kmin; k <= krangemax; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
		vtx = get_vtx(i,j,k);
		vtx->position[3] = vtx->position[2] + dt * (psi_2 * vtx->velocity[1] + 
							    gamma_3 * vtx->velocity[2]);
	    }
	}
    }
    return SUCCESS;
}


/// \brief Calculate shock speed at interface for 2D. 
/// See Ian Johnston's thesis for an explanation.
///
int Block::calc_boundary_vertex_velocity(FV_Interface &IFace1, FV_Interface &IFace2,     
                                         FV_Vertex &vtx, Vector3 trv, int time_level )
{   
    double w1, w2;
    Vector3 ws1, ws2, vp;
    vp = vtx.position[time_level];
    velocity_weighting_factor(IFace1, vp, w1, ws1);
    velocity_weighting_factor(IFace2, vp, w2, ws2);
    if ( (w1 + w2) < 1e-3 ) {
	w1 = w2 = 1.0;
    }
    Vector3 wv = (w1*ws1 + w2*ws2) / (w1 + w2);
    vtx.vel = dot(wv, trv) * trv; // Constrain vertex velocity to body radial direction.
    return SUCCESS;				       			           
}

/// \brief Calculate shock speed at interface for 3D. 
/// See Ian Johnston's thesis for an explanation.
/// 
int Block::calc_boundary_vertex_velocity(FV_Interface &IFace1, FV_Interface &IFace2,     
                                         FV_Interface &IFace3, FV_Interface &IFace4,
                                         FV_Vertex &vtx, Vector3 trv, int time_level)
{   
    double w1, w2, w3, w4;
    Vector3 ws1, ws2, ws3, ws4, vp;
    vp = vtx.position[time_level];
    velocity_weighting_factor(IFace1, vp, w1, ws1);
    velocity_weighting_factor(IFace2, vp, w2, ws2);
    velocity_weighting_factor(IFace3, vp, w3, ws3);
    velocity_weighting_factor(IFace4, vp, w4, ws4);
    if ( (w1 + w2 + w3 + w4) < 1e-3 ) {
	w1 = w2 = w3 = w4 = 1.0;
    }
    Vector3 wv = (w1*ws1 + w2*ws2 + w3*ws3 + w4*ws4) / (w1 + w2 + w3 + w4);
    vtx.vel = dot(wv, trv) * trv; // Constrain vertex velocity to body radial direction.
    return SUCCESS;					       			           
}

int Block::velocity_weighting_factor(FV_Interface &IFace, Vector3 vp, double &w, Vector3 &ws)
{
    double M = dot(IFace.fs->vel, unit(vp - IFace.pos)) / IFace.fs->gas->a;
    if ( vabs(IFace.vel) == vabs(IFace.vel) ) // If not NaN.
	ws = IFace.vel;
    else {
	ws = 0.0;
    }
    if ( M > 1 ) {
        w = M;
    } else if ( M != M ) { // Check for NaN.
	w = 0.0;
    } else {
        w = 0.125*( pow(M+1, 2) + (M+1)*fabs(M+1) );
    }
    return SUCCESS;
}

int Block::set_geometry_velocities(int dimensions, int time_level)
{
    if ( dimensions == 2 ) { 
	set_vertex_velocities2D(time_level);
	set_interface_velocities2D(time_level);
    } else {
	set_vertex_velocities3D(time_level);
	set_interface_velocities3D(time_level);
    }
    return SUCCESS;
}

/// \brief Set vertex velocities based on previously calculated boundary interface velocities in 2D. 
///
///  Based on Ian Johnston's thesis, see for explanation.
///  Assumes inflow at west boundary and wall at east boundary.
///
int Block::set_vertex_velocities2D( int time_level )
{
    // Only works with one block in the i-direction. Supports multiple blocks
    // in the j-direction.
    int i, j, k;
    FV_Interface *IFaceU, *IFaceD;
    FV_Vertex *svtx, *wvtx, *vtx;
    Vector3 trv;
    double length;

    i = imin;
    k = kmin;
    // Set boundary vertex velocities.
    // Ghost cell geometry will be invalid, but NaNs will be caught by the weighting function.
    for (j = jmin; j <= jmax+1; ++j) {
	IFaceD = get_ifi(i,j-1,k);
	IFaceU = get_ifi(i,j,k);
	vtx = get_vtx(i,j,k);
	wvtx = get_vtx(imax,j,k);
	// Direction vector from vertex to body.
	trv = unit(wvtx->position[time_level] - vtx->position[time_level]); 
	calc_boundary_vertex_velocity(*IFaceD, *IFaceU, *vtx, trv, time_level);
    } // for j
    // // Set first and last two boundary vertex velocities
    // if ( bcp[SOUTH]->type_code != ADJACENT ) { // If not adjacent to another block on south side.
    // 	// First
    // 	vtx = get_vtx(imin,jmin,k);
    // 	wvtx = get_vtx(imax,jmin,k);
    // 	IFaceU = get_ifi(imin,jmin,k);
    // 	trv = unit(wvtx->position[time_level] - vtx->position[time_level]); 	    
    // 	vtx->vel = dot(IFaceU->vel, trv) * trv; // Constrain vertex velocity to body radial direction.
    // }
    // if ( bcp[NORTH]->type_code != ADJACENT ) { // If not adjacent to another block on north side.
    // 	// Last
    // 	vtx = get_vtx(imin,jmax+1,k);
    // 	wvtx = get_vtx(imax,jmax+1,k);
    // 	IFaceD = get_ifi(imin,jmax,k);
    // 	trv = unit(wvtx->position[time_level] - vtx->position[time_level]);
    // 	vtx->vel = dot(IFaceD->vel, trv) * trv; // Constrain vertex velocity to body radial direction.
    // }
    // Set interior vertex velocities.
    // Velocities are set as linear functions of position between
    // the shock boundary and the wall.
    for (j = jmin; j <= jmax+1; ++j) {
	svtx = get_vtx(imin,j,k); // Shock boundary vertex.
	wvtx = get_vtx(imax+1,j,k); // Wall vertex.
	length = vabs(svtx->position[time_level] - wvtx->position[time_level]);
	for (i = imin; i <= imax; ++i) {
	    vtx = get_vtx(i,j,k);
	    vtx->velocity[time_level] = (vabs(vtx->position[time_level] - 
					      wvtx->position[time_level]) 
					 / length ) * svtx->vel;
	}
    }
    return SUCCESS;
}

/// \brief Set vertex velocities based on previously calculated boundary interface velocities in 3D.
///
///  Based on Ian Johnston's thesis, see for explanation.
///  Assumes inflow at west boundary and wall at east boundary.
///
int Block::set_vertex_velocities3D( int time_level )
{
    // Only works with one block in the i-direction. Supports multiple blocks
    // in the j-direction.
    int i, j, k;
    FV_Interface *IFace1, *IFace2, *IFace3, *IFace4;
    FV_Vertex *svtx, *wvtx, *vtx;
    Vector3 trv;
    double length;
    i = imin;
    // Set boundary vertex velocities.
    // Ghost cell geometry will be invalid, but NaNs will be caught by the weighting function.
    for (k = kmin; k <= kmax+1; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    IFace1 = get_ifi(i,j,k);
	    IFace2 = get_ifi(i,j-1,k);
	    IFace3 = get_ifi(i,j,k-1);
	    IFace4 = get_ifi(i,j-1,k-1);
	    vtx = get_vtx(i,j,k);
	    wvtx = get_vtx(imax,j,k);
	    // Direction vector from vertex to body.
	    trv = unit(wvtx->position[time_level] - vtx->position[time_level]); 
	    calc_boundary_vertex_velocity(*IFace1, *IFace2, *IFace3, *IFace4, *vtx, trv, time_level);
	} // for j
    } // for k
   // // Set vertex velocities on edges of domain as ghost cell will not have valid data.
   //  if ( bcp[SOUTH]->type_code != ADJACENT ) { // If not adjacent to another block on south side.
   // 	for (k = kmin; k <= kmax+1; ++k) {
   // 	    vtx = get_vtx(imin,jmin,kmin);
   // 	    wvtx = get_vtx(imax,jmin,kmin);
   // 	    IFace1 = get_ifi(imin,jmin,kmin);
   // 	    IFace1 = get_ifi(imin,jmin,kmin);
   // 	    trv = unit(wvtx->position[time_level] - vtx->position[time_level]); 	    
   // 	    calc_boundary_vertex_velocity(*IFace1, *IFace2, *vtx, trv, time_level);
   // 	}
   //  }
   //  if ( bcp[NORTH]->type_code != ADJACENT ) { // If not adjacent to another block on north side.
   // 	for (k = kmin; k <= kmax+1; ++k) {
   // 	    vtx = get_vtx(imin,jmax+1,kmin);
   // 	    wvtx = get_vtx(imax,jmax+1,kmin);
   // 	    IFace1 = get_ifi(imin,jmax,kmin);
   // 	    trv = unit(wvtx->position[time_level] - vtx->position[time_level]);
   // 	    calc_boundary_vertex_velocity(*IFace1, *IFace2, *vtx, trv, time_level);
   // 	}
   //  }
    // Set interior vertex velocities.
    // Velocities are set as linear functions of position between
    // the shock boundary and the wall.
    for (k = kmin; k <= kmax+1; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    svtx = get_vtx(imin,j,k); // Shock boundary vertex.
	    wvtx = get_vtx(imax+1,j,k); // Wall vertex.
	    length = vabs(svtx->position[time_level] - wvtx->position[time_level]);
	    for (i = imin; i <= imax; ++i) {
		vtx = get_vtx(i,j,k);
		vtx->velocity[time_level] = (vabs(vtx->position[time_level] - 
						  wvtx->position[time_level]) 
					     / length ) * svtx->vel;
	    }
	}
    }
    return SUCCESS;
}


/// \brief Function used to test GCL adherence.
///
int Block::set_gcl_test_vertex_velocities2D( int time_level )
{
    // Only works with one block in the i-direction. Supports multiple blocks
    // in the j-direction.
    int i, j, k;
    FV_Vertex *vtx;
    k = kmin;
    // Set boundary vertex velocities.
    // Ghost cell geometry will be invalid, but NaNs will be caught by the weighting function.
    for (j = jmin; j <= jmax+1; ++j) {
	for (i = imin; i <= imax+1; ++i) {
	    vtx = get_vtx(i,j,k);
	    vtx->velocity[time_level] = 1.0 * vtx->position[time_level];
	}
    } // for j
    return SUCCESS;
}

/// \brief Function used to test GCL adherence. 
///
int Block::set_gcl_test_vertex_velocities3D( int time_level )
{
    // Only works with one block in the i-direction. Supports multiple blocks
    // in the j-direction.
    int i, j, k;
    FV_Vertex *vtx;
    // Set boundary vertex velocities.
    // Ghost cell geometry will be invalid, but NaNs will be caught by the weighting function.
    for (k = kmin; k <= kmax+1; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
		vtx = get_vtx(i,j,k);
		vtx->velocity[time_level] = 1.0 * vtx->position[time_level];
	    }
	}
    }
    return SUCCESS;
}

/// \brief Set vertex velocities based on previously calculated boundary interface velocities in 2D. 
///
///  Based on Ian Johnston's thesis, see for explanation.
///  Assumes inflow at west boundary and wall at east boundary.
///
int Block::set_gcl_test_random_vertex_velocities2D( int time_level )
{
    // Only works with one block in the i-direction. Supports multiple blocks
    // in the j-direction.
    int i, j, k;
    FV_Vertex *vtx;
    FV_Cell *cell;
    k = kmin;
    // Set boundary vertex velocities.
    // Ghost cell geometry will be invalid, but NaNs will be caught by the weighting function.
    srand ( 1 );
    for (j = jmin+1; j <= jmax; ++j) {
	for (i = imin+1; i <= imax; ++i) {
	    vtx = get_vtx(i,j,k);
	    cell = get_cell(i,j,k);
	    vtx->velocity[time_level].x = (double)(rand() % (int)( cell->fs->gas->a / 50.0 ));
	    vtx->velocity[time_level].y = (double)(rand() % (int)( cell->fs->gas->a / 50.0 ));
	    vtx->velocity[time_level].z = 0.0;
	    cout << vtx->velocity[time_level] << endl;
	    cout << vabs(vtx->velocity[time_level]) << endl;
	    cout << cell->fs->gas->a / 50.0 << endl;
	}
    } // for j
    return SUCCESS;
}


int Block::set_gcl_interface_properties(int dimensions, int time_level, double dt)
{
    if ( dimensions == 2 ) { 
	set_gcl_interface_properties2D(time_level, dt);

    } else {
	set_gcl_interface_properties3D(time_level, dt);
    }
    return SUCCESS;
}

/// \brief Set interface velocities and area to average value over timestep for GCL adherence. 
///
int Block::set_gcl_interface_properties2D( int time_level, double dt )
{
    int i, j, k;
    FV_Vertex *vtx1, *vtx2;
    FV_Interface *IFace;
    Vector3 vpm1, vpm2;
    double xA, xB, yA, yB;
    int tl_old = 0;
    if ( get_Torder_flag() == 3 ) {
	if ( time_level == 1 ) tl_old = 1;
	else if ( time_level == 2 ) tl_old = 2;
	dt /= 3.0;
    }
    k = kmin;
    for (j = jmin; j <= jmax; ++j) {
	for (i = imin; i <= imax+1; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i,j+1,k);
	    IFace = get_ifi(i,j,k);   
	    vpm1 = 0.5 * ( vtx1->position[tl_old] + vtx1->position[time_level+1] );
	    vpm2 = 0.5 * ( vtx2->position[tl_old] + vtx2->position[time_level+1] );
	    IFace->pos = 0.5 * (vpm1 + vpm2);
	    IFace->vel = 0.5 * (vtx1->position[time_level+1] + vtx2->position[time_level+1] - 
	    			vtx1->position[tl_old] - vtx2->position[tl_old]) / dt;
            xA = vpm1.x;
	    yA = vpm1.y;
            xB = vpm2.x;
	    yB = vpm2.y;	 
	    // Interface area at midpoint.   
	    IFace->ar[time_level] = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA)); 
	    if (get_axisymmetric_flag() == 1) {
		IFace->Ybar = 0.5 * (yA + yB);
                IFace->ar[time_level] *= IFace->Ybar;
            }
	}
    }
    for (j = jmin; j <= jmax+1; ++j) {
	for (i = imin; i <= imax; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i+1,j,k);
	    IFace = get_ifj(i,j,k);
	    vpm1 = 0.5 * ( vtx1->position[tl_old] + vtx1->position[time_level+1] );
	    vpm2 = 0.5 * ( vtx2->position[tl_old] + vtx2->position[time_level+1] );
	    IFace->pos = 0.5 * (vpm1 + vpm2);
	    IFace->vel = 0.5 * (vtx1->position[time_level+1] + vtx2->position[time_level+1] - 
	    			vtx1->position[tl_old] - vtx2->position[tl_old]) / dt;
            xA = vpm2.x;
	    yA = vpm2.y;
            xB = vpm1.x;
	    yB = vpm1.y;
	    // Interface area at midpoint.   
	    IFace->ar[time_level] = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA)); 
	    if (get_axisymmetric_flag() == 1) {
		IFace->Ybar = 0.5 * (yA + yB);
                IFace->ar[time_level] *= IFace->Ybar;
            }
	}
    }
    return SUCCESS;
}

/// \brief Set interface velocities and area to average value over timestep for GCL adherence. 
///
int Block::set_gcl_interface_properties3D( int time_level, double dt )
{
    int i, j, k;
    FV_Vertex *vtx1, *vtx2, *vtx3, *vtx4;
    FV_Interface *IFace;
    Vector3 vpm1, vpm2, vpm3, vpm4, p1, p2, p3, p4;
    int tl_old = 0;
    if ( get_Torder_flag() == 3 ) {
	if ( time_level == 1 ) tl_old = 1;
	else if ( time_level == 2 ) tl_old = 2;
	dt /= 3.0;
    }
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i,j+1,k);
		vtx3 = get_vtx(i,j,k+1);
		vtx4 = get_vtx(i,j+1,k+1);
		IFace = get_ifi(i,j,k);   
		vpm1 = 0.5 * ( vtx1->position[tl_old] + vtx1->position[time_level+1] );
		vpm2 = 0.5 * ( vtx2->position[tl_old] + vtx2->position[time_level+1] );
		vpm3 = 0.5 * ( vtx3->position[tl_old] + vtx3->position[time_level+1] );
		vpm4 = 0.5 * ( vtx4->position[tl_old] + vtx4->position[time_level+1] );
		IFace->pos = 0.25 * (vpm1 + vpm2 + vpm3 + vpm4);
		IFace->vel = 0.25 * (vtx1->position[time_level+1] + vtx2->position[time_level+1] +
				     vtx3->position[time_level+1] + vtx4->position[time_level+1] - 
				     vtx1->position[tl_old] - vtx2->position[tl_old] - 
				     vtx3->position[tl_old] - vtx4->position[tl_old]) / dt;
		p1 = vpm1;
		p4 = vpm2;
		p2 = vpm3;
		p3 = vpm4;
		// Interface area at midpoint.
		IFace->ar[time_level] = vabs(0.25 * cross(p2-p1+p3-p4, p4-p1+p3-p2)); 
	    }
	}
    }
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    for (i = imin; i <= imax; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i+1,j,k);
		vtx3 = get_vtx(i,j,k+1);
		vtx4 = get_vtx(i+1,j,k+1);
		IFace = get_ifj(i,j,k);
		vpm1 = 0.5 * ( vtx1->position[tl_old] + vtx1->position[time_level+1] );
		vpm2 = 0.5 * ( vtx2->position[tl_old] + vtx2->position[time_level+1] );
		vpm3 = 0.5 * ( vtx3->position[tl_old] + vtx3->position[time_level+1] );
		vpm4 = 0.5 * ( vtx4->position[tl_old] + vtx4->position[time_level+1] );
		IFace->pos = 0.25 * (vpm1 + vpm2 + vpm3 + vpm4);
		IFace->vel = 0.25 * (vtx1->position[time_level+1] + vtx2->position[time_level+1] +
				     vtx3->position[time_level+1] + vtx4->position[time_level+1] - 
				     vtx1->position[tl_old] - vtx2->position[tl_old] - 
				     vtx3->position[tl_old] - vtx4->position[tl_old]) / dt;
		p1 = vpm1;
		p4 = vpm2;
		p2 = vpm3;
		p3 = vpm4;
		// Interface area at midpoint.		
		IFace->ar[time_level] = vabs(0.25 * cross(p2-p1+p3-p4, p4-p1+p3-p2)); 
	    }
	}
    }
    for (k = kmin; k <= kmax+1; ++k) {
	for (j = jmin; j <= jmax; ++j) {
	    for (i = imin; i <= imax; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i+1,j,k);
		vtx3 = get_vtx(i,j+1,k);
		vtx4 = get_vtx(i+1,j+1,k);
		IFace = get_ifk(i,j,k);
		vpm1 = 0.5 * ( vtx1->position[tl_old] + vtx1->position[time_level+1] );
		vpm2 = 0.5 * ( vtx2->position[tl_old] + vtx2->position[time_level+1] );
		vpm3 = 0.5 * ( vtx3->position[tl_old] + vtx3->position[time_level+1] );
		vpm4 = 0.5 * ( vtx4->position[tl_old] + vtx4->position[time_level+1] );
		IFace->pos = 0.25 * (vpm1 + vpm2 + vpm3 + vpm4);
		IFace->vel = 0.25 * (vtx1->position[time_level+1] + vtx2->position[time_level+1] +
				     vtx3->position[time_level+1] + vtx4->position[time_level+1] - 
				     vtx1->position[tl_old] - vtx2->position[tl_old] - 
				     vtx3->position[tl_old] - vtx4->position[tl_old]) / dt;
		p1 = vpm1;
		p4 = vpm2;
		p2 = vpm3;
		p3 = vpm4;
		// Interface area at midpoint.		
		IFace->ar[time_level] = vabs(0.25 * cross(p2-p1+p3-p4, p4-p1+p3-p2)); 
	    }
	}
    }
    return SUCCESS;
}

/// \brief Set interface velocities as average of adjacent vertex velocities. 
///
int Block::set_interface_velocities2D( int time_level )
{
    int i, j, k;
    FV_Vertex *vtx1, *vtx2;
    FV_Interface *IFace;
    k = kmin;
    for (j = jmin; j <= jmax; ++j) {
	for (i = imin; i <= imax+1; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i,j+1,k);
	    IFace = get_ifi(i,j,k);
	    IFace->vel = (vtx1->velocity[time_level] + vtx2->velocity[time_level]) / 2.0;
	}
    }
    for (j = jmin; j <= jmax+1; ++j) {
	for (i = imin; i <= imax; ++i) {
	    vtx1 = get_vtx(i,j,k);
	    vtx2 = get_vtx(i+1,j,k);
	    IFace = get_ifj(i,j,k);
	    IFace->vel = (vtx1->velocity[time_level] + vtx2->velocity[time_level]) / 2.0;
	}
    }
    return SUCCESS;
}

/// \brief Set interface velocities as average of adjacent vertex velocities. 
///
int Block::set_interface_velocities3D( int time_level )
{
    int i, j, k;
    FV_Vertex *vtx1, *vtx2, *vtx3, *vtx4;
    FV_Interface *IFace;
    
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i,j+1,k);
		vtx3 = get_vtx(i,j,k+1);
		vtx4 = get_vtx(i,j+1,k+1);
		IFace = get_ifi(i,j,k);
		IFace->vel = (vtx1->velocity[time_level] + vtx2->velocity[time_level] + 
			      vtx3->velocity[time_level] + vtx4->velocity[time_level]) / 2.0;
	    }
	}
    }
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax+1; ++j) {
	    for (i = imin; i <= imax; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i,j,k+1);
		vtx3 = get_vtx(i+1,j,k);
		vtx4 = get_vtx(i+1,j,k+1);
		IFace = get_ifj(i,j,k);
		IFace->vel = (vtx1->velocity[time_level] + vtx2->velocity[time_level] + 
			      vtx3->velocity[time_level] + vtx4->velocity[time_level]) / 2.0;
	    }
	}
    }
    for (k = kmin; k <= kmax+1; ++k) {
	for (j = jmin; j <= jmax; ++j) {
	    for (i = imin; i <= imax; ++i) {
		vtx1 = get_vtx(i,j,k);
		vtx2 = get_vtx(i+1,j,k);
		vtx3 = get_vtx(i,j+1,k);
		vtx4 = get_vtx(i+1,j+1,k);
		IFace = get_ifk(i,j,k);
		IFace->vel = (vtx1->velocity[time_level] + vtx2->velocity[time_level] + 
			      vtx3->velocity[time_level] + vtx4->velocity[time_level]) / 2.0;
	    }
	}
    }
    return SUCCESS;
}

// int Block::diffuse_vertex_velocities(double mu, int npass, int dimensions, int time_level)
// /// \brief Filter the cell-centred primary variables.
// ///
// /// This filtering is done on a block-by-block basis.
// /// Valid flow data are needed in the ghost cells at the edges.
// /// \param alpha : filter coefficient (closer to 1.0, more fudging)
// /// \param npass : this many passes of the simple averager
// //
// {
//     FV_Vertex *vtx, *vtxN, *vtxS;
//     int i,j;
//     Vector3 diffuse[nnj+5];
//     i = imin;
//     if ( bcp[SOUTH]->type_code != ADJACENT ) { // If not adjacent to another block on south side.
// 	// First
// 	get_vtx(i,jmin-1)->velocity[time_level] = get_vtx(i,jmin)->velocity[time_level];
// 	// Second
// 	get_vtx(i,jmin-2)->velocity[time_level] = get_vtx(i,jmin+1)->velocity[time_level];
//     }
//     if ( bcp[NORTH]->type_code != ADJACENT ) { // If not adjacent to another block on south side.
// 	// Last
// 	get_vtx(i,jmax+2)->velocity[time_level] = get_vtx(i,jmax+1)->velocity[time_level];
// 	// Second last
// 	get_vtx(i,jmax+3)->velocity[time_level] = get_vtx(i,jmax)->velocity[time_level];
//     }
//     for (j = jmin; j <= jmax+1; ++j) {
// 	vtx = get_vtx(i,j);
// 	vtxN = get_vtx(i,j+1);
// 	vtxS = get_vtx(i,j-1);
// 	diffuse[j] = ( 1 - 4 * mu ) * (vtx)->velocity[time_level] +
// 	    mu * ( (vtxN)->velocity[time_level] + (vtxS)->velocity[time_level]);
//     }
//      for (j = jmin; j <= jmax+1; ++j) {
// 	vtx = get_vtx(i,j);
// 	vtx->velocity[time_level] = diffuse[j];
//     }
//     return SUCCESS;
// } // end of diffuse_vertex_velocities()

// int Block::anti_diffuse_vertex_velocities(double mu, int npass, int dimensions, int time_level)
// /// \brief Filter the cell-centred primary variables.
// ///
// /// This filtering is done on a block-by-block basis.
// /// Valid flow data are needed in the ghost cells at the edges.
// /// \param alpha : filter coefficient
// /// \param npass : this many passes
// //
// {
//     FV_Vertex *vtx, *vtxN, *vtxS;
//     int i, j;
//     Vector3 m2, m1, p1, p2;
//     Vector3 vel;
//     Vector3 adflux[nnj+5];
//     // Apply the anti-diffusion.
//     i = imin;
//     for (j = jmin; j <= jmax+2; ++j) {
// 	p2 = get_vtx(i,j+1)->velocity[time_level];
// 	p1 = get_vtx(i,j)->velocity[time_level];
// 	m1 = get_vtx(i,j-1)->velocity[time_level];
// 	m2 = get_vtx(i,j-2)->velocity[time_level];
// 	vel.x = calc_anti_diffusive_flux(m2.x, m1.x, p1.x, p2.x, mu);
// 	vel.y = calc_anti_diffusive_flux(m2.y, m1.y, p1.y, p2.y, mu);
// 	vel.z = calc_anti_diffusive_flux(m2.z, m1.z, p1.z, p2.z, mu);
// 	adflux[j] = vel;
//     }
//     for (j = jmin; j <= jmax+1; ++j) {
// 	vtx = get_vtx(i,j);
// 	vtx->velocity[time_level] += adflux[j] - adflux[j+1];
//     }
//     return SUCCESS;
// } // end of anti_diffuse_vertex_velocities()

int Block::compute_boundary_flux(FV_Interface *IFaceL, FV_Interface *IFaceR, double omegaz)
{
    ConservedQuantities &F = *(IFaceR->F);
    FlowState *IFace_flow_state = IFaceL->fs;

    // Transform to interface frame of reference.
    IFace_flow_state->vel -= IFaceR->vel;
    IFaceR->vel.transform_to_local(IFaceR->n, IFaceR->t1, IFaceR->t2);
    IFace_flow_state->vel.transform_to_local(IFaceR->n, IFaceR->t1, IFaceR->t2);
    // also transform the magnetic field
    if (get_mhd_flag() == 1) {
	IFace_flow_state->B.transform_to_local(IFaceR->n, IFaceR->t1, IFaceR->t2);
    }
    set_interface_flux(*IFaceR, IFace_flow_state);
    if ( omegaz != 0.0 ) {
	// Rotating frame.
	double x = IFaceR->pos.x;
	double y = IFaceR->pos.y;
	double rsq = x*x + y*y;
	// The conserved quantity is rothalpy,
	// so we need to take -(u**2)/2 off the total energy Shock.
	// Note that rotating frame velocity u = omegaz * r.
	F.total_energy -= F.mass * 0.5*omegaz*omegaz*rsq;
    }

    // Transform fluxes back from interface frame of reference to local frame of reference.
    /* Flux of Total Energy */
    F.total_energy += 0.5 * F.mass * pow(vabs(IFaceR->vel),2) + dot(F.momentum, IFaceR->vel);
    /* Flux of momentum */
    F.momentum += F.mass * IFaceR->vel;

    // Rotate momentum fluxes back to the global frame of reference.
    F.momentum.transform_to_global(IFaceR->n, IFaceR->t1, IFaceR->t2);
    // also transform the interface velocities
    IFaceR->vel.transform_to_global(IFaceR->n, IFaceR->t1, IFaceR->t2);	  
    // also transform the magnetic field
    if (get_mhd_flag() == 1) {
	F.B.transform_to_global(IFaceR->n, IFaceR->t1, IFaceR->t2);
    }
    
    return SUCCESS;
} // end Block::compute_boundary_flux()

/** \brief Computes the (pressure and shear) forces applied by the gas 
 *         to the bounding surface.
 *
 * We make use of geometric quantities stored at 
 * the cell interfaces.  
 * (Area is area per unit radian for axisymmetric calculations.)
 */
void Block::compute_x_forces( char *text_string, int ibndy, int dimensions )
{
    double fx_p, fx_v, x1, y1, cosX, cosY, area;
    double xc, yc, d, vt, mu;
    int i, j, ivisc;
    FV_Cell *cell;
    FV_Interface *IFace;
    
    if ( dimensions == 3 ) {
	printf( "X-Force calculations not implemented for 3D geometries, yet." );
	exit( NOT_IMPLEMENTED_ERROR );
    }

    fx_p = 0.0;
    fx_v = 0.0;
    ivisc = get_viscous_flag();

    if ( ibndy == NORTH ) {
	j = jmax;
	for ( i = imin; i <= imax; ++i ) {
	    cell = get_cell(i,j);
	    IFace = cell->iface[NORTH];
	    cosX = IFace->n.x;
	    cosY = IFace->n.y;
	    area = IFace->area;
	    mu   = IFace->fs->gas->mu;
	    fx_p += cell->fs->gas->p * area * cosX;
	    if ( ivisc ) {
		/* pieces needed to reconstruct the local velocity gradient */
		x1 = cell->vtx[0]->pos.x; y1 = cell->vtx[0]->pos.y;
		xc = cell->pos.x;   yc = cell->pos.y;
		d = -(xc - x1) * cosX + -(yc - y1) * cosY;
		vt = cell->fs->vel.x * cosY - cell->fs->vel.y * cosX;
		/* x-component of the shear force, assuming no-slip wall */
		fx_v += mu * vt / d * area * cosY;
	    }
	}
    } else if ( ibndy == SOUTH ) {
	j = jmin;
	for ( i = imin; i <= imax; ++i ) {
	    cell = get_cell(i,j);
	    IFace = cell->iface[SOUTH];
	    cosX = IFace->n.x;
	    cosY = IFace->n.y;
	    area = IFace->area;
	    mu   = IFace->fs->gas->mu;
	    fx_p -= cell->fs->gas->p * area * cosX;
	    if ( ivisc ) {
		/* pieces needed to reconstruct the local velocity gradient */
		x1 = cell->vtx[0]->pos.x; y1 = cell->vtx[0]->pos.y;
		xc = cell->pos.x;   yc = cell->pos.y;
		d = (xc - x1) * cosX + (yc - y1) * cosY;
		vt = cell->fs->vel.x * cosY - cell->fs->vel.y * cosX;
		/* x-component of the shear force, assuming no-slip wall */
		fx_v += mu * vt / d * area * cosY;
	    }
	}
    } else if ( ibndy == EAST ) {
	i = imax;
	for ( j = jmin; j <= jmax; ++j ) {
	    cell = get_cell(i,j);
	    IFace = cell->iface[EAST];
	    cosX = IFace->n.x;
	    cosY = IFace->n.y;
	    area = IFace->area;
	    mu   = IFace->fs->gas->mu;
	    fx_p += cell->fs->gas->p * area * cosX;
	    if ( ivisc ) {
		/* pieces needed to reconstruct the local velocity gradient */
		x1 = cell->vtx[1]->pos.x; y1 = cell->vtx[1]->pos.y;
		xc = cell->pos.x;   yc = cell->pos.y;
		d = -(xc - x1) * cosX + -(yc - y1) * cosY;
		vt = -(cell->fs->vel.x) * cosY + cell->fs->vel.y * cosX;
		/* x-component of the shear force, assuming no-slip wall */
		fx_v -= mu * vt / d * area * cosY;
	    }
	}
    } else if ( ibndy == WEST ) {
	i = imin;
	for ( j = jmin; j <= jmax; ++j ) {
	    cell = get_cell(i,j);
	    IFace = cell->iface[WEST];
	    cosX = IFace->n.x;
	    cosY = IFace->n.y;
	    area = IFace->area;
	    mu   = IFace->fs->gas->mu;
	    fx_p -= cell->fs->gas->p * area * cosX;
	    if ( ivisc ) {
		/* pieces needed to reconstruct the local velocity gradient */
		x1 = cell->vtx[0]->pos.x; y1 = cell->vtx[0]->pos.y;
		xc = get_cell(i,j)->pos.x;  yc = get_cell(i,j)->pos.y;
		d = (xc - x1) * cosX + (yc - y1) * cosY;
		vt = -(cell->fs->vel.x) * cosY + cell->fs->vel.y * cosX;
		/* x-component of the shear force, assuming no-slip wall */
		fx_v -= mu * vt / d * area * cosY;
	    }
	}
    }   /* end if: boundary selection */

    if ( get_axisymmetric_flag() == 1 ) {
	fx_p *= (2.0 * 3.1415927);
	fx_v *= (2.0 * 3.1415927);
    }

    sprintf( text_string, "FX_P %e FX_V %e ", fx_p, fx_v );
} // end compute_x_forces()


/// \brief Assemble the x-force numbers for each block into a single
///        (string) report and send it to the logfile. 
int Block::print_forces( FILE *fp, double t, int dimensions )
{
    char msg_text[512], small_text[132];

    if ( bcp[NORTH]->xforce_flag == 1 ) {
	sprintf( small_text, "XFORCE: TIME %e BLOCK %d BNDY %d ", t, id, NORTH );
	strcpy( msg_text, small_text );
	this->compute_x_forces( small_text, NORTH, dimensions );
	strcat( msg_text, small_text );
	fprintf( fp, "%s\n",  msg_text );
    }
    if ( bcp[EAST]->xforce_flag == 1 ) {
	sprintf( small_text, "XFORCE: TIME %e BLOCK %d BNDY %d ", t, id, EAST );
	strcpy( msg_text, small_text );
	this->compute_x_forces( small_text, EAST, dimensions );
	strcat( msg_text, small_text );
	fprintf( fp, "%s\n",  msg_text );
    }
    if ( bcp[SOUTH]->xforce_flag == 1 ) {
	sprintf( small_text, "XFORCE: TIME %e BLOCK %d BNDY %d ", t, id, SOUTH );
	strcpy( msg_text, small_text );
	this->compute_x_forces( small_text, SOUTH, dimensions );
	strcat( msg_text, small_text );
	fprintf( fp, "%s\n",  msg_text );
    }
    if ( bcp[WEST]->xforce_flag == 1 ) {
	sprintf( small_text, "XFORCE: TIME %e BLOCK %d BNDY %d ", t, id, WEST );
	strcpy( msg_text, small_text );
	this->compute_x_forces( small_text, WEST, dimensions );
	strcat( msg_text, small_text );
	fprintf( fp, "%s\n",  msg_text );
    }

    return SUCCESS;
} // end print_forces()


int Block::read_grid(std::string filename, int dimensions, int zip_file)
/// \brief Read the grid from a disc file as a set of cell vertices.
/// \returns 0 if successful but 1 if it hits the end of the grid file prematurely.
{
#   define NCHAR 132
    char line[NCHAR];
    char *gets_result;
    FV_Vertex *vp;
    int i, j, k;
    int retries = 10;
    FILE *fp = NULL;
    gzFile zfp = NULL;
    if (id == 0) printf("read_grid(): Start block %d.\n", id);
    retries = 10;
    if (zip_file) filename += ".gz";
    while (retries > 0 && zfp == NULL && fp == NULL) {
	if (zip_file) {
	    zfp = gzopen(filename.c_str(), "r");
	} else {
	    fp = fopen(filename.c_str(), "r");
	}
	if (zfp == NULL && fp == NULL) {
	    --retries;
	    cerr << "read_grid(): Could not open " << filename 
		 << "; " << retries << " retries to go." << endl;
	    sleep(2);
	}
    }
    if (zfp == NULL && fp == NULL) {
	cerr << "read_grid(): Could not open " << filename << "; BAILING OUT" << endl;
	return FILE_ERROR;
    }
    if (zip_file) {
	gets_result = gzgets(zfp, line, NCHAR);
    } else {
	gets_result = fgets(line, NCHAR, fp);
    }
    if (gets_result == NULL) {
	printf("read_grid(): Empty grid file, block %d.\n", id);
	return BAD_INPUT_ERROR;
    }
    sscanf(line, "%d %d %d", &i, &j, &k);
    if (dimensions == 3) {
	if ( i != nni+1 || j != nnj+1 || k != nnk+1 ) {
	    printf("read_grid(): Mismatch in cell numbers, block %d\n", id);
	    printf("    i=%d nni+1=%d j=%d nnj+1=%d k=%d nnk+1=%d\n", 
		   i, nni+1, j, nnj+1, k, nnk+1);
	    return BAD_INPUT_ERROR;
	}
	for ( k = kmin; k <= kmax+1; ++k ) {
	    for ( j = jmin; j <= jmax+1; ++j ) {
		for ( i = imin; i <= imax+1; ++i ) {
		    if (zip_file) {
			gets_result = gzgets(zfp, line, NCHAR);
		    } else {
			gets_result = fgets(line, NCHAR, fp);
		    }
		    if (gets_result == NULL) {
			printf("read_grid(): Premature end of file, block %d, vertex[%d,%d,%d]\n", id, i, j, k);
			return BAD_INPUT_ERROR;
		    }
		    vp = get_vtx(i,j,k);
		    sscanf(line, "%lf %lf %lf", &(vp->pos.x), &(vp->pos.y), &(vp->pos.z));
		} // for i
	    } // for j
	} // for k
    } else {
	// 2-dimensional case.
	if ( i != nni+1 || j != nnj+1 || k != 1 ) {
	    printf( "read_grid(): Mismatch in cell numbers, block %d\n", id );
	    printf( "    i=%d nni+1=%d j=%d nnj+1=%d k=%d nnk=%d\n", 
		    i, nni+1, j, nnj+1, k, nnk);
	    return BAD_INPUT_ERROR;
	}
	for ( j = jmin; j <= jmax+1; ++j ) {
	    for ( i = imin; i <= imax+1; ++i ) {
		if (zip_file) {
		    gets_result = gzgets(zfp, line, NCHAR);
		} else {
		    gets_result = fgets(line, NCHAR, fp);
		}
		if (gets_result == NULL) {
		    printf("read_grid(): Premature end of file, block %d, vertex[%d,%d]\n", id, i, j);
		    return BAD_INPUT_ERROR;
		}
		vp = get_vtx(i,j);
		sscanf(line, "%lf %lf", &(vp->pos.x), &(vp->pos.y));
		vp->pos.z = 0.0;
	    } // for i
	} // for j
    }
    if (zip_file) {
	gzclose(zfp);
    } else {
	fclose(fp);
    }
    return SUCCESS;
#   undef NCHAR
} /* end of Block::read_grid() */


/// \brief Read the flow solution (i.e. the primary variables at the 
///        cell centers) from a disk file.
/// Returns a status flag.
int Block::read_solution(std::string filename, double *sim_time,
			 int dimensions, int zip_file)
{
#   define NCHAR 4000
    char line[NCHAR];
    char *gets_result;
    int i, j, k;
    int retries = 10;
    FILE *fp = NULL;
    gzFile zfp = NULL;
    if (id == 0) printf("read_solution(): Start block %d.\n", id); 
    if (zip_file) filename += ".gz";
    while (retries > 0 && zfp == NULL && fp == NULL) {
	if (zip_file) {
	    zfp = gzopen(filename.c_str(), "r");
	} else {
	    fp = fopen(filename.c_str(), "r");
	}
	if (zfp == NULL && fp == NULL) {
	    --retries;
	    cerr << "read_solution(): Could not open " << filename 
		 << "; " << retries << " retries to go." << endl;
	    sleep(2);
	}
    }
    if (zfp == NULL && fp == NULL) {
	cerr << "read_solution(): Could not open " << filename << "; BAILING OUT" << endl;
	return FILE_ERROR;
    }

    cout << "read_solution(): " << filename << ";" << endl;

    if (zip_file) {
	gets_result = gzgets(zfp, line, NCHAR);
    } else {
	gets_result = fgets(line, NCHAR, fp);
    }
    if (gets_result == NULL) {
	printf("read_solution(): Empty flow field file while looking for sim_time value.\n");
	return BAD_INPUT_ERROR;
    }
    sscanf(line, "%lf", sim_time);
    if (id == 0) printf("read_solution(): Time = %e\n", *sim_time);
    if (zip_file) {
	gets_result = gzgets(zfp, line, NCHAR);
    } else {
	gets_result = fgets(line, NCHAR, fp);
    }
    if (gets_result == NULL) {
	printf("read_solution(): Empty flow field file while looking for line of variable names.\n");
	return BAD_INPUT_ERROR;
    }
    // The line just read should be the list of variable names, double-quoted.
    if (zip_file) {
	gets_result = gzgets(zfp, line, NCHAR);
    } else {
	gets_result = fgets(line, NCHAR, fp);
    }
    if ( gets_result == NULL ) {
	printf("read_solution(): Empty flow field file while looking for numbers of cells.\n");
	return BAD_INPUT_ERROR;
    }
    sscanf(line, "%d %d %d", &i, &j, &k);
    if ( i != nni || j != nnj || k != ((dimensions == 3) ? nnk : 1) ) {
	printf("read_solution(): block %d, mismatch in cell numbers\n", id);
	printf("    This misalignment could be caused by a having a different number\n");
	printf("    of fields for each cell's entry.\n");
	return BAD_INPUT_ERROR;
    }
    for ( k = kmin; k <= kmax; ++k ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( i = imin; i <= imax; ++i ) {
		// The new format for Elmer3 puts all cell data onto one line.
		if (zip_file) {
		    gets_result = gzgets(zfp, line, NCHAR);
		} else {
		    gets_result = fgets(line, NCHAR, fp);
		}
		if (gets_result == NULL) {
		    printf("read_solution(): Empty flow field file while reading cell data.\n");
		    return BAD_INPUT_ERROR;
		}
		get_cell(i,j,k)->scan_values_from_string(line);
	    }
	}
    }
    if (zip_file) {
	gzclose(zfp);
    } else {
	fclose(fp);
    }
    return SUCCESS;
#   undef NCHAR
} // end of Block::read_solution()


int Block::write_solution( std::string filename, double sim_time, int dimensions, int zip_file )
/// \brief Write the flow solution (i.e. the primary variables at the
///        cell centers) for a single block.
///
/// This is "almost-Tecplot" POINT format.
{
    FILE *fp;
    gzFile zfp;
    string str;
    if (id == 0) {
	printf("write_solution(): At t = %e, start block = %d.\n", sim_time, id);
    }
    if (zip_file) {
	fp = NULL;
	filename += ".gz";
	if ((zfp = gzopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_solution(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	gzprintf(zfp, "%20.12e\n", sim_time);
	gzprintf(zfp, "%s\n", variable_list_for_cell().c_str());
	gzprintf(zfp, "%d %d %d\n", nni, nnj, nnk);
    } else {
	zfp = NULL;
	if ((fp = fopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_solution(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	fprintf(fp, "%20.12e\n", sim_time);
	fprintf(fp, "%s\n", variable_list_for_cell().c_str());
	fprintf(fp, "%d %d %d\n", nni, nnj, nnk);
    }
    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		str = get_cell(i,j,k)->write_values_to_string();
		if (zip_file) {
		    gzputs(zfp, str.c_str()); gzputc(zfp, '\n');
		} else {
		    fputs(str.c_str(), fp); fputc('\n', fp);
		}
	    } // i-loop
	} // j-loop
    } // k-loop
    if (zip_file) {
	gzclose(zfp);
    } else {
	fclose(fp);
    }
    return SUCCESS;
} // end of Block::write_solution()

int Block::write_block( std::string filename, double sim_time, int dimensions, int zip_file )
/// \brief Write the flow solution (i.e. the primary variables at the
///        cell centers) for a single block.
///
/// This is "almost-Tecplot" POINT format.
{
    FV_Vertex *vtx;
    FILE *fp;
    gzFile zfp;
    string str;
    int krangemax;
    if (id == 0) {
	printf("write_block(): At t = %e, start block = %d.\n", sim_time, id);
    }
    if (zip_file) {
	fp = NULL;
	filename += ".gz";
	if ((zfp = gzopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_block(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	if ( dimensions == 2 ) gzprintf(zfp, "%d %d %d  # ni nj nk\n", nni+1, nnj+1, nnk);
	else gzprintf(zfp, "%d %d %d  # ni nj nk\n", nni+1, nnj+1, nnk+1);
    } else {
	zfp = NULL;
	if ((fp = fopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_block(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	if ( dimensions == 2 ) fprintf(fp, "%d %d %d  # ni nj nk\n", nni+1, nnj+1, nnk);
	else fprintf(fp, "%d %d %d  # ni nj nk\n", nni+1, nnj+1, nnk+1);
    }
    if ( dimensions == 2 ) krangemax = kmax;
    else krangemax = kmax+1;
    for ( int k = kmin; k <= krangemax; ++k ) {
	for ( int j = jmin; j <= jmax+1; ++j ) {
	    for ( int i = imin; i <= imax+1; ++i ) {
		vtx = get_vtx(i,j,k);
		if (zip_file) {
		    gzprintf(zfp, "%20.12e %20.12e %20.12e\n", vtx->pos.x, vtx->pos.y, vtx->pos.z);
		} else {
		    fprintf(fp, "%20.12e %20.12e %20.12e\n", vtx->pos.x, vtx->pos.y, vtx->pos.z);
		}
	    } // i-loop
	} // j-loop
    } // k-loop
    if (zip_file) {
	gzclose(zfp);
    } else {
	fclose(fp);
    }
    return SUCCESS;
} // end of Block::write_block()

/// \brief Read the BGK discrete velocity distribution values
/// Returns a status flag.
int Block::read_BGK(std::string filename, double *sim_time,
		    int dimensions, int zip_file)
{
#   define NCHAR 4000
    char line[NCHAR];
    char *gets_result;
    int i, j, k;
    int retries = 10;
    FILE *fp = NULL;
    gzFile zfp = NULL;
    if (id == 0) printf("read_BGK(): Start block %d.\n", id); 
    if (zip_file) filename += ".gz";
    while (retries > 0 && zfp == NULL && fp == NULL) {
	if (zip_file) {
	    zfp = gzopen(filename.c_str(), "r");
	} else {
	    fp = fopen(filename.c_str(), "r");
	}
	if (zfp == NULL && fp == NULL) {
	    --retries;
	    cerr << "read_BGK(): Could not open " << filename 
		 << "; " << retries << " retries to go." << endl;
	    sleep(2);
	}
    }
    if (zfp == NULL && fp == NULL) {
	cerr << "read_BGK(): Could not open " << filename << "; BAILING OUT" << endl;
	return FILE_ERROR;
    }
    if (zip_file) {
	gets_result = gzgets(zfp, line, NCHAR);
    } else {
	gets_result = fgets(line, NCHAR, fp);
    }
    if (gets_result == NULL) {
	printf("read_BGK(): Empty flow field file while looking for sim_time value.\n");
	return BAD_INPUT_ERROR;
    }
    sscanf(line, "%lf", sim_time);
    if (id == 0) printf("read_BGK(): Time = %e\n", *sim_time);
    if (zip_file) {
	gets_result = gzgets(zfp, line, NCHAR);
    } else {
	gets_result = fgets(line, NCHAR, fp);
    }
    if (gets_result == NULL) {
	printf("read_BGK(): Empty flow field file while looking for line of variable names.\n");
	return BAD_INPUT_ERROR;
    }
    // The line just read should be the list of variable names, double-quoted.
    if (zip_file) {
	gets_result = gzgets(zfp, line, NCHAR);
    } else {
	gets_result = fgets(line, NCHAR, fp);
    }
    if ( gets_result == NULL ) {
	printf("read_BGK(): Empty flow field file while looking for numbers of cells.\n");
	return BAD_INPUT_ERROR;
    }
    sscanf(line, "%d %d %d", &i, &j, &k);
    if ( i != nni || j != nnj || k != ((dimensions == 3) ? nnk : 1) ) {
	printf("read_BGK(): block %d, mismatch in cell numbers\n", id);
	printf("    This misalignment could be caused by a having a different number\n");
	printf("    of fields for each cell's entry.\n");
	return BAD_INPUT_ERROR;
    }
    for ( k = kmin; k <= kmax; ++k ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    for ( i = imin; i <= imax; ++i ) {
		// The new format for Elmer3 puts all cell data onto one line.
		if (zip_file) {
		    gets_result = gzgets(zfp, line, NCHAR);
		} else {
		    gets_result = fgets(line, NCHAR, fp);
		}
		if (gets_result == NULL) {
		    printf("read_BGK(): Empty flow field file while reading cell data.\n");
		    return BAD_INPUT_ERROR;
		}
		get_cell(i,j,k)->scan_BGK_from_string(line);
	    }
	}
    }
    if (zip_file) {
	gzclose(zfp);
    } else {
	fclose(fp);
    }
    return SUCCESS;
#   undef NCHAR
} // end of Block::read_BGK()

/// \brief Initialise the discrete samples of the BGK velocity distribution function 
/// from the macroscopic variables already loaded into each cell assuming equilibrium
// conditions in the cell.
int Block::initialise_BGK_equilibrium( void )
{
    for (int k = kmin; k <= kmax; ++k ) {
	for (int j = jmin; j <= jmax; ++j ) {
	    for (int i = imin; i <= imax; ++i ) {
		get_cell(i,j,k)->fs->BGK_equilibrium();
	    }
	}
    }
    return SUCCESS;
} // end of Block::initialise_BGK_equilibrium()

/// \brief Write the BGK discrete velocity distribution values for a single block
///
/// This is "almost-Tecplot" POINT format.
int Block::write_BGK( std::string filename, double sim_time, int dimensions, int zip_file )
{
    FILE *fp;
    gzFile zfp;
    string str;
    if (id == 0) {
	printf("write_BGK(): At t = %e, start block = %d.\n", sim_time, id);
    }
    if (zip_file) {
	fp = NULL;
	filename += ".gz";
	if ((zfp = gzopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_BGK(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	gzprintf(zfp, "%20.12e\n", sim_time);
	gzprintf(zfp,"\"pos.x\" \"pos.y\" \"pos.z\" \"volume\"");
	for (int ii = 0; ii < get_velocity_buckets(); ++ii) {
	    gzprintf(zfp," \"G[%d]\" \"H[%d]\"",ii,ii);
	}
	gzprintf(zfp, "\n%d %d %d\n", nni, nnj, nnk);
    } else {
	zfp = NULL;
	if ((fp = fopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_BGK(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	fprintf(fp, "%20.12e\n", sim_time);
	fprintf(fp,"\"pos.x\" \"pos.y\" \"pos.z\" \"volume\"");
	
	for (int ii = 0; ii < get_velocity_buckets(); ++ii) {
	    fprintf(fp," \"G[%d]\" \"H[%d]\"",ii,ii);
	}
	fprintf(fp, "\n%d %d %d\n", nni, nnj, nnk);
    }
    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int j = jmin; j <= jmax; ++j ) {
	    for ( int i = imin; i <= imax; ++i ) {
		str = get_cell(i,j,k)->write_BGK_to_string();
		if (zip_file) {
		    gzputs(zfp, str.c_str()); gzputc(zfp, '\n');
		} else {
		    fputs(str.c_str(), fp); fputc('\n', fp);
		}
	    } // i-loop
	} // j-loop
    } // k-loop
    if (zip_file) {
	gzclose(zfp);
    } else {
	fclose(fp);
    }
    return SUCCESS;
} // end of Block::write_BGK()


int Block::write_history( std::string filename, double sim_time, int write_header )
/// \brief Write out the flow solution in a (small) subset of cells.
///
/// This us usually done at a different (often smaller) time interval 
/// to the full flow solution.
/// Note that, after writing the header, the file is opened in append mode 
/// so that the data may accumulate.
{
    int i, j, k;
    FILE *fp;
    string str;
    if ( write_header ) {
	if ((fp = fopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_history(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	fprintf(fp, "# \"time\" \"i\" \"j\" \"k\" ");
	fprintf(fp, "%s\n", variable_list_for_cell().c_str());
    } else {
	if ((fp = fopen(filename.c_str(), "a")) == NULL) {
	    cerr << "write_history(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	for ( int ih = 0; ih < hncell; ++ih ) {
	    i = hicell[ih] + imin;
	    j = hjcell[ih] + jmin;
	    k = hkcell[ih] + kmin;
	    fprintf( fp, "%e %d %d %d ", sim_time, hicell[ih], hjcell[ih], hkcell[ih]);
	    str = get_cell(i,j,k)->write_values_to_string();
	    fputs( str.c_str(), fp );
	    fputc( '\n', fp );
	}
    }
    fclose(fp);
    return SUCCESS;
} // end of Block::write_history()


int Block::count_invalid_cells( int dimensions )
/// \brief Returns the number of cells that contain invalid data.
///
/// This data can be identified by the density of internal energy 
/// being on the minimum limit or the velocity being very large.
//
// To do: We should probably make this function more 3D friendly, however,
// it should not be invoked (ever) if the code is working well!
{
    FV_Cell *cell;
    int  number_of_invalid_cells;
    Gas_model *gmodel = get_gas_model_ptr();
    FV_Cell dummy_cell(gmodel);
    
    number_of_invalid_cells = 0;

    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int i = imin; i <= imax; ++i ) {
	    for ( int j = jmin; j <= jmax; ++j ) {
		cell = get_cell(i,j,k);
		if ( cell->check_flow_data() != 1 ) {
		    ++number_of_invalid_cells;
		    if ( get_bad_cell_complain_flag() ) {
			printf("count_invalid_cells: block_id = %d, cell[%d,%d,%d]\n", 
			       id, i, j, k);
			cell->print();
		    }
		    
#                   if ADJUST_INVALID_CELL_DATA == 1
		    // We shall set the cell data to something that
		    // is valid (and self consistent).
		    FV_Cell *other_cell;
#                   if PREFER_COPY_FROM_LEFT == 1
		    if ( get_bad_cell_complain_flag() ) {
			printf( "Adjusting cell data by copying data from left.\n" );
		    }
		    other_cell = get_cell(i-1,j,k);
		    cell->copy_values_from(other_cell, COPY_FLOW_STATE);
#                   else
		    if ( get_bad_cell_complain_flag() ) {
			printf( "Adjusting cell data to a local average.\n" );
		    }
		    // Presently, only searches around in the i,j plane
		    std::vector<FV_Cell *> neighbours;
		    other_cell = get_cell(i-1,j,k);
		    if ( other_cell->check_flow_data() ) neighbours.push_back(other_cell);
		    other_cell = get_cell(i+1,j,k);
		    if ( other_cell->check_flow_data() ) neighbours.push_back(other_cell);
		    other_cell = get_cell(i,j-1,k);
		    if ( other_cell->check_flow_data() ) neighbours.push_back(other_cell);
		    other_cell = get_cell(i,j+1,k);
		    if ( other_cell->check_flow_data() ) neighbours.push_back(other_cell);
		    if ( neighbours.size() == 0 ) {
			printf( "It seems that there were no valid neighbours, I give up.\n" );
			exit( BAD_CELLS_ERROR );
		    }
		    cell->replace_flow_data_with_average(neighbours);
#                   endif
		    cell->encode_conserved(omegaz);
		    cell->decode_conserved(omegaz);
		    if ( get_bad_cell_complain_flag() ) {
			printf("after flow-data replacement: block_id = %d, cell[%d,%d,%d]\n", 
			       id, i, j, k);
			cell->print();
		    }
#                   endif
		} // end of if invalid data...
	    } // j loop
	} // i loop
    } // k loop
    return number_of_invalid_cells;
} // end of count_invalid_cells()


int Block::init_residuals( int dimensions )
/// \brief Initialization of data for later computing residuals.
{
    FV_Cell *cellp;
    mass_residual = 0.0;
    mass_residual_loc = Vector3(0.0, 0.0, 0.0);
    energy_residual = 0.0;
    energy_residual_loc = Vector3(0.0, 0.0, 0.0);
    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int j = jmin; j <= jmax; ++j ) {
	    for (int  i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		cellp->rho_at_start_of_step = cellp->fs->gas->rho;
		cellp->rE_at_start_of_step = cellp->U->total_energy;
	    } // for i
	} // for j
    } // for k
    return SUCCESS;
} // end of check_residual()


int Block::compute_residuals( int dimensions )
/// \brief Compute the residuals using previously stored data.
///
/// The largest residual of density for all cells was the traditional way
/// mbcns/Elmer estimated the approach to steady state.
/// However, with the splitting up of the increments for different physical
/// processes, this residual calculation code needed a bit of an update.
/// Noting that the viscous-stress, chemical and radiation increments
/// do not affect the mass within a cell, we now compute the residuals 
/// for both mass and (total) energy for all cells, the record the largest
/// with their location. 
{
    double local_residual;
    FV_Cell *cellp;
    mass_residual = 0.0;
    mass_residual_loc = Vector3(0.0, 0.0, 0.0);
    energy_residual = 0.0;
    energy_residual_loc = Vector3(0.0, 0.0, 0.0);
    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int j = jmin; j <= jmax; ++j ) {
	    for (int  i = imin; i <= imax; ++i ) {
		cellp = get_cell(i,j,k);
		// printf ( "blk[%d] cell(%d,%d,%d) rho %g rho_old %g\n",
		// 	 id, i, j, k, cellp->gas->U.rho, cellp->U_old->rho );
		local_residual = ( cellp->fs->gas->rho - cellp->rho_at_start_of_step ) / cellp->fs->gas->rho;
		local_residual = FABS(local_residual);
		if ( local_residual > mass_residual ) {
		    mass_residual = local_residual;
		    mass_residual_loc.x = cellp->pos.x;
		    mass_residual_loc.y = cellp->pos.y;
		    mass_residual_loc.z = cellp->pos.z;
		}
		local_residual = ( cellp->U->total_energy - cellp->rE_at_start_of_step ) / cellp->U->total_energy;
		local_residual = FABS(local_residual);
		if ( local_residual > energy_residual ) {
		    energy_residual = local_residual;
		    energy_residual_loc.x = cellp->pos.x;
		    energy_residual_loc.y = cellp->pos.y;
		    energy_residual_loc.z = cellp->pos.z;
		}
	    } // for i
	} // for j
    } // for k
    return SUCCESS;
} // end of compute_residuals()


int Block::determine_time_step_size( double cfl_target, int dimensions )
/// \brief Compute the local time step limit for all cells in the block.
///
/// The overall time step is limited by the worst-case cell.
/// \returns 0 on success, DT_SEARCH_FAILED otherwise.
///
/// \verbatim
/// Some Definitions...
/// ----------------
/// dt_global  : global time step for the block
/// cfl_target : desired CFL number
/// cfl_min  : approximate minimum CFL number in the block
/// cfl_max  : approximate maximum CFL number in the block
/// dt_allow : allowable time step (i.e. the maximum dt that
///            satisfies both the CFL target and the viscous
///            time step limit)
/// cfl_allow : allowable CFL number, t_order dependent
/// \endverbatim
{
    global_data *gdp = get_global_data_ptr();
    bool first;
    double dt_local, cfl_local, signal, cfl_allow;
    if ( get_Torder_flag() < 3 ) {
	cfl_allow = 0.9;
    } else {
	cfl_allow = 1.1;
    }

    first = true;
    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int i = imin; i <= imax; ++i ) {
	    for ( int j = jmin; j <= jmax; ++j ) {
		signal = get_cell(i,j,k)->signal_frequency(dimensions);
		cfl_local = gdp->dt_global * signal; // Current (Local) CFL number
		dt_local = cfl_target / signal; // Recommend a time step size.
		if ( first ) {
		    cfl_min = cfl_local;
		    cfl_max = cfl_local;
		    dt_allow = dt_local;
		    first = false;
		} else {
		    if (cfl_local < cfl_min) cfl_min = cfl_local;
		    if (cfl_local > cfl_max) cfl_max = cfl_local;
		    if (dt_local < dt_allow) dt_allow = dt_local;
		}
	    } // for j
	} // for i
    } // for k
    if ( cfl_max > 0.0 && cfl_max < cfl_allow ) {
	return SUCCESS;
    } else {
	printf( "determine_time_step_size(): bad CFL number was encountered\n" );
	printf( "    cfl_max = %e for Block %d\n", cfl_max, id );
	printf( "    If this cfl_max value is not much larger than 1.0,\n" );
	printf( "    your simulation could probably be restarted successfully\n" );
	printf( "    with some minor tweaking." );
	printf( "    That tweaking should probably include a reduction\n");
	printf( "    in the size of the initial time-step, dt\n");
	printf( "    If this job is a restart/continuation of an old job, look in\n");
	printf( "    the old-job.finish file for the value of dt at termination.\n");
	return DT_SEARCH_FAILED;
    }
} // end of determine_time_step_size()


int Block::detect_shock_points( int dimensions )
/// \brief Detects shocks by looking for compression between adjacent cells.
///
/// The velocity component normal to the cell interfaces
/// is used as the indicating variable.
{
    FV_Cell *cL, *cR, *cell;
    FV_Interface *IFace;
    double uL, uR, aL, aR, a_min;

    // Change in normalised velocity to indicate a shock.
    // A value of -0.05 has been found suitable to detect the levels of
    // shock compression observed in the "sod" and "cone20" test cases.
    // It may need to be tuned for other situations, especially when
    // viscous effects are important.
    double tol = get_compression_tolerance();

    // First, work across North interfaces and
    // locate shocks using the (local) normal velocity.
    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int i = imin; i <= imax; ++i ) {
	    for ( int j = jmin-1; j <= jmax; ++j ) {
		cL = get_cell(i,j,k);
		cR = get_cell(i,j+1,k);
		IFace = cL->iface[NORTH];
		uL = cL->fs->vel.x * IFace->n.x + cL->fs->vel.y * IFace->n.y + cL->fs->vel.z * IFace->n.z;
		uR = cR->fs->vel.x * IFace->n.x + cR->fs->vel.y * IFace->n.y + cR->fs->vel.z * IFace->n.z;
		aL = cL->fs->gas->a;
		aR = cR->fs->gas->a;
		if (aL < aR)
		    a_min = aL;
		else
		    a_min = aR;
		IFace->fs->S = ((uR - uL) / a_min < tol);
	    } // j loop
	} // i loop
    } // for k
    
    // Second, work across East interfaces and
    // locate shocks using the (local) normal velocity.
    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int i = imin-1; i <= imax; ++i ) {
	    for ( int j = jmin; j <= jmax; ++j ) {
		cL = get_cell(i,j,k);
		cR = get_cell(i+1,j,k);
		IFace = cL->iface[EAST];
		uL = cL->fs->vel.x * IFace->n.x + cL->fs->vel.y * IFace->n.y + cL->fs->vel.z * IFace->n.z;
		uR = cR->fs->vel.x * IFace->n.x + cR->fs->vel.y * IFace->n.y + cR->fs->vel.z * IFace->n.z;
		aL = cL->fs->gas->a;
		aR = cR->fs->gas->a;
		if (aL < aR)
		    a_min = aL;
		else
		    a_min = aR;
		IFace->fs->S = ((uR - uL) / a_min < tol);
	    } // j loop
	} // i loop
    } // for k
    
    if ( dimensions == 3 ) {
	// Third, work across Top interfaces.
	for ( int i = imin; i <= imax; ++i ) {
	    for ( int j = jmin; j <= jmax; ++j ) {
		for ( int k = kmin-1; k <= kmax; ++k ) {
		    cL = get_cell(i,j,k);
		    cR = get_cell(i,j,k+1);
		    IFace = cL->iface[TOP];
		    uL = cL->fs->vel.x * IFace->n.x + cL->fs->vel.y * IFace->n.y + cL->fs->vel.z * IFace->n.z;
		    uR = cR->fs->vel.x * IFace->n.x + cR->fs->vel.y * IFace->n.y + cR->fs->vel.z * IFace->n.z;
		    aL = cL->fs->gas->a;
		    aR = cR->fs->gas->a;
		    if (aL < aR)
			a_min = aL;
		    else
			a_min = aR;
		    IFace->fs->S = ((uR - uL) / a_min < tol);
		} // for k
	    } // j loop
	} // i loop
    } // if ( dimensions == 3 )
    
    // Finally, mark cells as shock points if any of their
    // interfaces are shock points.
    for ( int k = kmin; k <= kmax; ++k ) {
	for ( int i = imin; i <= imax; ++i ) {
	    for ( int j = jmin; j <= jmax; ++j ) {
		cell = get_cell(i,j,k);
		cell->fs->S = cell->iface[EAST]->fs->S || cell->iface[WEST]->fs->S ||
		    cell->iface[NORTH]->fs->S || cell->iface[SOUTH]->fs->S ||
		    ( dimensions == 3 && (cell->iface[BOTTOM]->fs->S || cell->iface[TOP]->fs->S) );
	    } // j loop
	} // i loop
    } // for k
    
    return SUCCESS;
} // end of detect_shock_points()


/// \brief Part of the spatial filter function.
/// 
/// \param PROP : name of the variable to be filtered.
#define SIMPLE_FILTER(PROP)						\
    for (i = imin; i <= imax; ++i) {					\
        for (j = jmin; j <= jmax; ++j) {				\
            cell = get_cell(i,j);					\
	    cN = get_cell(i,j+1);					\
	    cE = get_cell(i+1,j);					\
	    cS = get_cell(i,j-1);					\
	    cW = get_cell(i-1,j);					\
            (cell)->fs->PROP = (1.0 - alpha) * (cell)->fs->PROP +	\
		alpha * 0.25 * ((cN)->fs->PROP + (cE)->fs->PROP + (cS)->fs->PROP + (cW)->fs->PROP); \
        }								\
    }

/// \brief Flux Corrected Transport filter. 
/// Applies diffusion and anti-diffusion step
/// as a filter.
/// Part of the spatial filter function. Based on 
/// technical report by M. N. Macrossan, in turn based
/// on presentation by Kaplan and Oran.
/// 
/// \param PROP : name of the variable to be filtered.
#define FCT_DIFFUSION(PROP)				\
    CALC_DIFFUSION(PROP)				\
    APPLY_DIFFUSION(PROP)				\

#define FCT_ANTI_DIFFUSION(PROP)			\
    CALC_ANTI_DIFFUSION(PROP)				\
    APPLY_ANTI_DIFFUSION(PROP)				\

#define CALC_DIFFUSION(PROP)				\
    for (i = imin; i <= imax; ++i) {			\
        for (j = jmin; j <= jmax; ++j) {		\
            cell = get_cell(i,j);			\
	    cN = get_cell(i,j+1);			\
	    cE = get_cell(i+1,j);			\
	    cS = get_cell(i,j-1);			\
	    cW = get_cell(i-1,j);			\
            diffuse[i][j] = ( 1 - 4 * mu ) * (cell)->fs->PROP + \
		mu * ( (cN)->fs->PROP + (cE)->fs->PROP + (cS)->fs->PROP + (cW)->fs->PROP ); \
	}						\
    }

#define APPLY_DIFFUSION(PROP)				\
    for (i = imin; i <= imax; ++i) {			\
        for (j = jmin; j <= jmax; ++j) {		\
            cell = get_cell(i,j);			\
            (cell)->fs->PROP = diffuse[i][j];		\
	}						\
    }

#define CALC_ANTI_DIFFUSION(PROP)			\
    for (i = imin; i <= imax+1; ++i) {			\
        for (j = jmin; j <= jmax; ++j) {		\
	    p2 = get_cell(i+1,j)->fs->PROP;		\
	    p1 = get_cell(i,j)->fs->PROP;		\
	    m1 = get_cell(i-1,j)->fs->PROP;		\
	    m2 = get_cell(i-2,j)->fs->PROP;		\
            iadflux[i][j] = calc_anti_diffusive_flux(m2, m1,		\
						     p1, p2, mu);	\
	}						\
    }							\
    for (i = imin; i <= imax; ++i) {			\
        for (j = jmin; j <= jmax+1; ++j) {		\
	    p2 = get_cell(i,j+1)->fs->PROP;		\
	    p1 = get_cell(i,j)->fs->PROP;		\
	    m1 = get_cell(i,j-1)->fs->PROP;		\
	    m2 = get_cell(i,j-2)->fs->PROP;		\
            jadflux[i][j] = calc_anti_diffusive_flux(m2, m1,		\
						     p1, p2, mu);	\
	}								\
    }

#define APPLY_ANTI_DIFFUSION(PROP)			\
    for (i = imin; i <= imax; ++i) {			\
        for (j = jmin; j <= jmax; ++j) {		\
            cell = get_cell(i,j);			\
	    (cell)->fs->PROP += iadflux[i][j] + jadflux[i][j] - \
		iadflux[i+1][j] - jadflux[i][j+1];		\
	}							\
    }

double Block::calc_anti_diffusive_flux(double m2, double m1, double p1, double p2, double mu)
{
    double delm = m1 - m2;
    double del =  p1 - m1;
    double delp = p2 - p1;
    double S = copysign(1.0, del);
    double delmin = min(S * delm, S * delp);

    return S * max(0.0, min(mu * fabs(del), delmin));
}

int Block::apply_spatial_filter_diffusion(double mu, int npass, int dimensions)
/// \brief Filter the cell-centred primary variables.
///
/// This filtering is done on a block-by-block basis.
/// Valid flow data are needed in the ghost cells at the edges.
/// \param alpha : filter coefficient (closer to 1.0, more fudging)
/// \param npass : this many passes of the simple averager
//
// To do: We should fix for 3D or remove. I think that there are no cases that need it. 
{
    FV_Cell *cell, *cN, *cE, *cS, *cW;
    int isp;
    int i, j;
    Gas_model *gm = get_gas_model_ptr();
    int nsp = gm->get_number_of_species();
    vector<vector<double> > diffuse;
    diffuse.resize(nni+4);
    for (size_t i = 0; i < diffuse.size(); ++i ) {
	diffuse[i].resize(nnj+4);
    }

    // Apply the "standard filter".
    //for ( ipass = 0; ipass < npass; ++ipass ) {
	//SIMPLE_FILTER(gas->rho)
	//SIMPLE_FILTER(gas->e[0])
	//SIMPLE_FILTER(vel.x)
        //SIMPLE_FILTER(vel.y)
        //for (isp = 0; isp < nsp; ++isp) {
        //    SIMPLE_FILTER(gas->massf[isp])
	FCT_DIFFUSION(gas->rho)
	FCT_DIFFUSION(gas->e[0])
	FCT_DIFFUSION(vel.x)
        FCT_DIFFUSION(vel.y)
        for (isp = 0; isp < nsp; ++isp) {
            FCT_DIFFUSION(gas->massf[isp])
        }
	//}

    // We should make the thermodynamic state consistent, at least.
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    cell = get_cell(i,j);
	    Gas_data *gas= cell->fs->gas;
	    gm->eval_thermo_state_rhoe(*gas);
	    if ( get_viscous_flag() ) gm->eval_transport_coefficients(*gas);
	    if ( get_diffusion_flag() ) gm->eval_diffusion_coefficients(*gas);
	} // j loop
    } // i loop
    return SUCCESS;
} // end of apply_spatial_filter_diffusion()

int Block::apply_spatial_filter_anti_diffusion(double mu, int npass, int dimensions)
/// \brief Filter the cell-centred primary variables.
///
/// This filtering is done on a block-by-block basis.
/// Valid flow data are needed in the ghost cells at the edges.
/// \param alpha : filter coefficient
/// \param npass : this many passes
//
// To do: We should fix for 3D or remove. I think that there are no cases that need it. 
{
    FV_Cell *cell;
    double m2, m1, p1, p2;
    int isp;
    int i, j;
    Gas_model *gm = get_gas_model_ptr();
    int nsp = gm->get_number_of_species();
    vector<vector<double> > iadflux;
    vector<vector<double> > jadflux;
    iadflux.resize(nni+5);
    for (size_t i = 0; i < iadflux.size(); ++i ) {
	iadflux[i].resize(nnj+4);
    }
    jadflux.resize(nni+4);
    for (size_t i = 0; i < jadflux.size(); ++i ) {
	jadflux[i].resize(nnj+5);
    }
    // Apply the anti-diffusion.
    FCT_ANTI_DIFFUSION(gas->rho)
    FCT_ANTI_DIFFUSION(gas->e[0])
    FCT_ANTI_DIFFUSION(vel.x)
    FCT_ANTI_DIFFUSION(vel.y)
    for (isp = 0; isp < nsp; ++isp) {
	FCT_ANTI_DIFFUSION(gas->massf[isp])
    }
    
    // We should make the thermodynamic state consistent, at least.
    for ( i = imin; i <= imax; ++i ) {
	for ( j = jmin; j <= jmax; ++j ) {
	    cell = get_cell(i,j);
	    Gas_data *gas= cell->fs->gas;
	    gm->eval_thermo_state_rhoe(*gas);
	    if ( get_viscous_flag() ) gm->eval_transport_coefficients(*gas);
	    if ( get_diffusion_flag() ) gm->eval_diffusion_coefficients(*gas);
	} // j loop
    } // i loop
    return SUCCESS;
} // end of apply_spatial_filter_anti_diffusion()

//-----------------------------------------------------------------------------

int find_nearest_cell( double x, double y, double z, 
		       int *jb_near, int *i_near, int *j_near, int *k_near )
/// \brief Given an (x,y,z) position, locate the nearest cell centre.  
///
/// @param x, y, z : coordinates of the desired point
/// @param i_near, j_near, k_near : pointers to indices of the cell centre are stored here
/// @returns 1 for a close match, 0 if there were no close-enough cells.
{
    global_data *gd = get_global_data_ptr();
    Block *bdp;
    int ig, jg, kg, jbg;
    double dx, dy, dz, nearest, distance;

    jbg = 0;
    bdp = get_block_data_ptr(jbg);
    ig = bdp->imin; jg = bdp->jmin; kg = bdp->kmin;
    dx = x - bdp->get_cell(ig,jg,kg)->pos.x;
    dy = y - bdp->get_cell(ig,jg,kg)->pos.y;
    dz = z - bdp->get_cell(ig,jg,kg)->pos.z;
    nearest = sqrt(dx*dx + dy*dy + dz*dz);

    for ( int jb = 0; jb < gd->nblock; jb++ ) {
	bdp = get_block_data_ptr(jb);
	for ( int k = bdp->kmin; k <= bdp->kmax; ++k ) {
	    for ( int i = bdp->imin; i <= bdp->imax; ++i ) {
		for ( int j = bdp->jmin; j <= bdp->jmax; ++j ) {
		    dx = x - bdp->get_cell(i,j,k)->pos.x;
		    dy = y - bdp->get_cell(i,j,k)->pos.y;
		    dz = z - bdp->get_cell(i,j,k)->pos.z;
		    distance = sqrt(dx*dx + dy*dy + dz*dz);
		    if (distance < nearest) {
			nearest = distance; ig = i; jg = j; kg = k; jbg = jb;
		    }
		} // end j loop
	    } // end i loop
	} // for k
    } // for jb
    bdp = get_block_data_ptr(jbg);
    *jb_near = jbg; *i_near = ig; *j_near = jg; *k_near = kg;
    if ( nearest > bdp->get_ifi(ig,jg,kg)->length || 
	 nearest > bdp->get_ifj(ig,jg,kg)->length || 
	 (gd->dimensions == 3 && (nearest > bdp->get_ifk(ig,jg,kg)->length)) ) {
        printf("We really did not get close to (%e, %e, %e)\n", x, y, z);
        printf("Nearest %e lengths %e %e %e \n", nearest, bdp->get_ifi(ig,jg,kg)->length, bdp->get_ifj(ig,jg,kg)->length, bdp->get_ifk(ig,jg,kg)->length);
	return 0;
    } else {
        return 1;
    }
} // end find_nearest_cell()


// Some global data for locate_cell().
// To shorten the search, starting_block is the first block searched
// when requested locate_cell() is requested to locate the cell
// containing a specified point.
int starting_block = 0;

int locate_cell(double x, double y, double z,
	        int *jb_found, int *i_found, int *j_found, int *k_found)
// Returns 1 if a cell containing the sample point (x,y,z) is found, else 0.
// The indices of the containing cell are recorded, if found.
//
// To consider: maybe we should use *jb_found as the starting_block value.
{
    global_data *gd = get_global_data_ptr();
    Block *bdp;
    Vector3 p;
    int i, j, k, jb;
    
    p.x = x; p.y = y; p.z = z;
    *i_found = 0; *j_found = 0; *k_found = 0; *jb_found = 0;
    
    // Search the blocks, 
    // starting from the block in which the last point was found.
    for (jb = starting_block; jb < gd->nblock; jb++) {
	bdp = get_block_data_ptr(jb);
        for (k = bdp->kmin; k <= bdp->kmax; ++k) {
	    for (i = bdp->imin; i <= bdp->imax; ++i) {
		for (j = bdp->jmin; j <= bdp->jmax; ++j) {
		    if ( bdp->get_cell(i,j,k)->point_is_inside(p, gd->dimensions) ) {
			*i_found = i; *j_found = j; *k_found = k; 
			*jb_found = jb;
			starting_block = jb; // remember for next time
			return 1;
		    }
		} // j-loop
	    } // i-loop
	} // k-loop
    } // jb-loop
    
    // If we reach this point, then the point may be in one of the other blocks.
    for (jb = 0; jb < starting_block; jb++) {
	bdp = get_block_data_ptr(jb);
        for (k = bdp->kmin; k <= bdp->kmax; ++k) {
	    for (i = bdp->imin; i <= bdp->imax; ++i) {
		for (j = bdp->jmin; j <= bdp->jmax; ++j) {
		    if ( bdp->get_cell(i,j,k)->point_is_inside(p, gd->dimensions) ) {
			*i_found = i; *j_found = j; *k_found = k; 
			*jb_found = jb;
			starting_block = jb; // remember for next time
			return 1;
		    }
		} // j-loop
	    } // i-loop
	} // k-loop
    } // jb-loop
    
    // if we arrive here, we have not located the containing cell
    return 0;
} // end locate_cell()
