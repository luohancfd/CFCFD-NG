/// \file block_bgk.cxx
/// \ingroup eilmer3
/// \brief Daryl Bond's BGK-specific functions that work on the block data.
///
/// \version 23-Mar-2013 extracted from block.cxx
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

/// \brief Read the BGK discrete velocity distribution values
/// Returns a status flag.
int Block::read_BGK(std::string filename, double *sim_time,
		    size_t dimensions, bool zip_file)
{
#   define NCHAR 4000
    char line[NCHAR];
    char *gets_result;
    unsigned int i, j, k;
    size_t retries = 10;
    FILE *fp = NULL;
    gzFile zfp = NULL;
    if (id == 0) printf("read_BGK(): Start block %d.\n", static_cast<int>(id)); 
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
    sscanf(line, "%u %u %u", &i, &j, &k);
    if ( i != nni || j != nnj || k != ((dimensions == 3) ? nnk : 1) ) {
	printf("read_BGK(): block %d, mismatch in cell numbers\n", static_cast<int>(id));
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
    for (size_t k = kmin; k <= kmax; ++k ) {
	for (size_t j = jmin; j <= jmax; ++j ) {
	    for (size_t i = imin; i <= imax; ++i ) {
		get_cell(i,j,k)->fs->BGK_equilibrium();
	    }
	}
    }
    return SUCCESS;
} // end of Block::initialise_BGK_equilibrium()


/// \brief Write the BGK discrete velocity distribution values for a single block
///
/// This is "almost-Tecplot" POINT format.
int Block::write_BGK( std::string filename, double sim_time, size_t dimensions, bool zip_file )
{
    FILE *fp;
    gzFile zfp;
    string str;
    if (id == 0) {
	printf("write_BGK(): At t = %e, start block = %d.\n", sim_time, static_cast<int>(id));
    }
    if (zip_file) {
	fp = NULL;
	filename += ".gz";
	if ((zfp = gzopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_BGK(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	gzprintf(zfp, "%20.16e\n", sim_time);
	gzprintf(zfp,"\"pos.x\" \"pos.y\" \"pos.z\" \"volume\"");
	for (size_t ii = 0; ii < get_velocity_buckets(); ++ii) {
	    gzprintf(zfp," \"G[%d]\" \"H[%d]\"", static_cast<int>(ii), static_cast<int>(ii));
	}
	gzprintf(zfp, "\n%d %d %d\n", static_cast<int>(nni), static_cast<int>(nnj), static_cast<int>(nnk));
    } else {
	zfp = NULL;
	if ((fp = fopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_BGK(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	fprintf(fp, "%20.16e\n", sim_time);
	fprintf(fp,"\"pos.x\" \"pos.y\" \"pos.z\" \"volume\"");
	
	for (size_t ii = 0; ii < get_velocity_buckets(); ++ii) {
	    fprintf(fp," \"G[%d]\" \"H[%d]\"", static_cast<int>(ii), static_cast<int>(ii));
	}
	fprintf(fp, "\n%d %d %d\n", static_cast<int>(nni), static_cast<int>(nnj), static_cast<int>(nnk));
    }
    for ( size_t k = kmin; k <= kmax; ++k ) {
	for ( size_t j = jmin; j <= jmax; ++j ) {
	    for ( size_t i = imin; i <= imax; ++i ) {
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

