/// \file block_io.cxx
/// \ingroup eilmer3
/// \brief Functions to read and write whole block data.
///
/// \version 23-Mar-2013 extracted from block.cxx.
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

int Block::read_grid(std::string filename, size_t dimensions,
		     int zip_file, size_t time_level)
/// \brief Read the grid from a disc file as a set of cell vertices.
/// \returns 0 if successful but 1 if it hits the end of the grid file prematurely.
{
#   define NCHAR 132
    char line[NCHAR];
    char *gets_result;
    FV_Vertex *vp;
    unsigned int i, j, k;
    size_t retries = 10;
    FILE *fp = NULL;
    gzFile zfp = NULL;
    if (id == 0) printf("read_grid(): Start block %d.\n", static_cast<int>(id));
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
	printf("read_grid(): Empty grid file, block %d.\n", static_cast<int>(id));
	return BAD_INPUT_ERROR;
    }
    
    sscanf(line, "%u %u %u", &i, &j, &k);
    if (dimensions == 3) {
	if ( i != nni+1 || j != nnj+1 || k != nnk+1 ) {
	    printf("read_grid(): Mismatch in cell numbers, block %d\n",
		   static_cast<int>(id));
	    printf("    i=%u nni+1=%u j=%u nnj+1=%u k=%u nnk+1=%u\n", 
		   i, static_cast<unsigned int>(nni+1),
		   j, static_cast<unsigned int>(nnj+1),
		   k, static_cast<unsigned int>(nnk+1));
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
			printf("read_grid(): Premature end of file, block %d, vertex[%u,%u,%u]\n",
			       static_cast<int>(id), i, j, k);
			return BAD_INPUT_ERROR;
		    }
		    vp = get_vtx(i,j,k);
		    sscanf(line, "%lf %lf %lf", &(vp->pos[time_level].x),
			   &(vp->pos[time_level].y), &(vp->pos[time_level].z));
		} // for i
	    } // for j
	} // for k
    } else {
	// 2-dimensional case.
	if ( i != nni+1 || j != nnj+1 || k != 1 ) {
	    printf( "read_grid(): Mismatch in cell numbers, block %d\n",
		    static_cast<int>(id) );
	    printf( "    i=%u nni+1=%u j=%u nnj+1=%u k=%u nnk=%u\n", 
		    i, static_cast<unsigned int>(nni+1),
		    j, static_cast<unsigned int>(nnj+1),
		    k, static_cast<unsigned int>(nnk));
	    cout << "   more debug:" << endl;
	    cout << "   i= " << i << " nni+1= " << (nni+1) << " (i != nni+1)= " << (i != nni+1) << endl;
	    cout << "   j= " << j << " nnj+1= " << (nnj+1) << " (j != nnj+1)= " << (j != nnj+1) << endl;
	    cout << "   k= " << k << " nnk=" << (nnk) << " (k != 1)= " << (k != 1) << endl;
	    cout << "   ( i != nni+1 || j != nnj+1 || k != 1 )= " << ( i != nni+1 || j != nnj+1 || k != 1 ) << endl;
	    cout << "   WTF!" << endl;
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
		    printf("read_grid(): Premature end of file, block %d, vertex[%d,%d]\n",
			   static_cast<int>(id), static_cast<int>(i), static_cast<int>(j));
		    return BAD_INPUT_ERROR;
		}
		vp = get_vtx(i,j);
		sscanf(line, "%lf %lf", &(vp->pos[time_level].x), &(vp->pos[time_level].y));
		vp->pos[time_level].z = 0.0;
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


int Block::write_grid(std::string filename, double sim_time, size_t dimensions,
		      int zip_file, size_t time_level)
/// \brief Write the grid for a single block.
///
/// This is "almost-Tecplot" POINT format.
{
    FV_Vertex *vtx;
    FILE *fp;
    gzFile zfp;
    string str;
    size_t krangemax;
    if (id == 0) {
	// Just one block to announce that the writing has started.
	printf("write_grid(): At t = %e, start block = %d.\n", sim_time, static_cast<int>(id));
    }
    if (zip_file) {
	fp = NULL;
	filename += ".gz";
	if ((zfp = gzopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_grid(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	if ( dimensions == 2 ) {
	    gzprintf(zfp, "%d %d %d  # ni nj nk\n", static_cast<int>(nni+1), 
		     static_cast<int>(nnj+1), static_cast<int>(nnk));
	} else {
	    gzprintf(zfp, "%d %d %d  # ni nj nk\n", static_cast<int>(nni+1),
		     static_cast<int>(nnj+1), static_cast<int>(nnk+1));
	}
    } else {
	zfp = NULL;
	if ((fp = fopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_grid(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	if ( dimensions == 2 ) {
	    fprintf(fp, "%d %d %d  # ni nj nk\n", static_cast<int>(nni+1), 
		    static_cast<int>(nnj+1), static_cast<int>(nnk));
	} else {
	    fprintf(fp, "%d %d %d  # ni nj nk\n", static_cast<int>(nni+1),
		    static_cast<int>(nnj+1), static_cast<int>(nnk+1));
	}
    }
    if ( dimensions == 2 ) krangemax = kmax;
    else krangemax = kmax+1;
    for ( size_t k = kmin; k <= krangemax; ++k ) {
	for ( size_t j = jmin; j <= jmax+1; ++j ) {
	    for ( size_t i = imin; i <= imax+1; ++i ) {
		vtx = get_vtx(i,j,k);
		if (zip_file) {
		    gzprintf(zfp, "%20.12e %20.12e %20.12e\n", vtx->pos[time_level].x,
			     vtx->pos[time_level].y, vtx->pos[time_level].z);
		} else {
		    fprintf(fp, "%20.12e %20.12e %20.12e\n", vtx->pos[time_level].x,
			    vtx->pos[time_level].y, vtx->pos[time_level].z);
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
} // end of Block::write_grid()


/// \brief Read the flow solution (i.e. the flow data at cell centers) from a file.
/// Returns a status flag.
int Block::read_solution(std::string filename, double *sim_time, size_t dimensions,
			 int zip_file, size_t time_level)
{
#   define NCHAR 4000
    char line[NCHAR];
    char *gets_result;
    unsigned int i, j, k;
    size_t retries = 10;
    FILE *fp = NULL;
    gzFile zfp = NULL;
    if (id == 0) printf("read_solution(): Start block %d.\n", static_cast<int>(id)); 
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
    sscanf(line, "%u %u %u", &i, &j, &k);
    if ( i != nni || j != nnj || k != ((dimensions == 3) ? nnk : 1) ) {
	printf("read_solution(): block %d, mismatch in cell numbers\n", static_cast<int>(id));
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


int Block::write_solution(std::string filename, double sim_time, size_t dimensions,
			  int zip_file, size_t time_level)
/// \brief Write the flow solution (i.e. the primary variables at the
///        cell centers) for a single block.
///
/// This is "almost-Tecplot" POINT format.
{
    FILE *fp;
    gzFile zfp;
    string str;
    if (id == 0) {
	printf("write_solution(): At t = %e, start block = %d.\n",
	       sim_time, static_cast<int>(id));
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
	gzprintf(zfp, "%d %d %d\n", static_cast<int>(nni), static_cast<int>(nnj),
		 static_cast<int>(nnk));
    } else {
	zfp = NULL;
	if ((fp = fopen(filename.c_str(), "w")) == NULL) {
	    cerr << "write_solution(): Could not open " << filename << "; BAILING OUT" << endl;
	    exit( FILE_ERROR );
	}
	fprintf(fp, "%20.12e\n", sim_time);
	fprintf(fp, "%s\n", variable_list_for_cell().c_str());
	fprintf(fp, "%d %d %d\n", static_cast<int>(nni), static_cast<int>(nnj),
		static_cast<int>(nnk));
    }
    for ( size_t k = kmin; k <= kmax; ++k ) {
	for ( size_t j = jmin; j <= jmax; ++j ) {
	    for ( size_t i = imin; i <= imax; ++i ) {
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


int Block::write_history(std::string filename, double sim_time, int write_header, size_t time_level)
/// \brief Write out the flow solution in a (small) subset of cells.
///
/// This us usually done at a different (often smaller) time interval 
/// to the full flow solution.
/// Note that, after writing the header, the file is opened in append mode 
/// so that the data may accumulate.
{
    size_t i, j, k;
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
	for ( size_t ih = 0; ih < hncell; ++ih ) {
	    i = hicell[ih] + imin;
	    j = hjcell[ih] + jmin;
	    k = hkcell[ih] + kmin;
	    fprintf( fp, "%e %d %d %d ", sim_time, static_cast<int>(hicell[ih]),
		     static_cast<int>(hjcell[ih]), static_cast<int>(hkcell[ih]));
	    str = get_cell(i,j,k)->write_values_to_string();
	    fputs( str.c_str(), fp );
	    fputc( '\n', fp );
	}
    }
    fclose(fp);
    return SUCCESS;
} // end of Block::write_history()


/** \brief Computes the (pressure and shear) forces applied by the gas 
 *         to the bounding surface.
 *
 * We make use of geometric quantities stored at 
 * the cell interfaces.  
 * (Area is area per unit radian for axisymmetric calculations.)
 */
void Block::compute_x_forces(char *text_string, int ibndy, size_t dimensions, size_t time_level)
{
    double fx_p, fx_v, x1, y1, cosX, cosY, area;
    double xc, yc, d, vt, mu;
    size_t i, j, ivisc;
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
	    area = IFace->area[time_level];
	    mu   = IFace->fs->gas->mu;
	    fx_p += cell->fs->gas->p * area * cosX;
	    if ( ivisc ) {
		/* pieces needed to reconstruct the local velocity gradient */
		x1 = cell->vtx[0]->pos[time_level].x; y1 = cell->vtx[0]->pos[time_level].y;
		xc = cell->pos[time_level].x;   yc = cell->pos[time_level].y;
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
	    area = IFace->area[time_level];
	    mu   = IFace->fs->gas->mu;
	    fx_p -= cell->fs->gas->p * area * cosX;
	    if ( ivisc ) {
		/* pieces needed to reconstruct the local velocity gradient */
		x1 = cell->vtx[0]->pos[time_level].x; y1 = cell->vtx[0]->pos[time_level].y;
		xc = cell->pos[time_level].x;   yc = cell->pos[time_level].y;
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
	    area = IFace->area[time_level];
	    mu   = IFace->fs->gas->mu;
	    fx_p += cell->fs->gas->p * area * cosX;
	    if ( ivisc ) {
		/* pieces needed to reconstruct the local velocity gradient */
		x1 = cell->vtx[1]->pos[time_level].x; y1 = cell->vtx[1]->pos[time_level].y;
		xc = cell->pos[time_level].x;   yc = cell->pos[time_level].y;
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
	    area = IFace->area[time_level];
	    mu   = IFace->fs->gas->mu;
	    fx_p -= cell->fs->gas->p * area * cosX;
	    if ( ivisc ) {
		/* pieces needed to reconstruct the local velocity gradient */
		x1 = cell->vtx[0]->pos[time_level].x; y1 = cell->vtx[0]->pos[time_level].y;
		xc = get_cell(i,j)->pos[time_level].x;  yc = get_cell(i,j)->pos[time_level].y;
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
int Block::print_forces(FILE *fp, double t, size_t dimensions, size_t time_level)
{
    char msg_text[512], small_text[132];
    if ( bcp[NORTH]->xforce_flag == 1 ) {
	sprintf(small_text, "XFORCE: TIME %e BLOCK %d BNDY %d ", t, static_cast<int>(id), NORTH);
	strcpy(msg_text, small_text);
	this->compute_x_forces(small_text, NORTH, dimensions, time_level);
	strcat(msg_text, small_text);
	fprintf(fp, "%s\n",  msg_text);
    }
    if ( bcp[EAST]->xforce_flag == 1 ) {
	sprintf(small_text, "XFORCE: TIME %e BLOCK %d BNDY %d ", t, static_cast<int>(id), EAST);
	strcpy(msg_text, small_text);
	this->compute_x_forces(small_text, EAST, dimensions, time_level);
	strcat(msg_text, small_text);
	fprintf(fp, "%s\n",  msg_text);
    }
    if ( bcp[SOUTH]->xforce_flag == 1 ) {
	sprintf(small_text, "XFORCE: TIME %e BLOCK %d BNDY %d ", t, static_cast<int>(id), SOUTH);
	strcpy(msg_text, small_text);
	this->compute_x_forces( small_text, SOUTH, dimensions, time_level);
	strcat(msg_text, small_text);
	fprintf(fp, "%s\n",  msg_text);
    }
    if ( bcp[WEST]->xforce_flag == 1 ) {
	sprintf(small_text, "XFORCE: TIME %e BLOCK %d BNDY %d ", t, static_cast<int>(id), WEST);
	strcpy(msg_text, small_text);
	this->compute_x_forces(small_text, WEST, dimensions, time_level);
	strcat(msg_text, small_text);
	fprintf(fp, "%s\n",  msg_text);
    }
    return SUCCESS;
} // end print_forces()
